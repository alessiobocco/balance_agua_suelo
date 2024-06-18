
CrearSeriesHibridas <- function(input.value, climate.data, script) {
  
  # Extraer el identificador de la estacion meteorologica
  station_id <- input.value$weather_id
  
  # Identificar cultivo 
  crop_type <- base::unlist(base::strsplit(input.value$crop, '-')) %>% dplyr::first() 
  
  # Se genera la primera fecha de siembra (en formato fecha)
  fecha_inicio = base::as.POSIXct(glue::glue('{input.value$last_year}-{input.value$simulation_date}'))
  fecha_fin <- fecha_inicio + lubridate::days(365)
  
  # Se genera la primera fecha de siembra (en formato fecha)
  fecha_inflexion = base::as.POSIXct(glue::glue('{input.value$last_year}-{input.value$planting_date}'))
  
  series_hibridas <- purrr::map_dfr(
    .x = seq(input.value$first_year, input.value$last_year, 1),
    .f = function(year_inflexion) {
      
      
      serie_hibrida <- pr_crear_serie(id_estacion = station_id, 
                                      datos = climate.data,
                                      fecha_inicio = fecha_inicio, 
                                      fecha_inflexion = fecha_inflexion,
                                      fecha_fin = fecha_fin, 
                                      year_inflexion = year_inflexion)
      
      serie_hibrida$datos %>%
        dplyr::mutate(year_inflexion,
                      date_inflexion = serie_hibrida$fechas)
      
    }
  )
  
  return(series_hibridas)
  
}

EstimarFenologiaTrigo <- function(input.value, climate.data, script, config) {
  
  # Estacion meteoreologica sobre la que iterar
  station_id <- input.value$weather_id
  
  # Filtrar datos meteorologicos
  weather <- climate.data %>%
    dplyr::filter(station_id == !!station_id) %>%
    dplyr::select(-date) %>%
    dplyr::rename(date = date_inflexion)
  
  # Identificar cultivo 
  crop_type <- base::unlist(base::strsplit(input.value$crop, '-')) %>% dplyr::first() 
  # Identificador de los modelos de cultivos
  dssat_models = config$dssat$models
  # Cargar datos de los cultivos
  cultivos.archivo <- glue::glue('{config$dir$shared}/{dssat_models[[crop_type]]}.CUL')
  ecotipos.archivo <- glue::glue('{config$dir$shared}/{dssat_models[[crop_type]]}.ECO')
  # Identificar material
  material <- input.value$cultivar_id
  # Leer ceficientes geneticos y seleccionar el genotipo deseado
  coeficientes_cultivo <- DSSAT::read_cul(cultivos.archivo) %>%
    dplyr::filter(`VAR#` == !!material)
  coeficientes_ecotipo <- DSSAT::read_eco(ecotipos.archivo) %>%
    dplyr::filter(`ECO#` == coeficientes_cultivo$`ECO#`)
  
  # Mensajes 
  script$info(glue::glue("Estimando fenología para la estación {station_id},
                         cultivo: {crop_type} y genotipo {material}"))
  
  fenologia <- purrr::map_dfr(
    .x = unique(weather$year_inflexion),
    .f = function(year_inflexion) {
      
      # Parametros 
      params <- list(
        weather = weather,
        sowing_date = paste0(year_inflexion, "-", input.value$planting_date), # Sowing date in YYYY-MM-DD
        latitude = -33.0, # Latitude of the site
        # SDEPTH = 3.0, # Sowing depth in cm (commented out)
        GDDE = 6.2, # Growing degree days per cm seed depth required for emergence, GDD/cm
        # DSGFT = 200, # GDD from End Ear Growth to Start Grain Filling period (commented out)
        VREQ  = 400, # Vernalization required for max.development rate (VDays)
        PHINT = coeficientes_cultivo$PHINT, # Phyllochron. Estimate for PHINT is 95 degree days
        P1V = coeficientes_cultivo$P1V, # Development genetic coefficients, vernalization. 1 for spring type, 5 for winter type
        P1D = coeficientes_cultivo$P1D/100, # Development genetic coefficients, Photoperiod (1 - 6, low - high sensitive to day length)
        P5 = coeficientes_cultivo$P5,
        P1 = coeficientes_ecotipo$P1,
        # P5 = 500, # Grain filling degree days e.g., 500 degree-days (commented out)
        P6 = 250 # Approximate the thermal time from physiological maturity to harvest (commented out)
      )
      
      fenologia_trigo <- determine_phenology_stages_wheat(initparams = params)
      
      fenologia_trigo %<>%
        dplyr::mutate(station_id = station_id,
                      year = year_inflexion, sowing_date = input.value$planting_date,
                      material = coeficientes_cultivo$`VAR-NAME`) %>%
        dplyr::select(station_id, year, sowing_date, stage = istage_old, 
                      material, date, DOY, DAP, SUMDTT)
      
      return(fenologia_trigo)
    }
  )
  
  # Devolver resultados
  return(fenologia)
  
}

BalanceAguaSuelo <- function(input.value, climate.data, parametros.cultivo, fenologia.serie, initial.conditions, soil_file, script, config) {
  
  station_id = input.value$weather_id
  
  script$info(glue::glue("Calculando balance de agua en el suelo para estacion {station_id}
                         y cultivo {input.value$crop}"))
  
  # Filtrar datos climaticos 
  climate.data.i <- climate.data %>%
    dplyr::filter(station_id == !!station_id) %>%
    dplyr::mutate(`Year-DOY` = paste0(lubridate::year(date_inflexion), "-", lubridate::yday(date_inflexion))) %>%
    dplyr::select(`Year-DOY`, date, Srad = rad, Tmin = tmin, Tmax = tmax, Rain = prcp) %>%
    dplyr::group_by(date) %>%
    dplyr::mutate(ETref = evapotranspiracion.hargreaves(
      tmax = Tmax, tmin = Tmin, lat = -32, date = date, Rs = Srad)) %>%
    dplyr::ungroup()
  
  if (input.value$crop == "WH") {
    # Equivalencias entre fenologia y el modelo de FAO
    initial <- c("Sowing", "Term Spklt")
    development <- c("Term Spklt", "End Ear Gr")
    mid <- c("End Ear Gr", "End Gr Fil")
    end <- c("End Gr Fil", "Harvest")
  }
  
  # Caracteristicas del suelo
  soil.properties <- agua_disponible_perfil(soil_file = soil_file, soil_id = input.value$soil_id)
  
  initial.conditions %<>%
    dplyr::filter(soil_id == input.value$soil_id) %>%
    tidyr::unnest_longer(col = c("initial_water", "SLB", "SH2O")) 
  
  # Estimacion del balance
  balance <- purrr::map_dfr(
    .x = seq(input.value$first_year, input.value$last_year, 1),
    .f = function(year) {
      
      script$info(glue::glue("Calculando balance de agua en el suelo para estacion {station_id},
                         cultivo {input.value$crop} y año {year}"))
      
      # Extraccion de datos de la corrida i
      fenologia.serie.i <- fenologia.serie %>%
        dplyr::filter(station_id == input.value$weather_id, sowing_date == input.value$planting_date & year == !!year)
      
      # Calculo de longitud de cada etapa en dias calendarios
      Lini = fenologia.serie.i %>%
        dplyr::filter(stage %in% initial) %>%
        dplyr::pull(date) %>% diff() %>% as.numeric()
      
      Ldev = fenologia.serie.i %>%
        dplyr::filter(stage %in% development) %>%
        dplyr::pull(date) %>% diff() %>% as.numeric()
      
      Lmid <- fenologia.serie.i %>%
        dplyr::filter(stage %in% mid) %>%
        dplyr::pull(date) %>% diff() %>% as.numeric()
      
      Lend <- fenologia.serie.i %>%
        dplyr::filter(stage %in% end) %>%
        dplyr::pull(date) %>% diff() %>% as.numeric() + 15
      
      # Altura del cultivo
      hmax <- parametros.cultivo %>%
        dplyr::filter(cultivo == input.value$crop) %>%
        dplyr::pull(hmax)
      
      # Coeficientes del cultivo
      Kcbini = parametros.cultivo %>%
        dplyr::filter(cultivo == input.value$crop) %>%
        dplyr::pull(Kcbini)
      
      Kcbmid = parametros.cultivo %>%
        dplyr::filter(cultivo == input.value$crop) %>%
        dplyr::pull(Kcbmid)
      
      Kcbend = parametros.cultivo %>%
        dplyr::filter(cultivo == input.value$crop) %>%
        dplyr::pull(Kcbend)
      
      # Definir fecha de comienzo de la simulacion
      startDate <- fenologia.serie.i %>%
        dplyr::slice(1) %>%
        dplyr::mutate(`Year-DOY` = paste0(lubridate::year(date), "-", DOY)) %>%
        dplyr::pull(`Year-DOY`)
      
      endDate <- fenologia.serie.i %>%
        dplyr::slice(nrow(.)) %>%
        dplyr::mutate(`Year-DOY` = paste0(lubridate::year(date), "-", DOY+30)) %>%
        dplyr::pull(`Year-DOY`)
      
      balance.anio.i <- purrr::map_dfr(
        .x = unlist(config$params$initial_conditions),
        .f = function(water.condition) {
          
          theta0 <- initial.conditions %>% 
            dplyr::filter(soil_id == input.value$soil_id, initial_water == !!water.condition) %>% 
            dplyr::arrange(SLB) %>%
            dplyr::mutate(horizon_depth = diff(c(0, SLB))) %>%
            dplyr::summarise(water_content = weighted.mean(SH2O, horizon_depth))
          
          thetaFC <- soil.properties %>%
            dplyr::summarise(water_content = weighted.mean(SDUL, horizon_depth))
          
          thetaWP <- soil.properties %>%
            dplyr::summarise(water_content = weighted.mean(SLLL, horizon_depth))
          
          
          # Preparar objeto con parametros
          par <- list(Kcbini = Kcbini,  # Kcb Initial (FAO-56 Table 17)
                      Kcbmid = Kcbmid,  # Kcb Mid (FAO-56 Table 17)
                      Kcbend = Kcbend,  # Kcb End (FAO-56 Table 17)
                      Lini = Lini,        # Length Stage Initial (days) (FAO-56 Table 11)
                      Ldev = Ldev,        # Length Stage Development (days) (FAO-56 Table 11)
                      Lmid = Lmid,        # Length Stage Mid (days) (FAO-56 Table 11)
                      Lend = Lend,        # Length Stage End (days) (FAO-56 Table 11)
                      hini = 0.0500,    # Plant Height Initial (m)
                      hmax = hmax,    # Plant Height Maximum (m) (FAO-56 Table 12)
                      thetaFC = thetaFC, # Vol. Soil Water Content, Field Capacity (cm3/cm3)
                      thetaWP = thetaWP, # Vol. Soil Water Content, Wilting Point (cm3/cm3)
                      theta0 = theta0,  # Vol. Soil Water Content, Initial (cm3/cm3)
                      Zrini = 0.6000,   # Rooting Depth Initial (m)
                      Zrmax = 1.7000,   # Rooting Depth Maximum (m) (FAO-56 Table 22)
                      pbase = 0.6500,   # Depletion Fraction (p) (FAO-56 Table 22)
                      Ze = 0.1143,      # Depth of surface evaporation layer (m) (FAO-56 Table 19 and Page 144)
                      REW = 9.0000)     # Total depth Stage 1 evaporation (mm) (FAO-56 Table 19)
          
          balance <- fao56(par = par, wth = climate.data.i, startDate = startDate, endDate = endDate)
          
          balance %<>%
            dplyr::mutate(station = input.value$weather_id, year = year,
                          initial_water = water.condition, crop = input.value$crop,
                          sowing_date = input.value$planting_date) %>%
            dplyr::select(station, crop, year, sowing_date, initial_water, Date, ETref, 
                          Kcb, Kcmax, Ks, Kr, Ke, Kc, ETc, TAW, RAW, Kcadj, ETcadj, Rain, E, `T`) %>%
            dplyr::mutate(Date = as.Date(Date))
          
          # Devolver resultados
          return(balance)  
        }
      )
      # Devolver resultados
      return(balance.anio.i)
    }
  )
  
}
