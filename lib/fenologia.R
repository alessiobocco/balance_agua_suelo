
# ---------------------------------------------------------------------------- #
# --- Funciones propias ----
# ---------------------------------------------------------------------------- #

# Funciones para el calculo de Fenologia de trigo
# Calculo de la temperatura de la corona de trigo
crown_temperatures_wheat <- function(snow_depth = 0,
                               tmin = NA,
                               tmax = NA) {
  # Crown temperatures are simulated according to the original routines in CERES-Wheat and the correspond
  # to air temperatures for non-freezing temperatures. The minimum and maximum crown temperatures (Tcmin and Tcmax)
  # are calculated according to the maximum and minimun air temperatures (tmax and tmin), respectively.
  #
  # Parameters:
  #     snow_depth (int): Snow depth in centimeters (cm). Default value is set to zero.
  #     tmin (float): Minimum Temperature (°C)
  #     tmax (float): Maximum Temperature (°C)
  #
  # Returns:
  #     Tcmin (float): Minimum Crown Temperature (°C)
  #     Tcmax (float): Maximum Crown Temperature (°C)
  #     Tcrown (float): Optimum Crown Temperature (°C)
  
  
  if (is.na(tmin) | is.na(tmax)) {
    print("Please check out your inputs.")
    base::stop()
  }
  
  Tcmax = NA
  Tcmin = NA
  
  calc_CrownTemp <- function(Temp, snow_depth = 0) {
    Tcrown = Temp
    snow_depth = min(snow_depth, 15)
    if (Temp < 0) {
      Tcrown = 2.0 + Temp * (0.4 + 0.0018 * (snow_depth - 15) ** 2)
    }
    return(Tcrown)
  }
  
  
  # Crown temperature for maximum development rate
  Tcmax = calc_CrownTemp(tmax, snow_depth)
  # Crown temperature when snow is present and tmin < 0.
  Tcmin = calc_CrownTemp(tmin, snow_depth)
  
  Tcrown = (Tcmax + Tcmin) / 2
  
  return(c(Tcmax, Tcmin, Tcrown))
  
}

# Calculo del tiempo termico en trigo
thermal_time_calculation <- function(snow_depth = 0,
                                     tmin = NA,
                                     tmax = NA,
                                     Tbase = 0,
                                     Topt = 26,
                                     Ttop = 34) {
  # The daily thermal time (daily_TT) or Growing degree days calculation
  #
  # It's calculated from the daily average of maximum and minimum crown temperatures,
  # and is adjusted by genetic and environments factors.
  #
  # Parameters:
  #   m (str): Name of the model. Default is 'CERES'. Options: CERES, NWHEAT, WHAPS
  # snow_depth (int): Snow depth in centimeters (cm). Default value is set to zero.
  # tmin (float): Minimum Temperature (°C)
  # tmax (float): Maximum Temperature (°C)
  # Tbase (float): Base temperature for development from ecotype database. Default 0°C
  # Topt (float): Optimum temperature for development from species database. Default 26°C
  # Ttop (float): Maximum temperature for development from species database. Default 34°C
  #
  # Returns:
  #   dTT (float): Thermal time or Growing degree days
  
  if (is.na(tmin) || is.na(tmax)) {
    cat("Check input parameters\n")
    return(NA)
  }
  
  # Calculate Crown Temperatures
  temperature_crown <-
    crown_temperatures_wheat(snow_depth = snow_depth,
                       tmin = tmin,
                       tmax = tmax)
  Tcmax <- temperature_crown[1]
  Tcmin <- temperature_crown[2]
  Tcrown <- temperature_crown[3]
  tcdif = Tcmax - Tcmin
  dTT = Tcrown - Tbase
  
  if (Tcmax <= Tbase) {
    return(0.0)
  }
  
  tcdif <- Tcmax - Tcmin
  if (tcdif == 0)
    tcdif <- 1.0
  
  if (Tcmax < Topt) {
    if (Tcmin < Tbase) {
      tcor <- (Tcmax - Tbase) / tcdif
      dTT <- (Tcmax - Tbase) / 2 * tcor
    } else {
      dTT <- Tcrown - Tbase
    }
  } else if (Tcmax >= Topt && Tcmax < Ttop) {
    if (Tcmin < Topt) {
      tcor <- (Tcmax - Topt) / tcdif
      dTT <-
        (Topt - Tbase) / 2 * (1 + tcor) + Tcmin / 2 * (1 - tcor)
    } else {
      dTT <- Topt - Tbase
    }
  } else {
    if (Tcmin < Topt) {
      tcor <- (Tcmax - Ttop) / tcdif
      dTT_temp <- (Topt + Ttop - Tcmax) * tcor + Topt * (1 - tcor)
      tcor <- (Topt - Tcmin) / tcdif
      dTT <- dTT_temp * (1 - tcor) + (Tcmin + Topt) / 2 * tcor
    } else {
      tcor <- (Tcmax - Ttop) / tcdif
      dTT <- (Topt + Ttop - Tcmax) * tcor + Topt * (1 - tcor)
    }
  }
  
  return(round(dTT, 2))
}

# Calculo de la longitud del dia
day_length <- function(DOY = 1, lat = 0.0, p = 0.0) {
  # Convertir latitud a radianes
  S1 <- sin(lat * pi / 180)
  C1 <- cos(lat * pi / 180)
  
  # Calcular la declinación solar
  DEC <- 0.4093 * sin(0.0172 * (DOY - 82.2))
  
  # Calcular la longitud del día
  DLV <- (-S1 * sin(DEC) - 0.1047) / (C1 * cos(DEC))
  DLV <- max(DLV, -0.87)
  TWILEN <- 7.639 * acos(DLV)
  
  return(TWILEN)
}

# Funcion para el calculo del efecto del fotoperiodo
photoperiod_factor <- function(P1D = 3.675, day_length = 20) {
  # Photoperiod factor
  # Phenology is affected by photoperiod between emergence and floral initiation, and
  # thermal time is affected by a photoperiod factor.
  # Parameters:
  #     P1D (float): The sensitive to photoperiod (P1D) which is cultivar-specific. (1 - 6, low- high sensitive to day length)
  #     day_length (float): Day length in hours
  # Returns:
  #     DF (float): A photoperiod factor
  
  DF = 1 - (0.002 * P1D) * (20 - day_length)**2
  return(DF)
}

# Funcion para el calculo de la acumulacion de grados dia para la vernalizacion
vernalization <- function(Tcrown, tmin, tmax, cumvd = 0) {
  
  #         Calculate damage to crop due to cold weather.
  # 
  #         Vernalization is a response to relatively cold temperatures
  #         in some species that must occur before reproductive growth will begin.
  #         For wheat, temperature above zero to about 8°C seem to be the most effective
  #         for vernalization (Ahrens & Loomis, 1963; Chujo, 1966).
  # 
  #         Vernalization affects phenology between emergence and floral initiation.
  #         Spring-type winter cereals have little sensitivity to vernalization, which is
  #         the principal difference between them and the winter types.
  # 
  #         In the model, if the number of vernalization days (cumvd) is less than 10 and
  #         the maximum temperature exceeds 30°C, then the number of vernalization days decreases
  #         by 0.5 days per degree above 30°C. If cumvd is greater than 10, no devernalization is calculated.
  # 
  #         Vernalization is simulated from daily average crown temperature (Tcrown), daily maximum (tmax) and
  #         minimum (tmin) temperatures using the original CEREES approach.
  # 
  #         Parameters:
  #             Tcrown (float): daily average crown temperature (°C)
  #             tmin (float): daily average minimum temperature (°C)
  #             tmax (float): daily average maximum temperature (°C)
  #             cumvd (float): the number of vernalization days of total vernalization
  # 
  #         Returns:
  #             dV (float): Vernalization
  
  # Vernalization
  if (tmin < 15 && tmax > 0.0) {
    vd1 <- 1.4 - 0.0778 * Tcrown
    vd2 =  0.5 + 13.44 / (tmax - tmin + 3)**2 * Tcrown # Extract by CERES Wheat 2.0 fortran code
    vd <- min(1.0, vd1, vd2)
    vd <- max(vd, 0.0)
    cumvd <- cumvd + vd
  }
  # Devernalization
  else if (tmax > 30 && cumvd < 10) {
    cumvd <- cumvd - 0.5 * (tmax - 30)
    cumvd <- max(cumvd, 0.0)
  }
  
  return(cumvd)
  
}

# Funcion para el calculo del efecto del fotoperiodo
vernalization_factor <- function(P1V = 1.00,
                                 dV = 50,
                                 ISTAGE = 1) {
  # Calculation of vernalization factor.
  #
  # Phenology is affected by vernalization between emergence and floral initiation, and
  # thermal time is affected by a vernalization factor.
  #
  # Parameters:
  #   P1V (float): The sensitive to vernalization (P1V) which is cultivar-specific. 1 for spring type, 5 for winter type
  # dV (float): The total vernalization.
  #
  # Returns:
  #   VF (float): A vernalization factor
  # Set genetic coefficients to appropriate units
  #VSEN = params['P1V'] * 0.0054545 + 0.0003
  # VF = 1 - VSEN * (50 - CUMVD)
  if (ISTAGE == 1) {
    #or ISTAGE==2
    VF = 1 - (0.0054545 * P1V  + 0.0003) * (50 - dV)
    VF = max(min(VF, 1.0), 0.0)
  } else {
    VF = 1.0
  }
  return(VF)
}

# Función para determinar las etapas fenológicas del trigo
# 
# Esta función calcula las fechas de las diferentes etapas fenológicas del trigo,
# utilizando datos climáticos, fecha de siembra, y otros parámetros importantes.
#
# Argumentos:
# initparams: Lista opcional de parámetros de inicialización para sobrescribir los valores predeterminados.
# useDefault: Lógico, si es TRUE se utilizan los parámetros predeterminados.
#
# Retorno:
# Un data frame con las fechas y detalles de las etapas fenológicas del trigo.

determine_phenology_stages_wheat <- function(initparams = NULL) {
  # Parámetros predeterminados
  params <- list(
    weather = NULL,  # Datos climáticos del sitio
    sowing_date = "",  # Fecha de siembra en formato YYYY-MM-DD
    latitude = -90.0,  # Latitud del sitio
    longitude = -180.0,  # Longitud del sitio
    genotype = "",  # Nombre del genotipo en la base de datos de pedigrees de IWIN
    TT_TBASE = 0.0,  # Temperatura base, 2.0 para estimar HI
    TT_TEMPERATURE_OPTIMUM = 26,  # Temperatura óptima de tiempo térmico
    TT_TEMPERATURE_MAXIMUM = 34,  # Temperatura máxima de tiempo térmico
    CIVIL_TWILIGHT = 0.0,  # Ángulo del sol con el horizonte, ej. 6.0 para crepúsculo civil
    HI = 0.0,  # Índice de dureza
    SNOW = 0,  # Caída de nieve
    SDEPTH = 3.0,  # Profundidad de siembra en cm
    GDDE = 6.2,  # Grados-día de crecimiento por cm de profundidad de semilla requeridos para la emergencia
    DSGFT = 200,  # Grados-día desde el fin del crecimiento del espigado hasta el inicio del llenado de grano
    VREQ  = 505.0,  # Vernalización requerida para la tasa de desarrollo máximo (días vernalización)
    PHINT = 95.0,  # Estimación de la filocrono
    DAYS_GERMIMATION_LIMIT = 40,  # Días límite para la germinación
    TT_EMERGENCE_LIMIT = 300,  # Límite de tiempo térmico para la emergencia
    TT_TDU_LIMIT = 500,  # Límite de unidades de desarrollo térmico
    ADAH = 7  # Días de antesis después del encañado
  )
  
  # Sobrescribir parámetros predeterminados con los parámetros iniciales si se proporcionan
  if (!is.null(initparams)) {
    params <- modifyList(params, initparams)
  }
  
  # Validaciones de los parámetros
  if (is.null(params$sowing_date) || params$latitude == -90.0 || is.null(params$weather)) {
    return(NULL)
  }
  
  # Limites de grados-día
  P2 <- params$PHINT * 3
  P3 <- params$PHINT * 2
  P4 <- params$DSGFT
  P5 <- 430 + params$P5 * 20
  
  # Inicialización de la lista de etapas de crecimiento
  growstages <- list(
    `7` = list(istage_old = 'Sowing', istage = 'Fallow', desc = 'No crop present to Sowing', date = '', DOY = '', AGE = '', DAP = '', SUMDTT = '', status = ''),
    `8` = list(istage_old = 'Germinate', istage = 'Sowing', desc = 'Sowing to Germination', date = '', DOY = '', AGE = '', DAP = '', SUMDTT = '', status = ''),
    `9` = list(istage_old = 'Emergence', istage = 'Germinate', desc = 'Emergence to End of Juvenile', date = '', DOY = '', AGE = '', DAP = '', SUMDTT = '', status = ''),
    `1` = list(istage_old = 'Term Spklt', istage = 'Emergence', desc = 'Emergence to End of Juvenile', date = '', DOY = '', AGE = '', DAP = '', SUMDTT = '', status = ''),
    `2` = list(istage_old = 'End Veg', istage = 'End Juveni', desc = 'End of Juvenile to End of Vegetative growth', date = '', DOY = '', AGE = '', DAP = '', SUMDTT = '', status = ''),
    `2.5` = list(istage_old = 'Anthesis', istage = 'Anthesis', desc = 'Anthesis', date = '', DOY = '', AGE = '', DAP = '', SUMDTT = '', status = ''),
    `3` = list(istage_old = 'End Ear Gr', istage = 'End Veg', desc = 'End of Vegetative Growth to End of Ear Grow', date = '', DOY = '', AGE = '', DAP = '', SUMDTT = '', status = ''),
    `4` = list(istage_old = 'Beg Gr Fil', istage = 'End Ear Gr', desc = 'End of Ear Growth to Start of Grain Filling', date = '', DOY = '', AGE = '', DAP = '', SUMDTT = '', status = ''),
    `5` = list(istage_old = 'End Gr Fil', istage = 'Beg Gr Fil', desc = 'Start of Grain Filling to Maturity', date = '', DOY = '', AGE = '', DAP = '', SUMDTT = '', status = ''),
    `6` = list(istage_old = 'Harvest', istage = 'Maturity', desc = 'End Gr Fil', date = '', DOY = '', AGE = '', DAP = '', SUMDTT = '', status = '')
  )
  
  # DETERMINAR fecha de SIEMBRA ----
  ISTAGE <- 7
  SOWING_DATE <- as.Date(params$sowing_date)
  DOY <- as.integer(format(SOWING_DATE, "%j"))
  growstages[[as.character(ISTAGE)]]$date <- as.character(SOWING_DATE)
  growstages[[as.character(ISTAGE)]]$DOY <- DOY
  growstages[[as.character(ISTAGE)]]$AGE <- 0
  growstages[[as.character(ISTAGE)]]$SUMDTT <- 0
  growstages[[as.character(ISTAGE)]]$DAP <- 0
  growstages[[as.character(ISTAGE)]]$status <- 1
  
  # DETERMINAR fecha de GERMINACIÓN ----
  ISTAGE <- 8
  SUMDTT <- 0.0
  DAP <- 0
  ndays <- 1 # La germinación de la semilla es un proceso rápido y se asume que ocurre en un día
  
  GERMINATION_DATE <- ""
  w <- subset(params$weather, date == (SOWING_DATE + ndays))
  tmin <- as.numeric(w$tmin[ndays])
  tmax <- as.numeric(w$tmax[ndays])
  
  # Cálculo del tiempo térmico
  DTT <- thermal_time_calculation(snow_depth = params$SNOW, tmin = tmin, tmax = tmax, 
                                  Tbase = params$TT_TBASE, Topt = params$TT_TEMPERATURE_OPTIMUM, 
                                  Ttop = params$TT_TEMPERATURE_MAXIMUM)
  SUMDTT <- SUMDTT + DTT
  GERMINATION_DATE <- as.Date(w$date[ndays])
  CROP_AGE <- as.integer(GERMINATION_DATE - SOWING_DATE)
  DAP <- DAP + CROP_AGE
  
  growstages[[as.character(ISTAGE)]]$date <- as.character(GERMINATION_DATE)
  growstages[[as.character(ISTAGE)]]$DOY <- as.integer(format(GERMINATION_DATE, "%j"))
  growstages[[as.character(ISTAGE)]]$AGE <- CROP_AGE
  growstages[[as.character(ISTAGE)]]$SUMDTT <- round(SUMDTT, 1)
  growstages[[as.character(ISTAGE)]]$DAP <- DAP
  growstages[[as.character(ISTAGE)]]$status <- 1
  
  # DETERMINAR fecha de EMERGENCIA ----
  ISTAGE <- 9
  P9 <- 40 + params[['GDDE']] * params[['SDEPTH']] # Valores por defecto
  SUMDTT = 0.0
  w <- subset(params$weather, as.Date(date) >= as.Date(GERMINATION_DATE))
  EMERGENCE_DATE <- NA
  
  for (i in seq_len(nrow(w))) {
    tmin <- as.numeric(w$tmin[i])
    tmax <- as.numeric(w$tmax[i])
    DTT <- thermal_time_calculation(snow_depth = params[['SNOW']], tmin = tmin, tmax = tmax, 
                                    Tbase = params[['TT_TBASE']], Topt = params[['TT_TEMPERATURE_OPTIMUM']], 
                                    Ttop = params[['TT_TEMPERATURE_MAXIMUM']])
    SUMDTT <- SUMDTT + DTT
    
    if (SUMDTT >= P9 || SUMDTT > params[['TT_EMERGENCE_LIMIT']]) {
      EMERGENCE_DATE <- as.Date(w$date[i])
      CROP_AGE <- as.numeric(EMERGENCE_DATE - GERMINATION_DATE)
      DAP <- DAP + CROP_AGE
      growstages[[as.character(ISTAGE)]]$date <- as.character(EMERGENCE_DATE)
      growstages[[as.character(ISTAGE)]]$DOY <- as.integer(format(EMERGENCE_DATE, "%j"))
      growstages[[as.character(ISTAGE)]]$AGE <- CROP_AGE
      growstages[[as.character(ISTAGE)]]$SUMDTT <- round(SUMDTT, 1)
      growstages[[as.character(ISTAGE)]]$DAP <- DAP
      break # Salimos del bucle una vez que encontramos la fecha de emergencia
    }
  }
  
  # DETERMINAR fecha de FIN DE JUVENIL ----
  ISTAGE <- 1
  isVernalization <- TRUE
  SUMDTT = SUMDTT - P9 
  CUMVD = 0
  TDU = 0
  DF = 0.001
  shoot_lag = 40 # Asumido ser alrededor de 40 °C d
  shoot_rate = 1.5 # 1.5 °C d por mm. Derivado de estudios donde se midió el tiempo térmico a la emergencia y donde se conocía la profundidad de siembra
  sowing_depth = params[['SDEPTH']] * 10.0 # mm o 3 cm como en CERES
  T_emer = shoot_lag + shoot_rate * sowing_depth
  TT_emergence = min(T_emer, P9)
  
  w <- subset(params$weather, as.Date(date) >= as.Date(EMERGENCE_DATE))
  END_JUVENILE_DATE <- NA
  
  for (i in seq_len(nrow(w))) {
    tmin <- as.numeric(w$tmin[i])
    tmax <- as.numeric(w$tmax[i])
    DTT <- thermal_time_calculation(snow_depth = params[['SNOW']], tmin = tmin, tmax = tmax, 
                                    Tbase = params[['TT_TBASE']], Topt = params[['TT_TEMPERATURE_OPTIMUM']], 
                                    Ttop = params[['TT_TEMPERATURE_MAXIMUM']])
    
    if (isVernalization) {
      crown_temperature <- crown_temperatures_wheat(snow_depth = params[['SNOW']], tmin = tmin, tmax = tmax)
      Tcmax = crown_temperature[1]
      Tcmin = crown_temperature[2]
      Tcrown = crown_temperature[3]
      CUMVD <- vernalization(Tcrown = Tcrown, tmin = tmin, tmax = tmax, cumvd = CUMVD)
      
      if (CUMVD < params[['VREQ']]) {
        VF <- vernalization_factor(P1V = params[['P1V']], dV = CUMVD, ISTAGE = ISTAGE)
        
        if (VF < 0.3) {
          VF = max(VF, 0)
          TDU <- TDU + DTT * min(VF, DF)
        } else {
          DOY <- as.integer(format(as.Date(w$date[i]), "%j"))
          TWILEN <- day_length(DOY = DOY, lat = params[['latitude']], p = params$CIVIL_TWILIGHT)
          DF <- photoperiod_factor(P1D = params[['P1D']], day_length = TWILEN)
          DF <- max(DF, 0.001)
          
          TDU <- TDU + DTT * min(VF, DF)
        }
        SUMDTT <- TDU 
      } else {
        isVernalization <- FALSE
      }
    } else {
      SUMDTT <- SUMDTT + DTT
    }
    
    if(DF < 0) break
    
    if (SUMDTT > P9) {
      END_JUVENILE_DATE <- as.Date(w$date[i])
      CROP_AGE <- as.integer(END_JUVENILE_DATE - as.Date(EMERGENCE_DATE))
      DAP <- DAP + CROP_AGE
      growstages[[as.character(ISTAGE)]]$AGE <- CROP_AGE
      growstages[[as.character(ISTAGE)]]$date <- as.character(END_JUVENILE_DATE)
      growstages[[as.character(ISTAGE)]]$DOY <- as.integer(format(END_JUVENILE_DATE, "%j"))
      growstages[[as.character(ISTAGE)]]$SUMDTT <- round(SUMDTT, 1)
      growstages[[as.character(ISTAGE)]]$DAP <- DAP
      break
    }
  }
  
  # DETERMINAR fecha de FIN DE VEGETACIÓN ----
  ISTAGE <- 1 # Nota: debe continuar con 1 como en la etapa anterior
  isVernalization <- TRUE
  VF <- 1.0
  DF = 0.001
  w <- subset(params$weather, as.Date(date) >= as.Date(END_JUVENILE_DATE))
  END_VEGETATION_DATE <- NA
  
  if (nrow(w) > 0) {
    for (i in seq_len(nrow(w))) {
      tmin <- as.numeric(w$tmin[i])
      tmax <- as.numeric(w$tmax[i])
      crown_temperature <- crown_temperatures_wheat(snow_depth = params[['SNOW']], tmin = tmin, tmax = tmax)
      Tcmax = crown_temperature[1]
      Tcmin = crown_temperature[2]
      Tcrown = crown_temperature[3]
      DTT <- thermal_time_calculation(snow_depth = params[['SNOW']], tmin = tmin, tmax = tmax, 
                                      Tbase = params[['TT_TBASE']], Topt = params[['TT_TEMPERATURE_OPTIMUM']], 
                                      Ttop = params[['TT_TEMPERATURE_MAXIMUM']])
      if (isVernalization) {
        CUMVD <- vernalization(Tcrown, tmin, tmax, CUMVD)
        
        if (CUMVD < params[['VREQ']]) {
          VF <- vernalization_factor(P1V = params[['P1V']], dV = CUMVD, ISTAGE = ISTAGE)
          
          if (VF < 0.3) {
            TDU <- TDU + DTT * min(VF, DF)
          } else {
            DOY <- as.integer(format(as.Date(w$date[i]), "%j"))
            TWILEN <- day_length(DOY = DOY, lat = params[['latitude']])
            DF <- photoperiod_factor(P1D = params[['P1D']], day_length = TWILEN)
            TDU <- TDU + DTT * min(VF, DF)
          }
          SUMDTT <- TDU
        } else {
          isVernalization <- FALSE
        }
      } else {
        SUMDTT <- SUMDTT + DTT
      }
      
      if (DF < 0) {break}
      
      if (SUMDTT > (params[['TT_TDU_LIMIT']] * (params[['PHINT']] / 95.0))) {
        END_VEGETATION_DATE <- as.Date(w$date[i])
        CROP_AGE <- as.integer(END_VEGETATION_DATE - as.Date(END_JUVENILE_DATE))
        DAP <- DAP + CROP_AGE
        growstages[[as.character(ISTAGE)]]$AGE <- as.integer(END_VEGETATION_DATE - as.Date(EMERGENCE_DATE))
        growstages[[as.character(ISTAGE)]]$date <- as.character(END_VEGETATION_DATE)
        growstages[[as.character(ISTAGE)]]$DOY <- as.integer(format(END_VEGETATION_DATE, "%j"))
        growstages[[as.character(ISTAGE)]]$SUMDTT <- round(SUMDTT, 1)
        growstages[[as.character(ISTAGE)]]$DAP <- DAP
        break
      }
    }
  }
  
  # DETERMINAR FIN DE CRECIMIENTO DEL ESPIGADO ----
  ISTAGE <- 2
  P2 <- params[['PHINT']] * 3
  SUMDTT <- 0.0
  w <- subset(params$weather, as.Date(date) >= END_VEGETATION_DATE)
  END_OF_EAR_GROWTH_DATE <- NA
  
  if (nrow(w) > 0) {
    for (i in seq_len(nrow(w))) {
      tmin <- as.numeric(w$tmin[i])
      tmax <- as.numeric(w$tmax[i])
      DTT <- thermal_time_calculation(snow_depth = params[['SNOW']], tmin = tmin, tmax = tmax,
                                      Tbase = params[['TT_TBASE']], Topt = params[['TT_TEMPERATURE_OPTIMUM']], 
                                      Ttop = params[['TT_TEMPERATURE_MAXIMUM']])
      SUMDTT <- SUMDTT + DTT
      
      if (SUMDTT >= P2) {
        END_OF_EAR_GROWTH_DATE <- as.Date(w$date[i])
        CROP_AGE <- as.integer(END_OF_EAR_GROWTH_DATE - END_VEGETATION_DATE)
        DAP <- DAP + CROP_AGE
        growstages[[as.character(ISTAGE)]]$AGE <- CROP_AGE
        growstages[[as.character(ISTAGE)]]$date <- as.character(END_OF_EAR_GROWTH_DATE)
        growstages[[as.character(ISTAGE)]]$DOY <- as.integer(format(END_OF_EAR_GROWTH_DATE, "%j"))
        growstages[[as.character(ISTAGE)]]$SUMDTT <- round(SUMDTT, 1)
        growstages[[as.character(ISTAGE)]]$DAP <- DAP
        break
      }
    }
  }
  
  # DETERMINAR ANTESIS ----
  ISTAGE <- 2.5
  ADAH <- params[['ADAH']]
  CROP_AGE <- ADAH
  ANTHESIS_DATE <- as.Date(END_OF_EAR_GROWTH_DATE) + ADAH
  growstages[[as.character(ISTAGE)]]$AGE <- CROP_AGE
  growstages[[as.character(ISTAGE)]]$date <- format(ANTHESIS_DATE, "%Y-%m-%d")
  growstages[[as.character(ISTAGE)]]$DOY <- as.integer(format(ANTHESIS_DATE, "%j"))
  growstages[[as.character(ISTAGE)]]$SUMDTT <- round(SUMDTT, 1) # Esto podría necesitar un recálculo.
  growstages[[as.character(ISTAGE)]]$DAP <- DAP + ADAH
  
  # DETERMINAR FIN DE CRECIMIENTO DEL PENACHO ----
  ISTAGE <- 3
  P3 <- params[['PHINT']] * 2
  SUMDTT <- 0.0
  w <- subset(params$weather, as.Date(date) > END_OF_EAR_GROWTH_DATE)
  END_OF_PANNICLE_GROWTH_DATE <- NA
  
  if (nrow(w) > 0) {
    for (i in seq_len(nrow(w))) {
      tmin <- as.numeric(w$tmin[i])
      tmax <- as.numeric(w$tmax[i])
      DTT <- thermal_time_calculation(snow_depth = params[['SNOW']], tmin = tmin, tmax = tmax,
                                      Tbase = params[['TT_TBASE']], Topt = params[['TT_TEMPERATURE_OPTIMUM']], 
                                      Ttop = params[['TT_TEMPERATURE_MAXIMUM']])
      SUMDTT <- SUMDTT + DTT
      
      if (SUMDTT >= P3) {
        END_OF_PANNICLE_GROWTH_DATE <- as.Date(w$date[i])
        CROP_AGE <- as.integer(END_OF_PANNICLE_GROWTH_DATE - END_OF_EAR_GROWTH_DATE)
        DAP <- DAP + CROP_AGE
        growstages[[as.character(ISTAGE)]]$AGE <- CROP_AGE
        growstages[[as.character(ISTAGE)]]$date <- as.character(END_OF_PANNICLE_GROWTH_DATE)
        growstages[[as.character(ISTAGE)]]$DOY <- as.integer(format(END_OF_PANNICLE_GROWTH_DATE, "%j"))
        growstages[[as.character(ISTAGE)]]$SUMDTT <- round(SUMDTT, 1)
        growstages[[as.character(ISTAGE)]]$DAP <- DAP
        break
      }
    }
  }
  
  # DETERMINAR INICIO DEL LLENADO DEL GRANO ----
  ISTAGE <- 4
  P4 <- params[['DSGFT']]
  SUMDTT <- 0.0
  w <- subset(params$weather, as.Date(date) >= END_OF_PANNICLE_GROWTH_DATE)
  BEGIN_GRAIN_FILLING_DATE <- NA
  
  if (nrow(w) > 0) {
    for (i in seq_len(nrow(w))) {
      tmin <- as.numeric(w$tmin[i])
      tmax <- as.numeric(w$tmax[i])
      DTT <- thermal_time_calculation(snow_depth = params[['SNOW']], tmin = tmin, tmax = tmax,
                                      Tbase = params[['TT_TBASE']], Topt = params[['TT_TEMPERATURE_OPTIMUM']], 
                                      Ttop = params[['TT_TEMPERATURE_MAXIMUM']])
      SUMDTT <- SUMDTT + DTT
      
      if (SUMDTT >= P4) {
        BEGIN_GRAIN_FILLING_DATE <- as.Date(w$date[i])
        CROP_AGE <- as.integer(BEGIN_GRAIN_FILLING_DATE - END_OF_PANNICLE_GROWTH_DATE)
        DAP <- DAP + CROP_AGE
        growstages[[as.character(ISTAGE)]]$AGE <- CROP_AGE
        growstages[[as.character(ISTAGE)]]$date <- as.character(BEGIN_GRAIN_FILLING_DATE)
        growstages[[as.character(ISTAGE)]]$DOY <- as.integer(format(BEGIN_GRAIN_FILLING_DATE, "%j"))
        growstages[[as.character(ISTAGE)]]$SUMDTT <- round(SUMDTT, 1)
        growstages[[as.character(ISTAGE)]]$DAP <- DAP
        break
      }
    }
  }
  
  # DETERMINAR FIN DEL LLENADO DEL GRANO ----
  ISTAGE <- 5
  P5 <- params[['P5']]
  SUMDTT <- 0
  w <- subset(params$weather, as.Date(date) >= BEGIN_GRAIN_FILLING_DATE)
  END_GRAIN_FILLING_DATE <- NA
  
  if (nrow(w) > 0) {
    for (i in seq_len(nrow(w))) {
      tmin <- as.numeric(w$tmin[i])
      tmax <- as.numeric(w$tmax[i])
      DTT <- thermal_time_calculation(snow_depth = params[['SNOW']], tmin = tmin, tmax = tmax,
                                      Tbase = params[['TT_TBASE']], Topt = params[['TT_TEMPERATURE_OPTIMUM']], 
                                      Ttop = params[['TT_TEMPERATURE_MAXIMUM']])
      SUMDTT <- SUMDTT + DTT
      
      if (SUMDTT >= P5) {
        END_GRAIN_FILLING_DATE <- as.Date(w$date[i])
        CROP_AGE <- as.integer(END_GRAIN_FILLING_DATE - BEGIN_GRAIN_FILLING_DATE)
        DAP <- DAP + CROP_AGE
        growstages[[as.character(ISTAGE)]]$AGE <- CROP_AGE
        growstages[[as.character(ISTAGE)]]$date <- as.character(END_GRAIN_FILLING_DATE)
        growstages[[as.character(ISTAGE)]]$DOY <- as.integer(format(END_GRAIN_FILLING_DATE, "%j"))
        growstages[[as.character(ISTAGE)]]$SUMDTT <- round(SUMDTT, 1)
        growstages[[as.character(ISTAGE)]]$DAP <- DAP
        break
      }
    }
  }
  
  # DETERMINAR COSECHA ----
  ISTAGE <- 6
  SUMDTT <- 0.0
  estimateHarvest <- TRUE
  P6 <- params[['P6']]
  w <- subset(params$weather, as.Date(date) >= as.Date(END_GRAIN_FILLING_DATE))
  HARVEST <- NA
  
  if (nrow(w) > 0) {
    for (i in seq_len(nrow(w))) {
      tmin <- as.numeric(w$tmin[i])
      tmax <- as.numeric(w$tmax[i])
      DTT <- thermal_time_calculation(snow_depth = params[['SNOW']], tmin = tmin, tmax = tmax,
                                      Tbase = params[['TT_TBASE']], Topt = params[['TT_TEMPERATURE_OPTIMUM']], 
                                      Ttop = params[['TT_TEMPERATURE_MAXIMUM']])
      SUMDTT <- SUMDTT + DTT
      
      if (SUMDTT >= P6) {
        HARVEST <- as.Date(w$date[i])
        CROP_AGE <- as.integer(HARVEST - as.Date(END_GRAIN_FILLING_DATE))
        DAP <- DAP + CROP_AGE
        growstages[[as.character(ISTAGE)]]$AGE <- CROP_AGE
        growstages[[as.character(ISTAGE)]]$date <- as.character(HARVEST)
        growstages[[as.character(ISTAGE)]]$DOY <- as.integer(format(HARVEST, "%j"))
        growstages[[as.character(ISTAGE)]]$SUMDTT <- round(SUMDTT, 1)
        growstages[[as.character(ISTAGE)]]$DAP <- DAP
        growstages[[as.character(ISTAGE)]]$status <- 1
        break
      }
    }
  }
  
  # Transformar la lista en un data frame
  growstages <- do.call(rbind, lapply(growstages, function(x) data.frame(matrix(unlist(x), ncol = length(x), byrow = T), stringsAsFactors = FALSE)))
  
  # Dar nombres a las columnas del data frame
  names(growstages) <- c("istage_old", "istage", "desc", "date", "DOY", "AGE", "DAP", "SUMDTT", "status")
  
  # Convertir columnas a tipos de datos adecuados
  growstages$date <- as.Date(growstages$date)
  growstages$DOY <- as.numeric(growstages$DOY)
  growstages$AGE <- as.numeric(growstages$AGE)
  growstages$DAP <- as.numeric(growstages$DAP)
  growstages$SUMDTT <- as.numeric(growstages$SUMDTT)
  growstages$status <- as.factor(growstages$status)
  
  # Devolver resultados
  return(growstages)
}


# Funciones para el calculo de la Fenologia de maiz
crown_temperatures <- function(snow_depth = 0,
                               Tmin = NA,
                               Tmax = NA) {
  # Crown temperatures are simulated according to the original routines in CERES-Wheat and the correspond
  # to air temperatures for non-freezing temperatures. The minimum and maximum crown temperatures (Tcmin and Tcmax)
  # are calculated according to the maximum and minimun air temperatures (Tmax and Tmin), respectively.
  #
  # Parameters:
  #     snow_depth (int): Snow depth in centimeters (cm). Default value is set to zero.
  #     Tmin (float): Minimum Temperature (°C)
  #     Tmax (float): Maximum Temperature (°C)
  #
  # Returns:
  #     Tcmin (float): Minimum Crown Temperature (°C)
  #     Tcmax (float): Maximum Crown Temperature (°C)
  #     Tcrown (float): Optimum Crown Temperature (°C)
  
  
  if (is.na(Tmin) | is.na(Tmax)) {
    print("Please check out your inputs.")
    base::stop()
  }
  
  Tcmax = NA
  Tcmin = NA
  
  calc_CrownTemp <- function(Temp, snow_depth = 0) {
    Tcrown = Temp
    snow_depth = min(snow_depth, 15)
    if (Temp < 0) {
      Tcrown = 2.0 + Temp * (0.4 + 0.0018 * (snow_depth - 15) ** 2)
    }
    return(Tcrown)
  }
  
  
  # Crown temperature for maximum development rate
  Tcmax = calc_CrownTemp(Tmax, snow_depth)
  # Crown temperature when snow is present and TMIN < 0.
  Tcmin = calc_CrownTemp(Tmin, snow_depth)
  
  Tcrown = (Tcmax + Tcmin) / 2
  
  return(c(Tcmax, Tcmin, Tcrown))
  
}

# Función para calcular el tiempo térmico diario (DTT) basado en varios factores ambientales.
computeDTT <- function(TMAX, TMIN, TBASE, TOPT, ROPT, ISTAGE, LEAFNO, XS, SRAD, DAYL) {
  # ------------------------------------------------------------
  # Compute thermal time based on new method developed by J.T.R
  # at CYMMIT, 5/5/98.  TBASE, TOPT, and ROPT are read in 
  # from the species file.
  #------------------------------------------------------------
  
  # DOPT, Devlopment optimum temperature, is set to TOPT 
  #   during vegetative growth and to ROPT after anthesis
  
  
  # DOPT es la temperatura óptima de desarrollo, que cambia después de la antesis
  DOPT <- TOPT
  # Cambiar DOPT a ROPT después de la antesis (etapas 4 a 6)
  if (ISTAGE > 3 && ISTAGE <= 6) {
    DOPT <- ROPT
  }
  
  # Inicializar DTT a 0.0
  DTT <- 0.0
  
  # Si la temperatura máxima es menor que la base, DTT es 0
  if (TMAX < TBASE) {
    DTT <- 0.0
    # Si la temperatura mínima es mayor que DOPT, ajustar DTT
  } else if (TMIN > DOPT) {
    DTT <- DOPT - TBASE
  } 
  
  # Considerar condiciones especiales para plantas con menos de 10 hojas
  if (LEAFNO <= 10) {
    # Si hay nieve, ajustar DTT basado en la temperatura de la corona
    if (XS > 0.0) {
      TEMPCROWN <- crown_temperatures(snow_depth = XS, Tmin = TMIN, Tmax = TMAX)
      DTT <- (TEMPCROWN[2] + TEMPCROWN[1]) / 2.0 - TBASE
    } else {
      # Sin nieve, calcular la temperatura del suelo basada en la radiación solar y temperatura
      ACOEF <- 0.01061 * SRAD + 0.5902
      TDSOIL <- ACOEF * TMAX + (1 - ACOEF) * TMIN
      TNSOIL <- 0.36354 * TMAX + 0.63646 * TMIN
      
      # Ajustar DTT si la temperatura del suelo diurna es menor que la base
      if (TDSOIL < TBASE) {
        DTT <- 0.0
        
        return(DTT)
      } else {
        # Ajustar temperaturas del suelo para no exceder límites
        TNSOIL <- max(TNSOIL, TBASE)
        TDSOIL <- min(TDSOIL, DOPT)
        
        # Calcular temperatura media del suelo basada en el día
        TMSOIL <- TDSOIL * (DAYL / 24) + TNSOIL * ((24 - DAYL) / 24)
        
        # Ajustar DTT basado en la temperatura media del suelo
        if (TMSOIL < TBASE) {
          DTT <- (TBASE + TDSOIL) / 2 - TBASE
        } else {
          DTT <- (TNSOIL + TDSOIL) / 2 - TBASE
        }
        
        # Asegurar que DTT no exceda el óptimo
        DTT <- min(DTT, DOPT - TBASE)
        
        # Asegurar que DTT no sea negativo
        print(max(DTT, 0.0))
        return(max(DTT, 0.0))
      }
    }
  } else if (TMIN < TBASE || TMAX > DOPT) {
    # Ajustar DTT si las temperaturas están fuera del rango óptimo
    DTT <- 0.0
    # Calcular DTT a lo largo del día en incrementos de 1 hora
    for (I in 1:24) {
      TH <- (TMAX + TMIN) / 2 + (TMAX - TMIN) / 2 * sin(pi / 12 * I)
      TH <- max(min(TH, DOPT), TBASE)
      DTT <- DTT + (TH - TBASE) / 24
    }
  } else {
    # Calcular DTT para condiciones normales
    DTT <- (TMAX + TMIN) / 2 - TBASE
  }
  
  # Asegurar que DTT no sea negativo
  DTT <- max(DTT, 0.0)
  
  # Aquí podrías acumular DTT y CUMDTT si se manejan fuera de la función
  # SUMDTT <- SUMDTT + DTT
  # CUMDTT <- CUMDTT + DTT
  
  return(DTT)
}

# Determinacion de las etapas fenologicas de maiz
determine_phenology_stages_maize <- function(initparams = NULL) {
  # Default Parameters
  params <- list(
    weather = NULL,  # Weather data of the site
    sowing_date = "",  # Sowing date in YYYY-MM-DD
    latitude = -90.0,  # Latitude of the site
    longitude = -180.0,  # Longitude of the site
    genotype = "",  # Name of the grand parent in IWIN pedigrees database
    TT_TBASE = 0.0,  # Base Temperature, 2.0 to estimate HI
    TT_TEMPERATURE_OPTIMUM = 26,  # Thermal time optimum temperature
    TT_TEMPERATURE_MAXIMUM = 34,  # Thermal time maximum temperature
    CIVIL_TWILIGHT = 0.0,  # Sun angle with the horizon. e.g., p = 6.0 : civil twilight
    HI = 0.0,  # Hardiness Index
    SNOW = 0,  # Snow fall
    SDEPTH = 3.0,  # Sowing depth in cm
    GDDE = 6.2,  # Growing degree days per cm seed depth required for emergence, GDD/cm
    DSGFT = 200,  # GDD from End Ear Growth to Start Grain Filling period
    VREQ  = 505.0,  # Vernalization required for max. development rate (VDays)
    PHINT = 95.0,  # Phyllochron estimate
    P1V = 1.0,  # Vernalization coefficient
    P1D = 3.675,  # Photoperiod response coefficient
    P5 = 500,  # Grain filling degree days
    P6 = 250,  # Thermal time from maturity to harvest
    DAYS_GERMIMATION_LIMIT = 40,  # Days to germination threshold
    TT_EMERGENCE_LIMIT = 300,  # Thermal time to emergence threshold
    TT_TDU_LIMIT = 500,  # Thermal development units limit
    ADAH = 6,  # Anthesis days after heading
    bruteforce = FALSE,  # Use brute force algorithm or not
    brute_params = list(
      obsEmergenceDAP = NULL,  # Observed days after planting to emergence
      obsHeadingDAP = NULL,  # Observed days after planting to heading
      obsAnthesisDAP = NULL,  # Observed days after planting to anthesis
      obsMaturityDAP = NULL,  # Observed days after planting to maturity
      max_tries = 300,  # Maximum number of tries for brute force
      error_lim = 0.5,  # Error limit for brute force
      gdde_steps = 1.0,  # Step size for adjusting GDDE
      maxGDDE = 50,  # Maximum GDDE limit
      phint_steps = 1.0,  # Step size for adjusting PHINT
      maxPHINT = 150,  # Maximum PHINT limit
      adap_steps = 1,  # Step size for adjusting ADAH
      maxADAP = 10,  # Maximum ADAH limit
      p5_steps = 1,  # Step size for adjusting P5
      maxP5 = 2000  # Maximum P5 limit
    )
  )
  
  # Override default parameters with initial parameters if provided
  if (!is.null(initparams)) {
    params <- modifyList(params, initparams)
  }
  
  # Validation checks
  if ( is.null(params$sowing_date)) {
    cat("Sowing date not defined\n")
    return(NULL)
  }
  if (params$latitude == -90.0 || is.null(params$latitude)) {
    cat("Problem with location of the site. Check the geographic coordinates.\n")
    return(NULL)
  }
  if (is.null(params$weather)) {
    cat("Weather data is not available\n")
    return(NULL)
  }
  
  
  # ---------------------
  # GDD limits
  # ---------------------
  #P2 <- params$PHINT * 3
  #P3 <- params$PHINT * 2
  #P4 <- params$DSGFT # 200 # APSIM-Wheat = 120 # GDD from End Ear Growth to Start Grain Filling period
  #P5 <- params$P5 # 430 + params$P5 * 20 # DSSAT v4.8
  # P6 <- params$P5 # 250 (commented out)
  # P9 <- 40 + params$GDDE * params$SDEPTH (commented out)
  
  # Initializing growstages list
  growstages <- list(
    `7` = list(istage_old = 'Sowing', istage = 'Fallow', desc = 'No crop present to Sowing', date = '', DOY = '', AGE = '', DAP = '', SUMDTT = '', status = ''),
    `8` = list(istage_old = 'Germinate', istage = 'Sowing', desc = 'Sowing to Germination', date = '', DOY = '', AGE = '', DAP = '', SUMDTT = '', status = ''),
    `9` = list(istage_old = 'Emergence', istage = 'Germinate', desc = 'Emergence to End of Juvenile', date = '', DOY = '', AGE = '', DAP = '', SUMDTT = '', status = ''),
    `1` = list(istage_old = 'End. Juvenile', istage = 'Emergence', desc = 'Emergence to End of Juvenile', date = '', DOY = '', AGE = '', DAP = '', SUMDTT = '', status = ''),
    `2` = list(istage_old = 'Pannicle Initiation', istage = 'End Juvenile', desc = 'End of Juvenile to Pannicle initiation', date = '', DOY = '', AGE = '', DAP = '', SUMDTT = '', status = ''),
    `3` = list(istage_old = 'End leaf growth', istage = 'Pannicle Initiation', desc = 'End of Pannicle Initiation to End leaf growth', date = '', DOY = '', AGE = '', DAP = '', SUMDTT = '', status = ''),
    `4` = list(istage_old = 'End pannicle growth', istage = 'End leaf growth', desc = 'End of End leaf growth to End pannicle growth', date = '', DOY = '', AGE = '', DAP = '', SUMDTT = '', status = ''),
    `5` = list(istage_old = 'Grain fill', istage = 'End pannicle growth', desc = 'End of End pannicle growth to Grain fill', date = '', DOY = '', AGE = '', DAP = '', SUMDTT = '', status = ''),
    `6` = list(istage_old = 'Maturity', istage = 'Grain fill', desc = 'Start of Grain Filling to Maturity', date = '', DOY = '', AGE = '', DAP = '', SUMDTT = '', status = '')
  )
  
  # Custom error class equivalent in R
  StageFailed <- function(message, istage, err) {
    structure(list(message = message, istage = istage, err = err), class = 'StageFailed')
  }
  
  # S3 method for printing StageFailed
  print.StageFailed <- function(x) {
    cat(x$message, "Stage (", x$istage, ") - Error: ", x$err, "\n")
  }
  
  # -------------------------------------------------------------------------- #
  # DETERMINE SOWING DATE ----
  # -------------------------------------------------------------------------- #
  ISTAGE <- 7
  tryCatch({
    SOWING_DATE <- as.Date(params$sowing_date)
    DOY <- as.integer(format(SOWING_DATE, "%j"))
    CUMPH   = 0.514
    
    
    growstages[[as.character(ISTAGE)]]$date <- as.character(SOWING_DATE)
    growstages[[as.character(ISTAGE)]]$DOY <- DOY
    growstages[[as.character(ISTAGE)]]$AGE <- 0
    growstages[[as.character(ISTAGE)]]$SUMDTT <- 0
    growstages[[as.character(ISTAGE)]]$DAP <- 0
    growstages[[as.character(ISTAGE)]]$status <- 1
    # cat("Sowing date:", SOWING_DATE, "\n") # Uncomment for debugging
  }, error = function(e) {
    stop(StageFailed("Problem initializing the determination of phenological stage. Please check your input parameters such as sowing date or latitude of the site", ISTAGE, e))
  })
  
  # --------------------------------------------------------------------------
  
  # -------------------------------------------------------------------------- #
  # DETERMINE GERMINATION DATE ----
  # -------------------------------------------------------------------------- #
  ISTAGE <- 8
  tryCatch({
    SUMDTT <- 0.0
    DAP <- 0
    ndays <- 1 # Seed germination is a rapid process and is assumed to occur in one day
    
    # Assuming 'weather' is a data frame with 'DATE', 'TMIN', and 'TMAX' columns
    GERMINATION_DATE <- ""
    w <- subset(params$weather, date == (SOWING_DATE + ndays))
    Tmin <- as.numeric(w$tmin)
    Tmax <- as.numeric(w$tmax)
    Srad <- as.numeric(w$rad)
    # Duracion del dia
    DOY <- as.integer(format(as.Date(w$date), "%j"))
    TWILEN <- day_length(DOY = DOY, lat = params[['latitude']])
    
    # Thermal time calculation
    DTT <- computeDTT(TMIN = Tmin, TMAX = Tmax, TBASE = params$TT_TBASE, TOPT = params$TT_TEMPERATURE_OPTIMUM,
                      ROPT = params$TT_TEMPERATURE_MAXIMUM, ISTAGE = 8, 
                      LEAFNO = 0, XS = 0, 
                      SRAD = Srad, DAYL = TWILEN)
    
    SUMDTT <- SUMDTT + DTT
    GERMINATION_DATE <- as.Date(w$date[ndays])
    CROP_AGE <- as.integer(GERMINATION_DATE - SOWING_DATE)
    DAP <- DAP + CROP_AGE
    
    growstages[[as.character(ISTAGE)]]$date <- as.character(GERMINATION_DATE)
    growstages[[as.character(ISTAGE)]]$DOY <- as.integer(format(GERMINATION_DATE, "%j"))
    growstages[[as.character(ISTAGE)]]$AGE <- CROP_AGE
    growstages[[as.character(ISTAGE)]]$SUMDTT <- round(SUMDTT, 1)
    growstages[[as.character(ISTAGE)]]$DAP <- DAP
    growstages[[as.character(ISTAGE)]]$status <- 1
    
    # Uncomment for debugging: cat("Germination date:", GERMINATION_DATE, "\n")
  }, error = function(err) {
    x <- StageFailed("Problem determining germination date.", ISTAGE, err)
    print(x)
    return(NULL)
  })
  # --------------------------------------------------------------------------
  
  PLTPOP = 5.5
  ROWSPC = 43
  RUE = 4.2
  # ------------------------------------------------------------------------- #
  # DETERMINE SEEDLING EMERGENCE DATE ----
  # ------------------------------------------------------------------------- #
  ISTAGE <- 9
  P9 <- 40 + params[['GDDE']] * params[['SDEPTH']] # Valores por defecto
  LFWT    = 0.2 # Leaf weight at emergence, g/plant         
  LEAFNO = as.integer(1) # Leaf number at emergence, #/plant      
  PLA     = 1 # Leaf area at emergence, cm2/plant
  SENLA  = 0.0
  LAI    = PLTPOP*PLA*0.0001 
  CUMPH   = 0.514
  tryCatch({
    SUMDTT = 0.0
    #print("Growing degree days from germination to emergence (P9): ",P9) 
    # The crop will die if germination has not occurred before a certain period (eg. 40 days)
    w <- subset(params$weather, date > GERMINATION_DATE)
    EMERGENCE_DATE <- NA
    
    
    for (i in seq_len(nrow(w))) {
      Tmin <- as.numeric(w$tmin[i])
      Tmax <- as.numeric(w$tmax[i])
      Srad <- as.numeric(w$rad[i])
      # Duracion del dia
      DOY <- as.integer(format(as.Date(w$date[i]), "%j"))
      TWILEN <- day_length(DOY = DOY, lat = params[['latitude']])
      # Asumimos que la función 'thermal_time_calculation' ya está definida y devuelve un valor numérico
      DTT <- computeDTT(TMIN = Tmin, TMAX = Tmax, 
                        TBASE = params[['TT_TBASE']], 
                        TOPT = params[['TT_TEMPERATURE_OPTIMUM']], 
                        ROPT = params[['TT_TEMPERATURE_MAXIMUM']],
                        ISTAGE = 9, LEAFNO = 0, XS = 0, 
                        SRAD = Srad, DAYL = TWILEN)
      SUMDTT <- SUMDTT + round(DTT, 2)
      
      #------------------------------------------------------------
      #                      STATE VARIABLES
      #------------------------------------------------------------
      
      PLA <- pmax(0, PLA)
      LAI <- pmax(0, LAI)
      PLTPOP <- pmax(0, PLTPOP)
      
      
      if (SUMDTT >= P9 || SUMDTT > params[['TT_EMERGENCE_LIMIT']]) {
        EMERGENCE_DATE <- as.Date(w$date[i]) 
        CROP_AGE <- as.numeric(EMERGENCE_DATE - GERMINATION_DATE)
        DAP <- DAP + CROP_AGE
        growstages[[as.character(ISTAGE)]]$date <- as.character(EMERGENCE_DATE)
        growstages[[as.character(ISTAGE)]]$DOY <- as.integer(format(EMERGENCE_DATE, "%j"))
        growstages[[as.character(ISTAGE)]]$AGE <- CROP_AGE
        growstages[[as.character(ISTAGE)]]$SUMDTT <- round(SUMDTT, 1)
        growstages[[as.character(ISTAGE)]]$DAP <- DAP
        break # Salimos del bucle una vez que encontramos la fecha de emergencia
      }
    }
    
  }, error = function(err) {
    x <- StageFailed("Problem determining emergence date.", ISTAGE, err)
    print(x)
    return(NULL)
  })
  # ------------------------------------------------------------------------- 
  
  # -------------------------------------------------------------------------------------- #
  # DETERMINE DURATION OF VEGETATIVE PAHSE (END JUVENILE DATE - END OF VEGETATION GROWTH ----
  # -------------------------------------------------------------------------------------- #
  ISTAGE <- 1
  
  tryCatch({
    SUMDTT = SUMDTT - P9 
    CUMVD = 0
    TLNO   = 30.0
    XN = 0
    XLFWT = 0
    #LAI = 0
    #CUMPH = 0
    
    w <- subset(params$weather, as.Date(date) >= as.Date(EMERGENCE_DATE))
    END_JUVENILE_DATE <- NA
    
    
    for (i in seq_len(nrow(w))) {
      
      # Tiempo termico
      Tmin <- as.numeric(w$tmin[i])
      Tmax <- as.numeric(w$tmax[i])
      Srad <- as.numeric(w$rad[i])
      # Duracion del dia
      DOY <- as.integer(format(as.Date(w$date[i]), "%j"))
      TWILEN <- day_length(DOY = DOY, lat = params[['latitude']])
      # Asumimos que la función 'thermal_time_calculation' ya está definida y devuelve un valor numérico
      DTT <- computeDTT(TMIN = Tmin, TMAX = Tmax, TBASE = params[['TT_TBASE']], 
                        TOPT = params[['TT_TEMPERATURE_OPTIMUM']], 
                        ROPT = params[['TT_TEMPERATURE_MAXIMUM']],
                        ISTAGE = 1, LEAFNO = LEAFNO, XS = 0, 
                        SRAD = Srad, DAYL = TWILEN)
      SUMDTT <- SUMDTT + DTT
      
      # Permitir que las primeras 5 hojas se expandan más rápido que el intervalo Phyllochron reduciendo PC
      PC <- 1.0
      if (CUMPH < 5.0) {
        PC <- 0.66 + 0.068 * CUMPH # Establece PC a 0.7 para la hoja 1 y 1.0 para la hoja 5
      }
      
      TI <- DTT / (params[["PHINT"]] * PC)
      CUMPH <- CUMPH + TI
      XN <- CUMPH + 1.0
      LEAFNO <- as.integer(XN)
      print(LEAFNO)
      
      # Crecimiento del área foliar, cm2/planta/día
      PLAG <- 3.0 * XN^2 * TI * 1 # No stress factors considered
      
      if (XN < 4.0) {
        # Permitir que las primeras 4 hojas se expandan un 30% más rápido
        PLAG <- 4.0 * XN * TI * 1 # No stress factors considered
      }
      
      PLA <- PLA + PLAG
      # XLFWT - nuevo peso de la hoja hoy, g/planta
      XLFWT <- max((PLA / 250)^1.25, LFWT)
      # GROLF - tasa de crecimiento de la hoja, g/planta/día
      GROLF <- XLFWT - LFWT
      
      LFWT <- LFWT + GROLF
      if (GROLF > 0) {
        # SLAN - senescencia normal de la hoja desde la emergencia, cm2/planta
        SLAN <- SUMDTT * PLA / 10000
      } else {
        SLAN = 0
      }
      # Ajuste para la respiración
      LFWT <- LFWT - SLAN / 600
      XLFWT <- XLFWT + LFWT
      #------------------------------------------------------------
      #               Compute Leaf Senescence Factors
      #------------------------------------------------------------
      
      # Senescence due to light competition
      SLFC <- 1.00        
      if (LAI > 4.0) {
        SLFC <- 1.0 - 0.008 * (LAI - 4.0)
      }
      
      # Senescence due to temperature
      SLFT <- 1.0
      if (Tmin <= 6.0) {
        SLFT <- max(0.0, 1.0 - 0.01 * (Tmin - 6.0)^2)
      }
      SLFT <- max(SLFT, 0.0)
      
      # Daily rate of leaf senescence
      PLAS <- (PLA - SENLA) * (1.0 - min(SLFC, SLFT))
      SENLA <- SENLA + PLAS
      SENLA <- max(SENLA, SLAN)
      SENLA <- min(SENLA, PLA)         
      LAI <- (PLA - SENLA) * PLTPOP * 0.0001
      print(LAI)
      
      #------------------------------------------------------------
      #                      STATE VARIABLES
      #------------------------------------------------------------
      
      LFWT <- pmax(0, LFWT)
      PLA <- pmax(0, PLA)
      LAI <- pmax(0, LAI)
      PLTPOP <- pmax(0, PLTPOP)
      
      XLAI <- LAI  # Leaf area index, m2/m2
      
      # Calculating canopy height (CANHT) based on an empirical equation
      #if (XLAI >= MAXLAI) {  # Keeps CANHT at the maximum value
      #  CANHT <- XLAI / (0.4238 * PLTPOP + 0.3424) * CANHT_POT
      #  if (CANHT > CANHT_POT) {
      #    CANHT <- CANHT_POT
      #  }
      #}
      #MAXLAI <- max(MAXLAI, XLAI)  # Maximum XLAI season
      
      if (SUMDTT > params[['P1']] ) {
        END_JUVENILE_DATE <- as.Date(w$date[i])
        CROP_AGE <- as.integer(END_JUVENILE_DATE - as.Date(EMERGENCE_DATE))
        DAP <- DAP + CROP_AGE
        growstages[[as.character(ISTAGE)]]$AGE <- CROP_AGE
        growstages[[as.character(ISTAGE)]]$date <- as.character(END_JUVENILE_DATE)
        growstages[[as.character(ISTAGE)]]$DOY <- as.integer(format(END_JUVENILE_DATE, "%j"))
        growstages[[as.character(ISTAGE)]]$SUMDTT <- round(SUMDTT, 1)
        growstages[[as.character(ISTAGE)]]$DAP <- DAP
        break
      }
    }
    
  }, error = function(err) {
    x <- StageFailed("Problem determining end juvenile date.", ISTAGE, err)
    print(x)
    return(NULL)
  })
  
  # --------------------------------------------------------------------------------------
  
  # -------------------------------------------------------------------------- #
  # DETERMINE END VEGETATION DATE - ISTAGE = 2 - End of Juvenile Stage to Tassel Initiation ----
  # -------------------------------------------------------------------------- #
  ISTAGE <- 2 
  
  tryCatch({
    
    CUMVD = 0
    SIND = 0
    PDTT = SUMDTT - params[["P1"]]
    w <- subset(params$weather, as.Date(date) > as.Date(END_JUVENILE_DATE))
    END_VEGETATION_DATE <- NA
    
    if (nrow(w) > 0) {
      for (i in seq_len(nrow(w))) {
        Tmin <- as.numeric(w$tmin[i])
        Tmax <- as.numeric(w$tmax[i])
        Srad <- as.numeric(w$rad[i])      
        # Duracion del dia
        DOY <- as.integer(format(as.Date(w$date[i]), "%j"))
        TWILEN <- day_length(DOY = DOY, lat = params[['latitude']])
        # Asumimos que la función 'thermal_time_calculation' ya está definida y devuelve un valor numérico
        DTT <- computeDTT(TMIN = Tmin, TMAX = Tmax, TBASE = params[['TT_TBASE']], 
                          TOPT = params[['TT_TEMPERATURE_OPTIMUM']], 
                          ROPT = params[['TT_TEMPERATURE_MAXIMUM']],
                          ISTAGE = ISTAGE, LEAFNO = LEAFNO, XS = 0, 
                          SRAD = Srad, DAYL = TWILEN)
        SUMDTT <- SUMDTT + DTT
        
        # Radiacion Daily Photosynthesis Rate
        PAR = Srad*.50    # PAR local variable
        
        LIFAC  = 1.5 - 0.768 * ((ROWSPC * 0.01)**2 * PLTPOP)**0.1 
        PCO2  = 1.014286
        
        # JIL 08/01/2006 Intercepted PAR (MJ/plant d)
        IPAR = PAR/PLTPOP * (1.0 - exp(-LIFAC * LAI))
        CARBO = IPAR * RUE * PCO2
        CARBO = pmax(0, CARBO)
        TAVGD = 0.25*Tmin+0.75*Tmax
        
        # Permitir que las primeras 5 hojas se expandan más rápido que el intervalo Phyllochron reduciendo PC
        PC <- 1.0
        if (CUMPH < 5.0) {
          PC <- 0.66 + 0.068 * CUMPH # Establece PC a 0.7 para la hoja 1 y 1.0 para la hoja 5
        }
        
        TI <- DTT / (params[["PHINT"]] * PC)
        CUMPH <- CUMPH + TI
        XN <- CUMPH + 1.0
        LEAFNO <- as.integer(XN)
        print(LEAFNO)
        
        # Cálculo del área foliar creciente
        PLAG <- 3.5 * XN^2 * TI * 1
        PLA <- PLA + PLAG
        # Calcula el nuevo peso de la hoja basado en el área de la hoja y ajusta el crecimiento de la hoja.
        XLFWT <- (PLA / 267.0)^1.25
        GROLF <- XLFWT - LFWT
        
        if (GROLF >= CARBO * 0.75) {
          GROLF <- CARBO * 0.75
          PLA <- (LFWT + GROLF)^0.8 * 267.0
        }
        
        # Actualiza el crecimiento de la raíz y el peso de la hoja, y calcula la senescencia normal de la hoja.
        LFWT <- LFWT + GROLF
        SLAN <- SUMDTT * PLA / 10000.0
        LFWT <- LFWT - SLAN / 600.0
        
        if (TWILEN > params[["P20"]]) {
          RATEIN <- 1.0 / (params[["DJTI"]]  + params[["P2"]] * (TWILEN - params[["P20"]]))
        } else {
          RATEIN <- 1.0 / params[["DJTI"]] 
        }
        
        PDTT <- 1.0
        SIND <- SIND + RATEIN * PDTT
        
        #------------------------------------------------------------
        #               Compute Leaf Senescence Factors
        #------------------------------------------------------------
        
        # Senescence due to light competition
        SLFC <- 1.00        
        if (LAI > 4.0) {
          SLFC <- 1.0 - 0.008 * (LAI - 4.0)
        }
        
        # Senescence due to temperature
        SLFT <- 1.0
        if (Tmin <= 6.0) {
          SLFT <- max(0.0, 1.0 - 0.01 * (Tmin - 6.0)^2)
        }
        SLFT <- max(SLFT, 0.0)
        
        # Daily rate of leaf senescence
        PLAS <- (PLA - SENLA) * (1.0 - min(SLFC, SLFT))
        SENLA <- SENLA + PLAS
        SENLA <- max(SENLA, SLAN)
        SENLA <- min(SENLA, PLA)         
        LAI <- (PLA - SENLA) * PLTPOP * 0.0001
        print(LAI)
        
        #------------------------------------------------------------
        #                      STATE VARIABLES
        #------------------------------------------------------------
        
        LFWT <- pmax(0, LFWT)
        PLA <- pmax(0, PLA)
        LAI <- pmax(0, LAI)
        PLTPOP <- pmax(0, PLTPOP)
        
        XLAI <- LAI  # Leaf area index, m2/m2
        
        
        if (SIND >= 1) {
          END_VEGETATION_DATE <- as.Date(w$date[i])
          CROP_AGE <- as.integer(END_VEGETATION_DATE - as.Date(END_JUVENILE_DATE))
          DAP <- DAP + CROP_AGE
          CROP_AGE_2 <- as.integer(END_VEGETATION_DATE - as.Date(EMERGENCE_DATE))
          growstages[[as.character(ISTAGE)]]$AGE <- CROP_AGE_2
          growstages[[as.character(ISTAGE)]]$date <- as.character(END_VEGETATION_DATE)
          growstages[[as.character(ISTAGE)]]$DOY <- as.integer(format(END_VEGETATION_DATE, "%j"))
          growstages[[as.character(ISTAGE)]]$SUMDTT <- round(SUMDTT, 1)
          growstages[[as.character(ISTAGE)]]$DAP <- DAP
          break
        }
      }
    }
    
  }, error = function(err) {
    x <- StageFailed("Problem determining end of vegetation growth date.", ISTAGE, err)
    print(x)
    return(NULL)
  })
  
  # -------------------------------------------------------------------------- 
  
  # ---------------------------------------------------------------------------------------------- #
  # DETERMINE END OF EAR GROWTH - End of Leaf Growth to Pannicle Initiation (End leaf growth) ----
  #----------------------------------------------------------------------------------------------- #
  ISTAGE <- 3
  TLNO <- SUMDTT/(params[["PHINT"]]*0.5) + 5.0           
  P3 <- ((TLNO + 0.5) * params[["PHINT"]]) - SUMDTT 
  BSGDD = 250.0        # Beginning of ear growth (gdd)
  
  w <- subset(params$weather, as.Date(date) > END_VEGETATION_DATE)
  END_OF_PANNICLE_GROWTH_DATE <- NA
  
  if (nrow(w) > 0) {
    # Reiniciar contador de GDD
    SUMDTT_2 = 0
    for (i in seq_len(nrow(w))) {
      Tmin <- as.numeric(w$tmin[i])
      Tmax <- as.numeric(w$tmax[i])
      Srad <- as.numeric(w$rad[i])      
      # Duracion del dia
      DOY <- as.integer(format(as.Date(w$date[i]), "%j"))
      TWILEN <- day_length(DOY = DOY, lat = params[['latitude']])
      # Asumimos que la función 'thermal_time_calculation' ya está definida y devuelve un valor numérico
      DTT <- computeDTT(TMIN = Tmin, TMAX = Tmax, TBASE = params[['TT_TBASE']], 
                        TOPT = params[['TT_TEMPERATURE_OPTIMUM']], 
                        ROPT = params[['TT_TEMPERATURE_MAXIMUM']],
                        ISTAGE = ISTAGE, LEAFNO = LEAFNO, XS = 0, 
                        SRAD = Srad, DAYL = TWILEN)
      SUMDTT_2 <- SUMDTT_2 + DTT
      SUMDTT <- SUMDTT + DTT
      
      # Radiacion Daily Photosynthesis Rate
      PAR = Srad*.50    # PAR local variable
      
      LIFAC  = 1.5 - 0.768 * ((ROWSPC * 0.01)**2 * PLTPOP)**0.1 
      PCO2  = 1.014286
      
      # JIL 08/01/2006 Intercepted PAR (MJ/plant d)
      IPAR = PAR/PLTPOP * (1.0 - exp(-LIFAC * LAI))
      CARBO = IPAR * RUE * PCO2
      CARBO = pmax(0, CARBO)
      TAVGD = 0.25*Tmin+0.75*Tmax
      
      # Permitir que las primeras 5 hojas se expandan más rápido que el intervalo Phyllochron reduciendo PC
      #PC <- 1.0
      #if (CUMPH < 5.0) {
      #  PC <- 0.66 + 0.068 * CUMPH # Establece PC a 0.7 para la hoja 1 y 1.0 para la hoja 5
      #}
      
      #TI <- DTT / (params[["PHINT"]] * PC)
      #CUMPH <- CUMPH + TI
      #XN <- CUMPH + 1.0
      #LEAFNO <- as.integer(XN)
      #print(LEAFNO)
      
      # Cálculo del área foliar creciente
      PLAG <- 3.5 * XN^2 * TI * 1
      PLA <- PLA + PLAG
      # Calcula el nuevo peso de la hoja basado en el área de la hoja y ajusta el crecimiento de la hoja.
      XLFWT <- (PLA / 267.0)^1.25
      GROLF <- XLFWT - LFWT
      
      if (GROLF >= CARBO * 0.75) {
        GROLF <- CARBO * 0.75
        PLA <- (LFWT + GROLF)^0.8 * 267.0
      }
      
      # Actualiza el crecimiento de la raíz y el peso de la hoja, y calcula la senescencia normal de la hoja.
      LFWT <- LFWT + GROLF
      SLAN <- SUMDTT * PLA / 10000.0
      LFWT <- LFWT - SLAN / 600.0
      
      #             Do not allow increase in leaf number during 
      #               expansion of final leaf.
      if (SUMDTT_2 > P3 - (2. * params[["PHINT"]])) {
        TI     = DTT/(params[["PHINT"]]*PC)
        CUMPH = CUMPH - TI
      } else {
        TI     = DTT/(params[["PHINT"]]*PC)
        CUMPH  = CUMPH + TI
      }
      
      XN     = CUMPH + 1.0
      LEAFNO = as.integer(XN)
      
      print(LEAFNO)
      
      #------------------------------------------------------------
      #               Compute Leaf Senescence Factors
      #------------------------------------------------------------
      
      # Senescence due to light competition
      SLFC <- 1.00        
      if (LAI > 4.0) {
        SLFC <- 1.0 - 0.008 * (LAI - 4.0)
      }
      
      # Senescence due to temperature
      SLFT <- 1.0
      if (Tmin <= 6.0) {
        SLFT <- max(0.0, 1.0 - 0.01 * (Tmin - 6.0)^2)
      }
      SLFT <- max(SLFT, 0.0)
      
      # Daily rate of leaf senescence
      PLAS <- (PLA - SENLA) * (1.0 - min(SLFC, SLFT))
      SENLA <- SENLA + PLAS
      SENLA <- max(SENLA, SLAN)
      SENLA <- min(SENLA, PLA)         
      LAI <- (PLA - SENLA) * PLTPOP * 0.0001
      print(LAI)
      
      #------------------------------------------------------------
      #                      STATE VARIABLES
      #------------------------------------------------------------
      
      LFWT <- pmax(0, LFWT)
      PLA <- pmax(0, PLA)
      LAI <- pmax(0, LAI)
      PLTPOP <- pmax(0, PLTPOP)
      
      XLAI <- LAI  # Leaf area index, m2/m2
      
      
      if (SUMDTT_2 > P3) {
        END_OF_PANNICLE_GROWTH_DATE <- as.Date(w$date[i])
        CROP_AGE <- as.integer(END_OF_PANNICLE_GROWTH_DATE - END_VEGETATION_DATE)
        DAP <- DAP + CROP_AGE
        growstages[[as.character(ISTAGE)]]$AGE <- CROP_AGE
        growstages[[as.character(ISTAGE)]]$date <- as.character(END_OF_PANNICLE_GROWTH_DATE)
        growstages[[as.character(ISTAGE)]]$DOY <- as.integer(format(END_OF_PANNICLE_GROWTH_DATE, "%j"))
        growstages[[as.character(ISTAGE)]]$SUMDTT <- round(SUMDTT_2, 1)
        growstages[[as.character(ISTAGE)]]$DAP <- DAP
        break
      }
    }
  } else {
    cat("Error reading weather data for end of ear growth phase\n")
  }
  #-----------------------------------------------------------------------------------------------
  
  # ---------------------------------------------------------------------------------------------- #
  # DETERMINE END OF PANNICLE GROWTH - End of Leaf Growth to Beginning Effective Gra ----
  # ---------------------------------------------------------------------------------------------- #
  ISTAGE = 4
  SUMDTT = SUMDTT_2 - P3
  IDURP  = 0
  
  w <- subset(params$weather, as.Date(date) > END_OF_PANNICLE_GROWTH_DATE)
  END_OF_LEAF_GROWTH_DATE <- NA
  
  if (nrow(w) > 0) {
    for (i in seq_len(nrow(w))) {
      Tmin <- as.numeric(w$tmin[i])
      Tmax <- as.numeric(w$tmax[i])
      Srad <- as.numeric(w$rad[i])      
      # Duracion del dia
      DOY <- as.integer(format(as.Date(w$date[i]), "%j"))
      TWILEN <- day_length(DOY = DOY, lat = params[['latitude']])
      # Asumimos que la función 'thermal_time_calculation' ya está definida y devuelve un valor numérico
      DTT <- computeDTT(TMIN = Tmin, TMAX = Tmax, TBASE = params[['TT_TBASE']], 
                        TOPT = params[['TT_TEMPERATURE_OPTIMUM']], 
                        ROPT = params[['TT_TEMPERATURE_MAXIMUM']],
                        ISTAGE = ISTAGE, LEAFNO = LEAFNO, XS = 0, 
                        SRAD = Srad, DAYL = TWILEN)
      SUMDTT <- SUMDTT + DTT
      
      SLAN   = PLA*(0.05+SUMDTT/200.0*0.05)
      LFWT   = LFWT  - SLAN/600.0
      
      #------------------------------------------------------------
      #               Compute Leaf Senescence Factors
      #------------------------------------------------------------
      
      # Senescence due to light competition
      SLFC <- 1.00        
      if (LAI > 4.0) {
        SLFC <- 1.0 - 0.008 * (LAI - 4.0)
      }
      
      # Senescence due to temperature
      SLFT <- 1.0
      if (Tmin <= 6.0) {
        SLFT <- max(0.0, 1.0 - 0.01 * (Tmin - 6.0)^2)
      }
      SLFT <- max(SLFT, 0.0)
      
      # Daily rate of leaf senescence
      PLAS <- (PLA - SENLA) * (1.0 - min(SLFC, SLFT))
      SENLA <- SENLA + PLAS
      SENLA <- max(SENLA, SLAN)
      SENLA <- min(SENLA, PLA)         
      LAI <- (PLA - SENLA) * PLTPOP * 0.0001
      print(LAI)
      
      #------------------------------------------------------------
      #                      STATE VARIABLES
      #------------------------------------------------------------
      
      LFWT <- pmax(0, LFWT)
      PLA <- pmax(0, PLA)
      LAI <- pmax(0, LAI)
      PLTPOP <- pmax(0, PLTPOP)
      
      XLAI <- LAI  # Leaf area index, m2/m2
      
      if (SUMDTT >= params[["DSGFT"]]) {
        END_OF_LEAF_GROWTH_DATE <- as.Date(w$date[i])
        CROP_AGE <- as.integer(END_OF_PANNICLE_GROWTH_DATE - END_OF_LEAF_GROWTH_DATE)
        DAP <- DAP + CROP_AGE
        growstages[[as.character(ISTAGE)]]$AGE <- CROP_AGE
        growstages[[as.character(ISTAGE)]]$date <- as.character(END_OF_LEAF_GROWTH_DATE)
        growstages[[as.character(ISTAGE)]]$DOY <- as.integer(format(END_OF_LEAF_GROWTH_DATE, "%j"))
        growstages[[as.character(ISTAGE)]]$SUMDTT <- round(SUMDTT, 1)
        growstages[[as.character(ISTAGE)]]$DAP <- DAP
        break
      }
    }
  } else {
    cat("Error reading weather data for end of pre-anthesis ear growth phase\n")
  }
  # ---------------------------------------------------------------------------------------------- 
  
  # ---------------------------------------------------------------------------------------------- #
  # DETERMINE BEGIN GRAIN FILLING - Grain fill - Start of Grain Filling to Maturity ----
  # ---------------------------------------------------------------------------------------------- #
  ISTAGE <- 5
  # When Silking phase ends and beginning of effective grain
  # filling begins.  Compute grains per plant, ears per pla
  # and barrenness
  
  w <- subset(params$weather, as.Date(date) > END_OF_LEAF_GROWTH_DATE)
  END_GRAIN_FILLING_DATE <- NA
  
  if (nrow(w) > 0) {
    for (i in seq_len(nrow(w))) {
      Tmin <- as.numeric(w$tmin[i])
      Tmax <- as.numeric(w$tmax[i])
      Srad <- as.numeric(w$rad[i])      
      # Duracion del dia
      DOY <- as.integer(format(as.Date(w$date[i]), "%j"))
      TWILEN <- day_length(DOY = DOY, lat = params[['latitude']])
      # Asumimos que la función 'thermal_time_calculation' ya está definida y devuelve un valor numérico
      DTT <- computeDTT(TMIN = Tmin, TMAX = Tmax, TBASE = params[['TT_TBASE']], 
                        TOPT = params[['TT_TEMPERATURE_OPTIMUM']], 
                        ROPT = params[['TT_TEMPERATURE_MAXIMUM']],
                        ISTAGE = ISTAGE, LEAFNO = LEAFNO, XS = 0, 
                        SRAD = Srad, DAYL = TWILEN)
      SUMDTT <- SUMDTT + DTT
      
      if (SUMDTT >= params$P5 * 0.95) {
        END_GRAIN_FILLING_DATE <- as.Date(w$date[i])
        CROP_AGE <- as.integer(END_GRAIN_FILLING_DATE - END_OF_LEAF_GROWTH_DATE)
        DAP <- DAP + CROP_AGE
        growstages[[as.character(ISTAGE)]]$AGE <- CROP_AGE
        growstages[[as.character(ISTAGE)]]$date <- as.character(END_GRAIN_FILLING_DATE)
        growstages[[as.character(ISTAGE)]]$DOY <- as.integer(format(END_GRAIN_FILLING_DATE, "%j"))
        growstages[[as.character(ISTAGE)]]$SUMDTT <- round(SUMDTT, 1)
        growstages[[as.character(ISTAGE)]]$DAP <- DAP
        break
      }
    }
  } else {
    cat("Error reading weather data for beginning of grain fill phase\n")
  }
  # ---------------------------------------------------------------------------------------------- 
  
  # ----------------------------------------------------------------------------------------------
  # DETERMINE END GRAIN FILLING - Maturity -----
  # ----------------------------------------------------------------------------------------------
  ISTAGE <- 6
  
  w <- subset(params$weather, as.Date(date) > END_GRAIN_FILLING_DATE)
  MATURITY_DATE <- NA
  
  if (nrow(w) > 0) {
    for (i in seq_len(nrow(w))) {
      Tmin <- as.numeric(w$tmin[i])
      Tmax <- as.numeric(w$tmax[i])
      Srad <- as.numeric(w$rad[i])      
      # Duracion del dia
      DOY <- as.integer(format(as.Date(w$date[i]), "%j"))
      TWILEN <- day_length(DOY = DOY, lat = params[['latitude']])
      # Asumimos que la función 'thermal_time_calculation' ya está definida y devuelve un valor numérico
      DTT <- computeDTT(TMIN = Tmin, TMAX = Tmax, TBASE = params[['TT_TBASE']], 
                        TOPT = params[['TT_TEMPERATURE_OPTIMUM']], 
                        ROPT = params[['TT_TEMPERATURE_MAXIMUM']],
                        ISTAGE = ISTAGE, LEAFNO = LEAFNO, XS = 0, 
                        SRAD = Srad, DAYL = TWILEN)
      SUMDTT <- SUMDTT + DTT
      
      if (DTT < 2) {
        SUMDTT = params[['P5']]
      }
      
      if (SUMDTT >= params[['P5']]) {
        MATURITY_DATE <- as.Date(w$date[i])
        CROP_AGE <- as.integer(MATURITY_DATE - END_GRAIN_FILLING_DATE)
        DAP <- DAP + CROP_AGE
        growstages[[as.character(ISTAGE)]]$AGE <- CROP_AGE
        growstages[[as.character(ISTAGE)]]$date <- as.character(MATURITY_DATE)
        growstages[[as.character(ISTAGE)]]$DOY <- as.integer(format(MATURITY_DATE, "%j"))
        growstages[[as.character(ISTAGE)]]$SUMDTT <- round(SUMDTT, 1)
        growstages[[as.character(ISTAGE)]]$DAP <- DAP
        break
      }
    }
  } else {
    cat("Error reading weather data for end of grain fill phase\n")
  }
  
  # Transformar la lista en un data frame
  growstages <- do.call(rbind, lapply(growstages, function(x) data.frame(matrix(unlist(x), ncol = length(x), byrow = T), stringsAsFactors = FALSE)))
  
  # Dar nombres a las columnas del data frame
  names(growstages) <- c("istage_old", "istage", "desc", "date", "DOY", "AGE", "DAP", "SUMDTT", "status")
  
  # Convertir columnas a tipos de datos adecuados
  growstages$date <- as.Date(growstages$date)
  growstages$DOY <- as.numeric(growstages$DOY)
  growstages$AGE <- as.numeric(growstages$AGE)
  growstages$DAP <- as.numeric(growstages$DAP)
  growstages$SUMDTT <- as.numeric(growstages$SUMDTT)
  growstages$status <- as.factor(growstages$status)
  
  # Devolver resultados
  return(growstages)
}

# -----------------------------------------------------------------------------




