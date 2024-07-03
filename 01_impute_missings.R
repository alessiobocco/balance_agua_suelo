
# -----------------------------------------------------------------------------#
# --- PASO 1. Cargar paquetes necesarios                                    ----
# -----------------------------------------------------------------------------#

#rm(list = ls()); gc()
Sys.setenv(TZ = "UTC")
list.of.packages <- c("dplyr", "doMC", "foreach", "iterators", 
                      "parallel", "doParallel", "randomForest", "xts", "zoo",
                      "lazyeval", "sirad", "missForest", "rgeos",
                      "gstat", "geosphere", "rgdal", "optparse")
for (pack in list.of.packages) {
  if (!require(pack, character.only = TRUE)) {
    stop(paste0("Paquete no encontrado: ", pack))
  }
}
rm(pack); gc()

# ------------------------------------------------------------------------------

# -----------------------------------------------------------------------------#
# --- PASO 2. Leer archivo de configuracion                                 ----
# -----------------------------------------------------------------------------#

normalize_dirnames <- function(dirnames) {
  if (is.atomic(dirnames)) 
    dirnames <- base::sub('/$', '', dirnames)
  if (!is.atomic(dirnames))
    for (nm in names(dirnames)) 
      dirnames[[nm]] <- normalize_dirnames(dirnames[[nm]])
  return (dirnames)
}

# a) YAML de configuracion
args <- base::commandArgs(trailingOnly = TRUE)
if (length(args) > 0) {
  archivo.config <- args[1]
} else {
  # No vino el archivo de configuracion por linea de comandos. Utilizo un archivo default
  archivo.config <- paste0(getwd(), "/configuracion.yml")
}
if (! file.exists(archivo.config)) {
  stop(paste0("El archivo de configuración ", archivo.config, " no existe\n"))
} else {
  cat(paste0("Leyendo archivo de configuración ", archivo.config, "...\n"))
  config <- yaml::yaml.load_file(archivo.config)
  config$dir <- normalize_dirnames(config$dir)
}

# b) YAML de parametros
if (length(args) > 1) {
  archivo.params <- args[2]
} else {
  # No vino el archivo de configuracion por linea de comandos. Utilizo un archivo default
  archivo.params <- paste0(getwd(), "/parametros.yml")
}
if (! file.exists(archivo.params)) {
  stop(paste0("El archivo de parámetros ", archivo.params, " no existe\n"))
} else {
  cat(paste0("Leyendo archivo de parámetros ", archivo.params, "...\n"))
  config$params <- yaml::yaml.load_file(archivo.params)
}

# c) YAML de configuración de la API
if (length(args) > 1) {
  archivo.config.api <- args[3]
} else {
  # No vino el archivo de configuracion por linea de comandos. Utilizo un archivo default
  archivo.config.api <- paste0(getwd(), "/configuracion_api.yml")
}
if (! file.exists(archivo.config.api)) {
  stop(paste0("El archivo de configuración ", archivo.config.api, " no existe\n"))
} else {
  cat(paste0("Leyendo archivo de configuración ", archivo.config.api, "...\n"))
  config$api <- yaml::yaml.load_file(archivo.config.api)$api
}

# d) Preparar carpetas de salida, crearlas si no existen
config$dir$outputs <- glue::glue('{config$dir$outputs}/{config$params$run_id}')
if (!fs::dir_exists(config$dir$outputs))
  fs::dir_create(config$dir$outputs)
config$dir$weather <- glue::glue('{config$dir$weather}/{config$params$run_id}')
if (!fs::dir_exists(config$dir$weather))
  fs::dir_create(config$dir$weather)

# ------------------------------------------------------------------------------

# -----------------------------------------------------------------------------#
# --- PASO 3. Cargar librerias                                              ----
# -----------------------------------------------------------------------------#

# a) Cargar funciones necesarias
source(glue::glue("./lib/helpers.R"), echo = FALSE)
source(glue::glue("./lib/crc-api.R"), echo = FALSE)
source(glue::glue("./lib/impute/Estimar.R"), echo = FALSE)
source(glue::glue("./lib/impute/Impute.R"), echo = FALSE)
source(glue::glue("./lib/impute/Params.R"), echo = FALSE)
source(glue::glue("./lib/impute/Utils.R"), echo = FALSE)

# b) Variable para almacenar los posibles errores
errors <- c()

# ------------------------------------------------------------------------------

# -----------------------------------------------------------------------------#
# --- PASO 4. Leer/procesar variables de entrada                            ----
# -----------------------------------------------------------------------------#

# Datos de entrada
años <- base::seq(from = config$params$weather$start_year, to = config$params$weather$end_year)
suelos <- purrr::pmap_dfr(
  .l = utils::tail(config$params$soils, -1) %>% purrr::transpose(), 
  .f = function(..., nombres_columnas) { 
    tibble::tibble(..., .name_repair = "minimal") %>% setNames(nombres_columnas) },
  nombres_columnas = config$params$soils[[1]]
) 

# Creamos un array con los OMM ID's de las estaciones sobre las que queremos estimar.
estacionesID <- suelos %>% dplyr::distinct(weather_id) %>% dplyr::pull()

# Control número de estaciones
stopifnot(length(estacionesID) > 0)

# ------------------------------------------------------------------------------

# -----------------------------------------------------------------------------#
# --- PASO 5. Iniciar proceso de imputación                                 ----
# -----------------------------------------------------------------------------#

httr::set_config( httr::config(ssl_verifypeer = FALSE) )

# Obtenemos los datos de las estaciones.
estaciones <- purrr::map_dfr(
  .x = config$params$weather$countries,
  .f = function(pais) {
    est <- ConsumirServicioJSON(url = paste0(config$api$url, glue::glue("/estaciones/{pais}")),
                                usuario = config$api$user, clave = config$api$pass) 
    est <- est %>%
      dplyr::mutate(pais_id = pais) 
    return (est)
  }
)


# Renombrar columnas y filtrar estaciones
estaciones <- estaciones %>%
  dplyr::rename(lon_dec = longitud, lat_dec = latitud, elev = elevacion) %>%
  dplyr::filter(omm_id %in% c(87244, 87328, 87344, 87345, 87349, 87466, 87453, 87467, 87534, 9987009)) 

# Calculamos las coordenadas en Gauss Kruger de cada estación.
GK.coords <- Gauss.Kruger.coordinates(estaciones)
# Agregamos las coordenadas x e y (GK) de cada estación.
estaciones <- data.frame(sp::coordinates(GK.coords), estaciones) %>%
  dplyr::select(x = 1, y = 2)

# Obtenemos los registros diarios de las estaciones con las que vamos a trabajar.
registrosDiarios <- purrr::map_dfr(
  .x = estaciones$omm_id,
  .f = function(omm_id) {
    fecha.desde           <- ConvertirFechaISO8601(as.Date(glue::glue("{config$params$weather$start_year}-01-01"), tz = UTC))
    fecha.hasta           <- ConvertirFechaISO8601(as.Date(glue::glue("{config$params$weather$end_year}-12-31"), tz = UTC))
    url.registros.diarios <- glue::glue("{config$api$url}/registros_diarios/{omm_id}/{fecha.desde}/{fecha.hasta}")
    registros.largo       <- ConsumirServicioJSON(url = url.registros.diarios,
                                                  usuario = config$api$user, clave = config$api$pass) %>%
      {
        if("estado" %in% names(.)) dplyr::filter(., estado == 'A') else .
        if("estado" %in% names(.)) dplyr::select(., -estado) else .
      }
    registros.ancho       <- tidyr::spread(registros.largo, key = variable_id, value = valor)
    return (registros.ancho)
  }
)

# Convertimos las fechas de la query en variables Date.
registrosDiarios$fecha <- as.Date(registrosDiarios$fecha)

# Check invalid temperatures to avoid DSSAT errors.
wrong_temps <- which(registrosDiarios$tmax <= registrosDiarios$tmin)
# Make invalid temperatures NA so the impute methods recalculate them.
registrosDiarios[wrong_temps, c('tmax', 'tmin')] <- NA

registrosImputados <- purrr::map_dfr(
  .x = estaciones$omm_id,
  .f = function(omm_id) {
    
    datosEstacion <- registrosDiarios[registrosDiarios$omm_id == omm_id,]
    indexesToWrite <- c()
    
    # Obtenemos los ID's y las coordenadas de cada estación vecina.max.distancia <- 150
    max.diff.elev <- 300
    max.vecinas   <- 15
    url.vecinas   <- glue::glue("{config$api$url}/estaciones_vecinas/{omm_id}")
    query.vecinas <- glue::glue("max_diferencia_elevacion={max.diff.elev}&max_vecinas={max.vecinas}")
    vecinos.data  <- ConsumirServicioJSON(url = glue::glue("{url.vecinas}?{query.vecinas}"),
                                          usuario = config$api$user, clave = config$api$pass) %>%
      {
        if("estado" %in% names(.)) dplyr::filter(., estado == 'A') else .
        if("estado" %in% names(.)) dplyr::select(., -estado) else .
      }
    
    # renombrar columnas
    vecinos.data <- vecinos.data %>%
      dplyr::filter(tipo == 'C', !omm_id %in% c(87146, 9987004, 9987005)) %>%
      dplyr::rename(lon_dec = longitud, lat_dec = latitud, elev = elevacion) %>%
      dplyr::select(omm_id, lat_dec, lon_dec)
    
    # Calculamos las coordenadas en GK.
    GK.coords <- Gauss.Kruger.coordinates(vecinos.data)
    
    # Agregamos las coordenadas x e y al data frame de vecinos.
    vecinos.data <- data.frame(sp::coordinates(GK.coords), omm_id=GK.coords@data$omm_id) %>%
      dplyr::select(x = 1, y = 2)
    
    # Traer registros de vecinos (solo cuando sea necesario).
    registrosVecinos <- NULL
    
    #
    # Imputar tmax, tmin y prcp
    #
    
    for (variable in c('tmax', 'tmin', 'prcp')) {
      missingIndexes <- which(is.na(datosEstacion[, variable]))
      
      writeLines(paste0("> Station: ", omm_id, ". Variable: ", variable, ". Missing: ", length(missingIndexes)))
      
      # Check if there are missing values for this station and variable, otherwise, skip it.
      if(length(missingIndexes) == 0) next;
      
      indexesToWrite <- c(indexesToWrite, missingIndexes)
      
      # Query the for neighbor's data. Done only if there are missing values.
      if(is.null(registrosVecinos)) {
        registrosVecinos <- purrr::map_dfr(
          .x = vecinos.data$omm_id,
          .f = function(omm_id) {
            fecha.desde           <- ConvertirFechaISO8601(as.Date(glue::glue("{config$params$weather$start_year}-01-01"), tz = UTC))
            fecha.hasta           <- ConvertirFechaISO8601(as.Date(glue::glue("{config$params$weather$end_year}-12-31"), tz = UTC))
            url.registros.diarios <- glue::glue("{config$api$url}/registros_diarios/{omm_id}/{fecha.desde}/{fecha.hasta}")
            registros.largo       <- ConsumirServicioJSON(url = url.registros.diarios,
                                                          usuario = config$api$user, clave = config$api$pass) %>% 
              {
                if("estado" %in% names(.)) dplyr::filter(., estado == 'A') else .
                if("estado" %in% names(.)) dplyr::select(., -estado) else .
              }

            registros.ancho       <- tidyr::spread(registros.largo, key = variable_id, value = valor) 
            return (registros.ancho)
          }
        )
        registrosVecinos$fecha <- as.Date(registrosVecinos$fecha)
        registrosVecinos <- tibble::as_tibble(registrosVecinos) %>% dplyr::left_join(vecinos.data, by='omm_id') 
      }
      
      if (variable == 'prcp') {
        datosEstacion <- impute_mf(omm_id, datosEstacion, variable, estaciones, missingIndexes, registrosVecinos, vecinos.data, config$max.procesos)
      } else {
        datosEstacion <- impute_idw(omm_id, datosEstacion, variable, estaciones, missingIndexes, registrosVecinos, vecinos.data, config$max.procesos)
      }
    }
    
    #
    # Imputar radiación
    #
    
    # Update radiation only.
    missing_dates <- as.Date(datosEstacion$fecha)
    
    writeLines(paste0("> Station: ", omm_id, ". Variable: rad. Missing: ", length(missing_dates)))
    
    codigo_pais <- toupper(estaciones[estaciones$omm_id == omm_id, 'pais_id'])
    
    if(!codigo_pais %in% names(srad_parameters)) codigo_pais <- 'AR'
    
    records <- datosEstacion %>% filter(fecha %in% missing_dates)
    
    estimado <- estimarRadiacion(estaciones=estaciones[estaciones$omm_id == omm_id, ],
                                 registrosDiarios=records,
                                 ap.cal = srad_parameters[[codigo_pais]]$ap,
                                 bc.cal = srad_parameters[[codigo_pais]]$bc,
                                 svk.cal = srad_parameters[[codigo_pais]]$svk)
    
    if(estimado$not_estimated > 0) {
      error_details <- paste0("Failed to estimate ", estimado$not_estimated, " radiation values for station ", omm_id)
      errors <<- c(errors, error_details)
      writeLines(error_details)
      return (NULL)
    } else {
      return (estimado$results %>% dplyr::select(omm_id, fecha, tmax, tmin, prcp, helio, nub, rad, metodo))
    }
  }
) %>%
  dplyr::mutate(rad = dplyr::if_else(rad < 2, 2, rad))

# ------------------------------------------------------------------------------

# -----------------------------------------------------------------------------#
# --- PASO 6. Se guardan los resultados obtenidos                           ----
# -----------------------------------------------------------------------------#

# Se guarda el resultado en un archivo .rds
filename <- stringi::stri_sub(config$params$weather$file, 1, -5)
saveRDS(registrosImputados, glue::glue("{config$dir$weather}/{filename}.rds"))

# Modificar formato para poder correr simulaciones
registrosImputadosOUT <- registrosImputados %>%
  dplyr::mutate(realization = 1,
                rad = dplyr::if_else(rad < 2, 2, rad)) %>%
  dplyr::select(station_id = omm_id, date = fecha, realization,
                tmax, tmin, prcp, rad) 


# Guardar como csv, si no existe un archivo con el mismo nombre
filename <- glue::glue("{config$dir$weather}/{config$params$weather$file}")
if (!fs::file_exists(filename)) {
  write.csv(registrosImputadosOUT, filename, row.names = F)
} else {
  warning(glue::glue("No se pudo guardar el archivo {filename}, ya existe ",
                     "un archivo con ese nombre en esa ubicación!"))
}

# ------------------------------------------------------------------------------


