#rm(list = objects())
# 
# require(dplyr)
# require(purrr)
# require(stringr)
# require(utils)

ReadOverview <- function(archivo.overview, num.stages) {
  contenido.overview <- base::readLines(archivo.overview)
  
  # Runs
  runs <- stringr::str_match(string = contenido.overview, pattern = "\\*RUN\\s+(\\d+)") %>%
    as.data.frame() %>%
    dplyr::filter(! is.na(V1)) %>%
    dplyr::mutate(RUN = as.numeric(V2)) %>%
    dplyr::select(RUN)
    
  # Starting dates
  starting.dates <- stringr::str_match(string = contenido.overview, pattern = "STARTING\\ DATE\\s+:\\s+(\\w{3})\\s+(\\d{1,2})\\s+(\\d{4})") %>%
    as.data.frame() %>%
    dplyr::filter(! is.na(V1)) %>%
    dplyr::mutate(STARTING_DATE = as.Date(sprintf("%s %s %s", V2, V3, V4), format = "%B %d %Y")) %>%
    dplyr::select(STARTING_DATE)
  
  # Planting dates
  #planting.dates <- stringr::str_match(string = contenido.overview, pattern = "PLANTING\\ DATE\\s+:\\s+(\\w{3})\\s+(\\d{1,2})\\s+(\\d{4})") %>%
  #  as.data.frame() %>%
  #  dplyr::filter(! is.na(V1)) %>%
  #  dplyr::mutate(PLANTING_DATE = as.Date(sprintf("%s %s %s", V2, V3, V4), format = "%B %d %Y")) %>%
  #  dplyr::select(PLANTING_DATE)
  
  # treatment \\Treatment-n\\d+\\_H2O-\\d+\\.*\\d*
  treatments <- stringr::str_match(string = contenido.overview, pattern = "\\TREATMENT\\s+\\d+\\s+:\\s+(\\Treatment-n\\d+\\_H2O-\\d+\\.*\\d*)") %>%
    as.data.frame() %>%
    dplyr::filter(! is.na(V1)) %>%
    dplyr::mutate(TREATMENT = V2) %>%
    dplyr::select(TREATMENT)
  
  # Metadatos
  metadata <- dplyr::bind_cols(runs, starting.dates, treatments)
  rm(starting.dates)
  
  # Overview
  lineas.overview <- stringr::str_match(string = contenido.overview,
                                        pattern = "\\*SIMULATED\\ CROP\\ AND\\ SOIL\\ STATUS\\ AT\\ MAIN\\ DEVELOPMENT\\ STAGES") %>%
    as.data.frame() %>%
    dplyr::mutate(es_inicio = ! is.na(V1)) %>%
    dplyr::pull(es_inicio) %>%
    which() + 6
  
  # Leer stages
  stages     <- purrr::map_dfr(
    .x = dplyr::pull(runs, RUN),
    .f = function(run) {
      linea         <- lineas.overview[run]
      rle.linea     <- rle(strsplit(stringr::str_trim(contenido.overview[linea], "right"), "")[[1]])
      offset        <- rle.linea$lengths[1]
      col.widths    <- rle.linea$lengths[which(rle.linea$values == "-")] + rle.linea$lengths[which(rle.linea$values != "-")]
      col.widths[1] <- col.widths[1] + offset
      col.names     <- strsplit(trimws(contenido.overview[linea-1]), split = "\\s+")[[1]]
      run.stages    <- utils::read.fwf(file = textConnection(contenido.overview), 
                                       widths = col.widths, col.names = col.names,
                                       skip = linea, n = num.stages)
      return (dplyr::bind_cols(dplyr::filter(metadata, RUN == run), run.stages) %>%
                dplyr::select(RUN, STARTING_DATE, TREATMENT, DATE, STAGE) %>%
                dplyr::mutate(RUN = as.numeric(RUN)))
    }
  )
  
  stages <- stages %>%
    dplyr::mutate(YDAY = lubridate::yday(as.Date(stages$DATE, format = '%d %B'))) %>%
    dplyr::mutate(YEAR = if_else(YDAY > lubridate::yday(STARTING_DATE) & YDAY < 366,
                                 lubridate::year(STARTING_DATE), lubridate::year(STARTING_DATE) + 1)) %>%
    dplyr::mutate(DATE = as.Date(paste(DATE, YEAR), format = '%d %B %Y')) %>%
    dplyr::select(-YDAY, -YEAR)

  return (stages)
}

#stages <- ReadOverview(archivo.overview = "/home/santiago/OVERVIEW.OUT", num.stages = 15)