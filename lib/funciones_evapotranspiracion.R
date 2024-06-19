# ---------------------------------------------------------------------------- #
# Funciones auxiliares
# ---------------------------------------------------------------------------- #
# tmax = evapotranspiracion.estacion$tmax[277]
# tmin = evapotranspiracion.estacion$tmin[277]
# tmed = evapotranspiracion.estacion$tmed[277]
# td = NA
# HRmed = evapotranspiracion.estacion$hr[277]
# HRmax = NA
# HRmin = NA
# Rs = NA
# presion = evapotranspiracion.estacion$pres_est[277]/10
# helio = evapotranspiracion.estacion$helio[277]
# u2 = evapotranspiracion.estacion$vmed[277]
# date = evapotranspiracion.estacion$fecha[277]

# Función de Pettit
pettitt.test <- function(x){    
  na.fail(x)
  n <- length(x)
  DNAME <- deparse(substitute(x))
  k <- 1:n
  r <- rank(x)
  Uk <- sapply(k, function(x) 2 * sum(r[1:x]) - x*(n+1))
  Uka <- abs(Uk)
  U <- max(Uka)
  K <- min(k[Uka == U])
  pval <- 2.0 * exp(( -6.0 * U^2) / (n^3 + n^2))
  
  if (is.ts(x)){
    fr <- frequency(x)
    st <- start(x)
    ed <- end(x)
    Uk <- ts(Uk, start=st, end = ed, frequency= fr)
  }
  
  attr(Uk, 'nm') <- "Uk"
  
  names(K) <- "probable change point at time K"
  retval <- list(nobs = n, 
    statistic = c("U*" = U),
    estimate= K,
    p.value =  pval,
    data.name= DNAME,
    alternative="two.sided",
    data = Uk,
    method = "Pettitt's test for single change-point detection")
  class(retval) <- c("htest", "cptest")
  return(retval)
}

# Calculo del anio agricola o hidrologico
fecha_anio_agricola <- function(dates, start_month=5) {
  
  # Convert dates into POSIXlt
  dates.posix = as.POSIXlt(dates)
  # Year offset
  offset = ifelse(dates.posix$mon >= start_month - 1, 1, 0)
  # Agri year
  adj.year = dates.posix$year + 1900 + offset - 1
  # Return the agri year
  adj.year
}

anio_agricola_fecha <- function(hydro_year, 
                                month,
                                start_month=5) {
  
  # Year offset
  offset = ifelse(month >= start_month, 0, 1)
  # Agri year
  adj.year = hydro_year + offset 
  # Return the agri year
  adj.year
}

# Evapotranspiracion
evapotranspiracion.penman.monteith.fao <- function(tmax = NA, tmin = NA, tmed = NA,
  td = NA, HRmed = NA, HRmax = NA, HRmin = NA, u2 = NA, helio = NA, Rs = NA,
  presion = NA, albedo = 0.23, a = NA, b = NA, z = NA, lat = NA,
  k = 0.16, date = NA, terms = FALSE) {
  
  # tmax: temperatura máxima diaria [°C]
  # tmin: temperatura mínima diaria [°C]
  # tmed: temperatura media diaria [°C]
  # td: temperatura del punto de rocío [°C]
  # HRmed: humedad relativa media diaria [%]
  # HRmax: humedad relativa diaria máxima [%]
  # HRmin: humedad relativa mínima [%]
  # u2: velocidad media del viento a 2 m de altura [m/s]
  # helio: heliofania, duración del brillo solar [horas]
  # Rs: radiación solar de onda corda [MJ.m^(-2).day^(-1)]
  # presion: presion atmosférica de la estacion [kPa]
  # albedo: albedo de ls superficie. 0.23 por defecto.
  # a y b: coeficientes de la ecuación de Armstrong para la estimación de radiación global
  # z: altura de la estacion [m]
  # lat: latitud de la estacion en grados
  # k: coeficiente de ajuste modelo de radiación de Hargreaves. 0.16 tierra firme, 0.19 zonas costeras [°C^-0.5]
  # date: fecha del día en formato de fecha
  # terms: devuelve los terminos radiativo y aerodinamico por seprado en mm
  
  ##################
  require(dplyr)
  require(sirad)
  require(lubridate)
  require(gtools)
  
  ##################
  
  # Si la temperatura media no está disponible, calcular
  if (is.na(tmed)) {
    
    tmed <- (tmin+tmax)/2  ## Eq. 9 Allen & al. 1998.
    
  }
  
  # Si la humedad relativa media no está disponible, calcular
  if (is.na(HRmed) && !is.na(HRmax) || !is.na(HRmin)) {
    
    HRmed = (HRmin+HRmax)/2
  }
  
  # ---- Término aerodinámico ----
  
  # Pendiente de la curva de presión de vapor
  ## Eq. 13 Allen & al. 1998 [kPa/°C]
  pendiente.curva.presion.vapor <- 4098*0.6108*exp(17.27*tmed/(tmed+237.3))/((tmed+237.3)^2)
  
  # Pendiente psicrometrica
  if (is.na(presion)) {
    ## Eq. 7 Allen & al. 1998 [kPa]
    presion <- 101.3*(((293-0.0065*z)/293)^5.26) 
  }
  
  # Calor específico a presión constante
  Cp <- 0.001013 ## Allen & al. 1998 p. 32
  # peso molecular del agua sobre el del aire seco
  epsilon <- 0.622 ## Allen & al. 1998 p. 32
  # Calor latente de vaporización
  lambda <- 2.45 ## Allen & al. 1998 p. 31
  # Constante psicrometrica
  constante.psicrometrica <- Cp*presion/(epsilon*lambda) ## Eq. 8 Allen & al. 1998.
  
  # Presión de vapor de saturacion
  eTmin <- 0.6108*exp(17.27*tmin/(tmin+237.3)) ## min saturation vapour pressure [kPa] Eq. 11 Allen & al. 1998.
  eTmax <- 0.6108*exp(17.27*tmax/(tmax+237.3)) ## max saturation vapour pressure [kPa] Eq. 11 Allen & al. 1998.
  es <- (eTmax+eTmin)/2 ## Eq. 12 Allen & al. 1998.
  
  # Presion de vapor actual
  if(!is.na(HRmax) && !is.na(HRmin)) {
    ea <- (eTmin*RHmax+eTmax*RHmin)/200 ## Eq. 17 Allen & al. 1998.
  } else if (!is.na(HRmax)) {
    ea <- eTmin*HRmax/100 ## Eq. 18 Allen & al. 1998
  } else if (!is.na(HRmed)) {
    ea <- HRmed*(eTmax+eTmin)/200 ## Eq. 19 Allen & al. 1998
  } else if (!is.na(td)) {
    ea <- 0.6108*exp(17.27*td/(td+237.3)) ## Eq. 48 Allen & al. 1998, but see limitations p 58. Should be checked against measured
  } else {
    ea <- NA
  }
  
  # Défict de presion de vapor
  deficit.vapor <- es - ea
  
  # Suma de términos
  termino.aerodinamico <- (constante.psicrometrica*900*u2*(deficit.vapor)/(tmed+273))/
    (pendiente.curva.presion.vapor+constante.psicrometrica*(1+0.34*u2))  ## Eq. 6 Allen & al. 1998.
  
  # ---- Término radiativo ----
  
  # Balance de radiación
  # Radiación estraterrestre
  # constante solar = 0.0820 MJ.m^(-2).min^(-1)
  Gsc <- 0.0820
  # convertir latitudes en grados a radianes
  phi <- sirad::radians(lat)           # Eq. 23 Allen & al. 1998.
  # dia juliano correspondiente al dia i
  dia.juliano <- lubridate::yday(date)
  # Inversa de la distancia entre la Tierra y el sol
  Dr <- 1+0.033*cos(2*pi*dia.juliano/365)      # Eq. 23 Allen & al. 1998.
  # inclinacion solar [rad]
  delta <- 0.409*sin((2*pi*dia.juliano/365)-1.39) # Eq. 24 Allen & al. 1998.
  # angulo solar a poniente [rad]
  Ws <- acos(-tan(phi)*tan(delta))        # Eq. 25 Allen & al. 1998.
  
  # Extraterrestrial radiation for daily periods [MJ.m^(-2).day^(-1)]
  eRad <- (24*60/pi)*Gsc*Dr*(Ws*sin(phi)*sin(delta)+cos(phi)*cos(delta)*sin(Ws)) # Eq. 21 Allen & al. 1998.
  
  # Radiación solar o global
  if (!is.na(Rs)) {
    Rs <- Rs
  } else if (!is.na(helio)) {
    helio.max <- 24/pi*Ws # Eq. 34 Allen & al. 1998.
    Rs <- (a+b*(helio/helio.max))*eRad # Eq. 35 Allen & al. 1998
  } else if (!is.na(tmax) && !is.na(tmin)) {
    Rs <- sqrt(tmax-tmin)*eRad*k # Eq. 50 Allen & al. 1998 
  } else {
    Rs <- NA
  }
  
  # Radiacion solar neta
  Rns <- (1-albedo)*Rs ## Eq. 38 Allen & al. 1998. albedo = 0.23 for the FAO hypothetical crop reference
  # Radiación solar en un día despejado
  Rso <- (0.75+2*10^(-5)*z)*eRad ## Eq. 37 Allen & al. 1998
  # Radiacion neta de onda larga
  SteBolCon <- 4.903*10^(-9) # Constante de Stefan-Boltzmann [4.903.10^(-9) MJ.K^(-4).m^(-2).day^(-1)]
  Rnl <- SteBolCon*(((tmax+273.16)^4+(tmin+273.16)^4)/2)*(0.34 - 0.14*sqrt(ea))*((1.35*ifelse((Rs/Rso)>1,1,Rs/Rso))-0.35) ## Eq. 39 Allen & al. 1998.
  # Radiacion neta sobre la superficie del cultivo
  Rn <- Rns-Rnl ## Eq. 40 Allen & al. 1998.
  # Flujo de calor del suelo
  G <-0 ## Eq. 42 Allen & al. 1998. G puede ser ignorado para calculo diarios
  
  # ---- Suma de terminos 
  termino.radiativo<- (0.408*pendiente.curva.presion.vapor*(Rn-G))/
    (pendiente.curva.presion.vapor+constante.psicrometrica*(1+0.34*u2))
  
  # ---- Evapotranspiración ----
  eto <- termino.radiativo + termino.aerodinamico
  
  if (terms) {
    
    terminos <- c(termino.aerodinamico, termino.radiativo)
    names(terminos) <- c('Eaero', 'Erad')
    
    return(terminos)
    
  } else {
    
    return(eto)
    
  }
}

# Evapotranspiracion
eto.penman.monteith.fao <- function(X, constantes) {
  
  # tmax: temperatura máxima diaria [°C]
  # tmin: temperatura mínima diaria [°C]
  # tmed: temperatura media diaria [°C]
  # td: temperatura del punto de rocío [°C]
  # HRmed: humedad relativa media diaria [%]
  # HRmax: humedad relativa diaria máxima [%]
  # HRmin: humedad relativa mínima [%]
  # u2: velocidad media del viento a 2 m de altura [m/s]
  # helio: heliofania, duración del brillo solar [horas]
  # Rs: radiación solar de onda corda [MJ.m^(-2).day^(-1)]
  # presion: presion atmosférica de la estacion [kPa]
  # albedo: albedo de ls superficie. 0.23 por defecto.
  # a y b: coeficientes de la ecuación de Armstrong para la estimación de radiación global
  # z: altura de la estacion [m]
  # lat: latitud de la estacion en grados
  # k: coeficiente de ajuste modelo de radiación de Hargreaves. 0.16 tierra firme, 0.19 zonas costeras [°C^-0.5]
  # date: fecha del día en formato de fecha
  # terms: devuelve los terminos radiativo y aerodinamico por seprado en mm
  tmax <- X$tmax
  tmin <- X$tmin
  u2 <- X$vmed
  HRmed <- X$hr
  Rs = X$rad
  date <- constantes$date
  lat <- constantes$lat
  albedo = 0.23
  k = 0.16
  z <- constantes$z
  # Calculo de temperatura media
  tmed <- (tmin+tmax)/2  ## Eq. 9 Allen & al. 1998.
  
  # ---- Término aerodinámico ----
  
  # Pendiente de la curva de presión de vapor
  ## Eq. 13 Allen & al. 1998 [kPa/°C]
  pendiente.curva.presion.vapor <- 4098*0.6108*exp(17.27*tmed/(tmed+237.3))/((tmed+237.3)^2)
  
  # Pendiente psicrometrica
  ## Eq. 7 Allen & al. 1998 [kPa]
  presion <- 101.3*(((293-0.0065*z)/293)^5.26) 
  
  
  # Calor específico a presión constante
  Cp <- 0.001013 ## Allen & al. 1998 p. 32
  # peso molecular del agua sobre el del aire seco
  epsilon <- 0.622 ## Allen & al. 1998 p. 32
  # Calor latente de vaporización
  lambda <- 2.45 ## Allen & al. 1998 p. 31
  # Constante psicrometrica
  constante.psicrometrica <- Cp*presion/(epsilon*lambda) ## Eq. 8 Allen & al. 1998.
  
  # Presión de vapor de saturacion
  eTmin <- 0.6108*exp(17.27*tmin/(tmin+237.3)) ## min saturation vapour pressure [kPa] Eq. 11 Allen & al. 1998.
  eTmax <- 0.6108*exp(17.27*tmax/(tmax+237.3)) ## max saturation vapour pressure [kPa] Eq. 11 Allen & al. 1998.
  es <- (eTmax+eTmin)/2 ## Eq. 12 Allen & al. 1998.
  
  # Presion de vapor actual
  ea <- HRmed*(eTmax+eTmin)/200 ## Eq. 19 Allen & al. 1998
  
  # Défict de presion de vapor
  deficit.vapor <- es - ea
  
  # Suma de términos
  termino.aerodinamico <- (constante.psicrometrica*900*u2*(deficit.vapor)/(tmed+273))/
    (pendiente.curva.presion.vapor+constante.psicrometrica*(1+0.34*u2))  ## Eq. 6 Allen & al. 1998.
  
  # ---- Término radiativo ----
  
  # Balance de radiación
  # Radiación estraterrestre
  # constante solar = 0.0820 MJ.m^(-2).min^(-1)
  Gsc <- 0.0820
  # convertir latitudes en grados a radianes
  phi <- sirad::radians(lat)           # Eq. 23 Allen & al. 1998.
  # dia juliano correspondiente al dia i
  #dia.juliano <- lubridate::yday(date)
  dia.juliano <- date
  # Inversa de la distancia entre la Tierra y el sol
  Dr <- 1+0.033*cos(2*pi*dia.juliano/365)      # Eq. 23 Allen & al. 1998.
  # inclinacion solar [rad]
  delta <- 0.409*sin((2*pi*dia.juliano/365)-1.39) # Eq. 24 Allen & al. 1998.
  # angulo solar a poniente [rad]
  Ws <- acos(-tan(phi)*tan(delta))        # Eq. 25 Allen & al. 1998.
  
  # Extraterrestrial radiation for daily periods [MJ.m^(-2).day^(-1)]
  eRad <- (24*60/pi)*Gsc*Dr*(Ws*sin(phi)*sin(delta)+cos(phi)*cos(delta)*sin(Ws)) # Eq. 21 Allen & al. 1998.
  
  # Radiación solar o global
  Rs <- Rs
  
  # Radiacion solar neta
  Rns <- (1-albedo)*Rs ## Eq. 38 Allen & al. 1998. albedo = 0.23 for the FAO hypothetical crop reference
  # Radiación solar en un día despejado
  Rso <- (0.75+2*10^(-5)*z)*eRad ## Eq. 37 Allen & al. 1998
  # Radiacion neta de onda larga
  SteBolCon <- 4.903*10^(-9) # Constante de Stefan-Boltzmann [4.903.10^(-9) MJ.K^(-4).m^(-2).day^(-1)]
  Rnl <- SteBolCon*(((tmax+273.16)^4+(tmin+273.16)^4)/2)*(0.34 - 0.14*sqrt(ea))*((1.35*ifelse((Rs/Rso)>1,1,Rs/Rso))-0.35) ## Eq. 39 Allen & al. 1998.
  # Radiacion neta sobre la superficie del cultivo
  Rn <- Rns-Rnl ## Eq. 40 Allen & al. 1998.
  # Flujo de calor del suelo
  G <- 0 ## Eq. 42 Allen & al. 1998. G puede ser ignorado para calculo diarios
  
  # ---- Suma de terminos 
  termino.radiativo<- (0.408*pendiente.curva.presion.vapor*(Rn-G))/
    (pendiente.curva.presion.vapor+constante.psicrometrica*(1+0.34*u2))
  
  # ---- Evapotranspiración ----
  eto <- termino.radiativo + termino.aerodinamico
  
  return(eto)
  
}


# evapotranspiracion.penman.monteith.fao(tmax = tmax, tmin = tmin, tmed = tmed, 
#  HRmed = HRmed, helio = helio, u2 = u2, presion = presion, date = date, z = 70, lat = -32,
#   a = a, b = b)


evapotranspiracion.hargreaves <- function(tmax = NULL, tmin = NULL, tmed = NULL,
  td = NULL, helio = NULL, Rs = NULL, z = NULL, lat = NULL, date = NULL) {
  
  # tmax: temperatura máxima diaria [°C]
  # tmin: temperatura mínima diaria [°C]
  # tmed: temperatura media diaria [°C]
  # z: altura de la estacion [m]
  # lat: latitud de la estacion en grados
  # date: fecha del día en formato de fecha

  ##################
  require(dplyr)
  require(sirad)
  require(lubridate)
  require(gtools)
  
  ##################
  
  # Si la temperatura media no está disponible, calcular
  if (gtools::invalid(tmed)) {
    
    tmed <- (tmin+tmax)/2  ## Eq. 9 Allen & al. 1998.
    
  }
  
  # Balance de radiación
  # Calor latente de vaporización
  lambda <- 2.45 ## Allen & al. 1998 p. 31
  # Radiación estraterrestre
  # constante solar = 0.0820 MJ.m^(-2).min^(-1)
  Gsc <- 0.0820
  # convertir latitudes en grados a radianes
  phi <- sirad::radians(lat)           # Eq. 23 Allen & al. 1998.
  # dia juliano correspondiente al dia i
  dia.juliano <- lubridate::yday(date)
  # Inversa de la distancia entre la Tierra y el sol
  Dr <- 1+0.033*cos(2*pi*dia.juliano/365)      # Eq. 23 Allen & al. 1998.
  # inclinacion solar [rad]
  delta <- 0.409*sin((2*pi*dia.juliano/365)-1.39) # Eq. 24 Allen & al. 1998.
  # angulo solar a poniente [rad]
  Ws <- acos(-tan(phi)*tan(delta))        # Eq. 25 Allen & al. 1998.
  
  # Extraterrestrial radiation for daily periods [MJ.m^(-2).day^(-1)]
  eRad <- (24*60/pi)*Gsc*Dr*(Ws*sin(phi)*sin(delta)+cos(phi)*cos(delta)*sin(Ws)) # Eq. 21 Allen & al. 1998.
  
  # ---- Evapotranspiración ----
  eto <- 0.0023*(tmed+17.8)*sqrt(tmax-tmin)*eRad/lambda ## Eq. 52 Allen & al. 1998; multiplied by 1/lambda to convert from MJ.m^(-2).day^(-1) to mm.day.day^(-1)
  
  return(eto)
  
}


evapotranspiracion.pristley.taylor <- function(tmax = NULL, tmin = NULL, tmed = NULL,
  td = NULL, HRmed = NULL, HRmax = NULL, HRmin = NULL, helio = NULL, Rs = NULL,
  presion = NULL, albedo = 0.23, a = NULL, b = NULL, z = NULL, lat = NULL,
  k = 0.16, alpha = 1.26, date = NULL) {
  
  # tmax: temperatura máxima diaria [°C]
  # tmin: temperatura mínima diaria [°C]
  # tmed: temperatura media diaria [°C]
  # td: temperatura del punto de rocío [°C]
  # HRmed: humedad relativa media diaria [%]
  # HRmax: humedad relativa diaria máxima [%]
  # HRmin: humedad relativa mínima [%]
  # helio: heliofania, duración del brillo solar [horas]
  # Rs: radiación solar de onda corda [MJ.m^(-2).day^(-1)]
  # presion: presion atmosférica de la estacion [kPa]
  # albedo: albedo de ls superficie. 0.23 por defecto.
  # a y b: coeficientes de la ecuación de Armstrong para la estimación de radiación global
  # z: altura de la estacion [m]
  # lat: latitud de la estacion en grados
  # k: coeficiente de ajuste modelo de radiación de Hargreaves. 0.16 tierra firme, 0.19 zonas costeras [°C^-0.5]
  # alpha: constante
  # date: fecha del día en formato de fecha

  ##################
  require(dplyr)
  require(sirad)
  require(lubridate)
  require(gtools)
  
  ##################
  
  # Si la temperatura media no está disponible, calcular
  if (gtools::invalid(tmed)) {
    tmed <- (tmin+tmax)/2  ## Eq. 9 Allen & al. 1998.
  }
  
  # Si la humedad relativa media no está disponible, calcular
  if (gtools::invalid(HRmed) & !gtools::invalid(HRmax) | !gtools::invalid(HRmin)) {
    HRmed = (HRmin+HRmax)/2
  }
  
  # Pendiente de la curva de presión de vapor
  ## Eq. 13 Allen & al. 1998 [kPa/°C]
  pendiente.curva.presion.vapor <- 4098*0.6108*exp(17.27*tmed/(tmed+237.3))/((tmed+237.3)^2)
  
  # Pendiente psicrometrica
  if (gtools::invalid(presion)) {
    ## Eq. 7 Allen & al. 1998 [kPa]
    presion <- 101.3*(((293-0.0065*z)/293)^5.26) 
  }
  
  # Calor específico a presión constante
  Cp <- 0.001013 ## Allen & al. 1998 p. 32
  # peso molecular del agua sobre el del aire seco
  epsilon <- 0.622 ## Allen & al. 1998 p. 32
  # Calor latente de vaporización
  lambda <- 2.45 ## Allen & al. 1998 p. 31
  # Constante psicrometrica
  constante.psicrometrica <- Cp*presion/(epsilon*lambda) ## Eq. 8 Allen & al. 1998.
  
  # Presión de vapor de saturacion
  eTmin <- 0.6108*exp(17.27*tmin/(tmin+237.3)) ## min saturation vapour pressure [kPa] Eq. 11 Allen & al. 1998.
  eTmax <- 0.6108*exp(17.27*tmax/(tmax+237.3)) ## max saturation vapour pressure [kPa] Eq. 11 Allen & al. 1998.
  es <- (eTmax+eTmin)/2 ## Eq. 12 Allen & al. 1998.
  
  # Presion de vapor actual
  if(!gtools::invalid(HRmax) & !gtools::invalid(HRmin)) {
    ea <- (eTmin*RHmax+eTmax*RHmin)/200 ## Eq. 17 Allen & al. 1998.
  } else if (!gtools::invalid(HRmax) & gtools::invalid(HRmin)) {
    ea <- eTmin*HRmax/100 ## Eq. 18 Allen & al. 1998
  } else if (!gtools::invalid(HRmed)) {
    ea <- HRmed*(eTmax+eTmin)/200 ## Eq. 19 Allen & al. 1998
  } else { (!gtools::invalid(td)) 
    ea <- 0.6108*exp(17.27*tmin/(tmin+237.3)) ## Eq. 48 Allen & al. 1998, but see limitations p 58. Should be checked against measured
  }
  
  # Deífict de presion de vapor
  deficit.vapor <- es - ea
  
  # Balance de radiación
  # Radiación estraterrestre
  # constante solar = 0.0820 MJ.m^(-2).min^(-1)
  Gsc <- 0.0820
  # convertir latitudes en grados a radianes
  phi <- sirad::radians(lat)           # Eq. 23 Allen & al. 1998.
  # dia juliano correspondiente al dia i
  dia.juliano <- lubridate::yday(date)
  # Inversa de la distancia entre la Tierra y el sol
  Dr <- 1+0.033*cos(2*pi*dia.juliano/365)      # Eq. 23 Allen & al. 1998.
  # inclinacion solar [rad]
  delta <- 0.409*sin((2*pi*dia.juliano/365)-1.39) # Eq. 24 Allen & al. 1998.
  # angulo solar a poniente [rad]
  Ws <- acos(-tan(phi)*tan(delta))        # Eq. 25 Allen & al. 1998.
  
  # Extraterrestrial radiation for daily periods [MJ.m^(-2).day^(-1)]
  eRad <- (24*60/pi)*Gsc*Dr*(Ws*sin(phi)*sin(delta)+cos(phi)*cos(delta)*sin(Ws)) # Eq. 21 Allen & al. 1998.
  
  # Radiación solar o global
  if (!gtools::invalid(Rs)) {
    Rs <- Rs
  } else if (!gtools::invalid(helio) & gtools::invalid(Rs)) {
    helio.max <- 24/pi*Ws # Eq. 34 Allen & al. 1998.
    Rs <- (a+b*(helio/helio.max))*eRad # Eq. 35 Allen & al. 1998
  } else { (!gtools::invalid(tmax) & !gtools::invalid(tmin) & gtools::invalid(Rs) | gtools::invalid(helio)) 
    Rs <- sqrt(tmax-tmin)*eRad*k # Eq. 50 Allen & al. 1998 
  }
  
  # Radiacion solar neta
  Rns <- (1-albedo)*Rs ## Eq. 38 Allen & al. 1998. albedo = 0.23 for the FAO hypothetical crop reference
  # Radiación solar en un día despejado
  Rso <- (0.75+2*10^(-5)*z)*eRad ## Eq. 37 Allen & al. 1998
  # Radiacion neta de onda larga
  SteBolCon <- 4.903*10^(-9) # Constante de Stefan-Boltzmann [4.903.10^(-9) MJ.K^(-4).m^(-2).day^(-1)]
  Rnl <- SteBolCon*(((tmax+273.16)^4+(tmin+273.16)^4)/2)*(0.34 - 0.14*sqrt(ea))*((1.35*ifelse((Rs/Rso)>1,1,Rs/Rso))-0.35) ## Eq. 39 Allen & al. 1998.
  # Radiacion neta sobre la superficie del cultivo
  Rn <- Rns-Rnl ## Eq. 40 Allen & al. 1998.
  # Flujo de calor del suelo
  G <-0 ## Eq. 42 Allen & al. 1998. G puede ser ignorado para calculo diarios
  
  # ---- Evapotranspiración ----
  eto <- alpha * (pendiente.curva.presion.vapor/(pendiente.curva.presion.vapor + constante.psicrometrica) * 
      Rn / lambda - G / lambda) # well-watered crop evapotranspiration in a semi-arid and windy location (S5.37)

    return(eto)
    
}

evapotranspiracion.turc <- function(tmax = NULL, tmin = NULL, tmed = NULL,
  HRmed = NULL, HRmax = NULL, HRmin = NULL, helio = NULL, Rs = NULL,
  presion = NULL, albedo = 0.23, a = NULL, b = NULL, z = NULL, lat = NULL,
  k = 0.16, humedo = TRUE, date = NULL) {
  
  # tmax: temperatura máxima diaria [°C]
  # tmin: temperatura mínima diaria [°C]
  # tmed: temperatura media diaria [°C]
  # td: temperatura del punto de rocío [°C]
  # HRmed: humedad relativa media diaria [%]
  # HRmax: humedad relativa diaria máxima [%]
  # HRmin: humedad relativa mínima [%]
  # u2: velocidad media del viento a 2 m de altura [m/s]
  # helio: heliofania, duración del brillo solar [horas]
  # Rs: radiación solar de onda corda [MJ.m^(-2).day^(-1)]
  # presion: presion atmosférica de la estacion [kPa]
  # albedo: albedo de ls superficie. 0.23 por defecto.
  # a y b: coeficientes de la ecuación de Armstrong para la estimación de radiación global
  # z: altura de la estacion [m]
  # lat: latitud de la estacion en grados
  # k: coeficiente de ajuste modelo de radiación de Hargreaves. 0.16 tierra firme, 0.19 zonas costeras [°C^-0.5]
  # date: fecha del día en formato de fecha
  # terms: devuelve los terminos radiativo y aerodinamico por seprado en mm
  
  ##################
  require(dplyr)
  require(sirad)
  require(lubridate)
  require(gtools)
  
  ##################
  
  # Si la temperatura media no está disponible, calcular
  if (gtools::invalid(tmed)) {
    tmed <- (tmin+tmax)/2  ## Eq. 9 Allen & al. 1998.
  }
  
  # Si la humedad relativa media no está disponible, calcular
  if (gtools::invalid(HRmed) & !gtools::invalid(HRmax) | !gtools::invalid(HRmin)) {
    HRmed = (HRmin+HRmax)/2
  }
  
  # Balance de radiación
  # Radiación estraterrestre
  # constante solar = 0.0820 MJ.m^(-2).min^(-1)
  Gsc <- 0.0820
  # convertir latitudes en grados a radianes
  phi <- sirad::radians(lat)           # Eq. 23 Allen & al. 1998.
  # dia juliano correspondiente al dia i
  dia.juliano <- lubridate::yday(date)
  # Inversa de la distancia entre la Tierra y el sol
  Dr <- 1+0.033*cos(2*pi*dia.juliano/365)      # Eq. 23 Allen & al. 1998.
  # inclinacion solar [rad]
  delta <- 0.409*sin((2*pi*dia.juliano/365)-1.39) # Eq. 24 Allen & al. 1998.
  # angulo solar a poniente [rad]
  Ws <- acos(-tan(phi)*tan(delta))        # Eq. 25 Allen & al. 1998.
  
  # Extraterrestrial radiation for daily periods [MJ.m^(-2).day^(-1)]
  eRad <- (24*60/pi)*Gsc*Dr*(Ws*sin(phi)*sin(delta)+cos(phi)*cos(delta)*sin(Ws)) # Eq. 21 Allen & al. 1998.
  
  # Radiación solar o global
  if (!gtools::invalid(Rs)) {
    Rs
  } else if (!gtools::invalid(helio) & gtools::invalid(Rs)) {
    helio.max <- 24/pi*Ws # Eq. 34 Allen & al. 1998.
    Rs <- (a+b*(helio/helio.max))*eRad # Eq. 35 Allen & al. 1998
  } else { (!gtools::invalid(tmax) & !gtools::invalid(tmin) & gtools::invalid(Rs) | gtools::invalid(helio)) 
    Rs <- sqrt(tmax-tmin)*eRad*k # Eq. 50 Allen & al. 1998 
  }
  
 if (humedo && HRmed >= 50) {
   eto <- 0.013 * (23.88 * Rs + 50) * tmed / (tmed + 15) * (1 + (50 - HRmed) / 70) # Turc reference crop evapotranspiration adjusted for non-humid conditions (RH < 50) by Alexandris et al., (S9.11)
   
 } else {
   eto <- 0.013 * (23.88 * Rs + 50) * tmed / (tmed + 15) # reference crop evapotranspiration by Turc (1961) (S9.10)
   
 }
    
    return(eto)
    
}

evapotranspiracion.thornthwaite <- function(tmed = NA, lat = NA, 
  indice.termal, alpha = NA, date = NA) {
  
  # tmed: temperatura media diaria [°C]
  # lat: latitud de la estacion en grados
  # indice.termal: indice termal de Thorthwaite
  # alpha: factor de corrección
  # date: fecha del día en formato de fecha
  
  require(meteor)
  
  # La temperaturea no puede ser menor a 0
  if (tmed < 0) {
    tmed <- 0
  }
  
  # Estimación de evapotranspiración mensual
  eto.mensual <- 16 * (10*tmed/indice.termal)^alpha
  
  # Corrección por la cantidad de días mensuales
  fotoperiodo <- meteor::photoperiod(doy = lubridate::yday(date), latitude = lat)
  # Factor de correción
  k <- fotoperiodo/360
 
  # Evapotranspiración diaria 
  eto <- eto.mensual * k # [mm.dia-1]
  
  return(eto)

}

