DAYLEN <- function(DOY, XLAT) {
  PI <- pi
  RAD <- PI / 180
  
  # Calcula la declinación solar
  DEC <- -23.45 * cos(2 * PI * (DOY + 10) / 365)
  
  # Limita SOC para latitudes por encima de los círculos polares
  SOC <- tan(RAD * DEC) * tan(RAD * XLAT)
  SOC <- min(max(SOC, -1), 1)
  
  # Calcula la duración del día, amanecer y atardecer
  DAYL <- 12 + 24 * asin(SOC) / PI
  SNUP <- 12 - DAYL / 2
  SNDN <- 12 + DAYL / 2
  
  return(list(DAYL = DAYL, DEC = DEC, SNDN = SNDN, SNUP = SNUP))
}

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
  SNUP <- 12 - TWILEN / 2
  SNDN <- 12 + TWILEN / 2
  
  return(list(TWILEN = TWILEN, SNUP = SNUP, SNDN = SNDN))
}

HMET <- function(DAYL, SNDN, SNUP, TMAX, TMIN, XLAT) {
  # Inicialización
  TAVG <- 0.0
  TDAY <- 0.0
  NDAY <- 0
  TS <- 24 # Number of hourly time steps per day
  TINCR <- 24.0 / TS
  
  # Vectores de salida
  TAIRHR <- numeric(TS)
  TGRO <- numeric(TS)

  # Bucle para calcular los datos meteorológicos por hora
  for (H in 1:TS) {
    HS <- H * TINCR
    
    TAIRHR[H] <- HTEMP(DAYL, HS, SNDN, SNUP, TMAX, TMIN)
    
    TAVG <- TAVG + TAIRHR[H]
    
    if (H >= SNUP && H <= SNDN) {
      TDAY <- TDAY + TAIRHR[H]
      NDAY <- NDAY + 1
    }
  }
  
  TAVG <- TAVG / TS
  TDAY <- TDAY / NDAY
  TGRODY <- TDAY
  TGROAV <- TAVG
  TGRO <- TAIRHR
  
  return(list(
    TAIRHR = TAIRHR,
    TAVG = TAVG,
    TDAY = TDAY,
    TGRO = TGRO,
    TGROAV = TGROAV,
    TGRODY = TGRODY
  ))
}

HTEMP <- function(DAYL, HS, SNDN, SNUP, TMAX, TMIN) {
  # Parámetros definidos en Parton y Logan
  A <- 2.0
  B <- 2.2
  C <- 1.0
  PI <- pi
  
  # Calcular las horas para temperatura mínima y máxima
  MIN <- SNUP + C
  MAX <- MIN + DAYL / 2 + A
  
  # Calcular la temperatura al atardecer y la temperatura mínima teórica
  T <- 0.5 * PI * (SNDN - MIN) / (MAX - MIN)
  TSNDN <- TMIN + (TMAX - TMIN) * sin(T)
  TMINI <- (TMIN - TSNDN * exp(-B)) / (1 - exp(-B))
  
  # Calcular el tiempo de decaimiento exponencial
  HDECAY <- 24.0 + C - DAYL
  
  # Calcular la temperatura para horas diurnas y nocturnas
  if (HS >= SNUP + C && HS <= SNDN) {
    T <- 0.5 * PI * (HS - MIN) / (MAX - MIN)
    TAIRHR <- TMIN + (TMAX - TMIN) * sin(T)
  } else {
    if (HS < SNUP + C) {
      T <- 24.0 + HS - SNDN
    } else if (HS > SNDN) {
      T <- HS - SNDN
    }
    ARG <- -B * T / HDECAY
    TAIRHR <- TMINI + (TSNDN - TMINI) * exp(ARG)
  }
  
  return(TAIRHR)
}


calculate_development_rates <- function(NPHS, TS = 24, DAS, NR1, TGRO, 
                                       TSELC, CTMP, TB, TO1, TO2, TM,
                                       FT, FUDAY, DLTYP, CSDVAR, CLDVAR, 
                                       THVAR, CSDVRR, CLDVRR, DAYL, PSENP) {
  
  for (J in 2:NPHS) {
    K <- TSELC[J]
    FT[J] <- 0
    for (I in 1:TS) {
      FTHR <- CURV(CTMP[J], TB[K], TO1[K], TO2[K], TM[K], TGRO[I])
      FT[J] <- FT[J] + FTHR / TS
    }
    
    if (DAS < NR1) {
      FUDAY[J] <- CURV(DLTYP[J], 1.0, CSDVAR, CLDVAR, THVAR, DAYL)
    } else {
      FUDAY[J] <- CURV(DLTYP[J], 1.0, CSDVRR, CLDVRR, THVAR, DAYL)
    }
    
  }
  
  return(list(FT = FT, FUDAY = FUDAY, FSW = FSW, FNSTR = FNSTR, FPSTR = FPSTR))
}








HRAD <- function(BETA, HS, ISINB, SNDN, SNUP, SRAD) {
  # Constantes
  PI <- pi
  RAD <- PI / 180.0
  
  # Inicializar la radiación solar por hora
  RADHR <- 0.0
  
  # Cálculos para horas diurnas
  if (HS > SNUP && HS < SNDN) {
    # Calcula el seno de BETA (elevación solar)
    SINB <- sin(RAD * BETA)
    
    # Calcula la radiación solar instantánea usando la ecuación de Spitters
    RADHR <- SINB * (1.0 + 0.4 * SINB) * SRAD * 1.0E6 / ISINB
  }
  
  return(RADHR)
}

FRACD <- function(BETA, CLOUDS, HS, RADHR, S0N, SNDN, SNUP) {
  # Constantes
  PI <- pi
  RAD <- PI / 180
  
  # Inicializar variables de salida
  FRDIFP <- 1.0
  FRDIFR <- 1.0
  AMTRH <- 0.0
  
  # Cálculos para horas diurnas
  if (HS > SNUP && HS < SNDN) {
    SINB <- sin(RAD * BETA)
    COSB <- cos(RAD * BETA)
    COS90B <- cos(RAD * (90.0 - BETA))
    
    # Calcular la transmisión atmosférica instantánea
    S0 <- S0N * SINB
    if (S0 > 0) {
      AMTRH <- RADHR / S0
    }
    
    # Calcular la fracción de difusión horaria
    FRDFH <- 0.156 + (0.86 / (1.0 + exp(11.1 * (AMTRH - 0.53))))
    FRDFH <- min(FRDFH, 1.0)
    
    # Calcular la irradiación global difusa instantánea
    SDF <- RADHR * FRDFH
    
    # Disminuir SDF por la parte circumsolar de la radiación difusa
    CORR <- (1.0 - CLOUDS) * COS90B^2 * COSB^3
    SDF <- SDF / (1.0 + CORR)
    
    # Calcular la fracción difusa instantánea de global y PAR
    if (RADHR > 0.0) {
      FRDIFR <- SDF / RADHR
    } else {
      FRDIFR <- 0.0
    }
    
    FRDIFP <- (1.0 + 0.3 * (1.0 - CLOUDS)) * FRDIFR
    FRDIFR <- min(max(FRDIFR, 0.0), 1.0)
    FRDIFP <- min(max(FRDIFP, 0.0), 1.0)
  }
  
  return(list(FRDIFP = FRDIFP, FRDIFR = FRDIFR, AMTRH = AMTRH))
}

HPAR <- function(HS, PAR, RADHR, SNDN, SNUP, SRAD) {
  # Parámetros
  PARQC <- 4.6  # umol/J
  
  # Inicializar la variable de salida
  PARHR <- numeric(length(HS))
  
  for (i in seq_along(HS)) {
    h <- HS[i]
    radhr <- RADHR[i]
    
    # Cálculos para horas diurnas
    if (h > SNUP && h < SNDN) {
      if (PAR > 1e-4) {
        PARHR[i] <- radhr * PAR / SRAD
      } else {
        # PARFAC según Lizaso et al. (2003)
        PARFAC <- (0.43 + 0.12 * exp(-SRAD / 2.8)) * PARQC
        PARHR[i] <- radhr * PARFAC
      }
    } else {
      # Noche
      PARHR[i] <- 0.0
    }
  }
  
  return(PARHR)
}

HWIND <- function(DAYL, HS, SNDN, SNUP, WINDAV) {
  # Parámetros iniciales
  PI <- 3.14159
  WNIGHT <- 0.3  # Porcentaje de velocidad del viento en la noche respecto al máximo
  
  # Calcular la velocidad máxima del viento del día
  WMAX <- WINDAV / (WNIGHT + (1.0 - WNIGHT) * DAYL / 48.0)
  
  # Inicializar la velocidad del viento por hora
  WINDHR <- numeric(length(HS))
  
  for (i in seq_along(HS)) {
    h <- HS[i]
    
    # Cálculos para horas diurnas
    if (h > SNUP && h < SNDN) {
      T <- 2.0 * PI * (h - SNUP) / (SNDN - SNUP) + 1.5 * PI
      WINDHR[i] <- WMAX * WNIGHT + (1.0 - WNIGHT) * WMAX / 2.0 * (1.0 + sin(T))
    } else {
      # Noche
      WINDHR[i] <- WMAX * WNIGHT
    }
  }
  
  return(WINDHR)
}

VPSAT <- function(T) {
  610.78 * exp(17.269 * T / (T + 237.3))
}

