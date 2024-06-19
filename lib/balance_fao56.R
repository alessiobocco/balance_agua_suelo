


fao56 <- function(par = NA, wth = NA, startDate = NA, endDate = NA, sol = NA) {
  
  # Fecha de comienzo de la simulacion
  if (is.na(startDate)) {
    startDate = wth %>%
      dplyr::arrange(1) %>%
      dplyr::slice(1) %>%
      dplyr::pull(1)
  }
  
  wth %<>%
    dplyr::filter(`Year-DOY` >= !!startDate)
  
  if (is.na(endDate)) {
    endDate <- wth %>%
      dplyr::slice(nrow(.)) %>%
      dplyr::pull(1)
  }
  
  # Suponiendo que startDate 
  tcurrent <- startDate
  tdelta <- 1  # Representa un cambio de un día
  
  # Inicializar el estado del modelo como una lista
  io <- list()
  io$i <- 0
  io$Kcbini <- par$Kcbini
  io$Kcbmid <- par$Kcbmid
  io$Kcbend <- par$Kcbend
  io$Lini <- par$Lini
  io$Ldev <- par$Ldev
  io$Lmid <- par$Lmid
  io$Lend <- par$Lend
  io$hini <- par$hini
  io$hmax <- par$hmax
  io$thetaFC <- par$thetaFC
  io$thetaWP <- par$thetaWP
  io$theta0 <- par$theta0
  io$Zrini <- par$Zrini
  io$Zrmax <- par$Zrmax
  io$pbase <- par$pbase
  io$Ze <- par$Ze
  io$REW <- par$REW
  
  # Calcular Total evaporable water (TEW, mm) - FAO-56 Equation 73
  io$TEW <- 1000 * (io$thetaFC - 0.50 * io$thetaWP) * io$Ze
  
  # Calcular Initial depth of evaporation (De, mm) - FAO-56 page 153
  io$De <- io$TEW  # Mismo cálculo que para TEW
  
  # Suponiendo que 'sol' es una lista o NULL en R representando el suelo
  if (is.na(sol)) {
    io$solmthd <- 'D'  # Método por defecto para suelo homogéneo
    
    # Inicialización de parámetros para suelo homogéneo
    io$Dr <- 1000 * (io$thetaFC - io$theta0) * io$Zrini
    io$Drmax <- 1000 * (io$thetaFC - io$theta0) * io$Zrmax
    io$TAW <- 1000 * (io$thetaFC - io$thetaWP) * io$Zrini
    
    # Valores por defecto para variables no consideradas por FAO-56
    io$TAWrmax <- -99.999
    io$Db <- -99.999
    io$TAWb <- -99.999
  } else {
    io$solmthd <- 'L'  # Método para perfil de suelo en capas
    
    # Asignar profundidades de capas y contenido de agua del suelo desde 'sol'
    io$lyr_dpths <- sol$sdata$index
    io$lyr_thFC <- sol$sdata$thetaFC
    io$lyr_thWP <- sol$sdata$thetaWP
    io$lyr_th0 <- sol$sdata$theta0
    
    # Inicializar variables
    io$Dr <- 0
    io$Drmax <- 0
    io$TAW <- 0
    io$TAWrmax <- 0
    
    # Iterar a través del perfil del suelo
    for (dpthmm in 1:(io$lyr_dpths[length(io$lyr_dpths)] * 10 + 1)) {
      # Encontrar el índice de la capa que contiene dpthmm
      lyr_idx <- which(dpthmm <= io$lyr_dpths * 10)[1]
      
      # Cálculos para Dr, Drmax, TAW, TAWrmax
      if (dpthmm <= io$Zrini * 1000) {
        diff <- (io$lyr_thFC[lyr_idx] - io$lyr_th0[lyr_idx])
        io$Dr <- io$Dr + diff
      }
      if (dpthmm <= io$Zrmax * 1000) {
        diff <- (io$lyr_thFC[lyr_idx] - io$lyr_th0[lyr_idx])
        io$Drmax <- io$Drmax + diff
        diff <- (io$lyr_thFC[lyr_idx] - io$lyr_thWP[lyr_idx])
        io$TAWrmax <- io$TAWrmax + diff
      }
    }
    
    # Cálculos para Db y TAWb
    io$Db <- io$Drmax - io$Dr
    io$TAWb <- io$TAWrmax - io$TAW
  }
  
  # Asignación de valores a campos de 'io'
  io$h <- io$hini
  io$Zr <- io$Zrini
  io$fw <- 1.0
  
  # Suponiendo que 'wth' es otro objeto en R con campos 'wndht' y 'rfcrp'
  io$wndht <- 3
  io$rfcrp <- "S"
  io$cons_p <- FALSE
  
  # Crear un vector con los nombres de las columnas
  cnames <- c('Year', 'DOY', 'DOW', 'Date', 'ETref', 'tKcb', 'Kcb',
              'h', 'Kcmax', 'fc', 'fw', 'few', 'De', 'Kr', 'Ke', 'E',
              'DPe', 'Kc', 'ETc', 'TAW', 'TAWrmax', 'TAWb', 'Zr', 'p',
              'RAW', 'Ks', 'Kcadj', 'ETcadj', 'T', 'DP', 'Dinc', 'Dr',
              'fDr', 'Drmax', 'fDrmax', 'Db', 'fDb',  'Rain')
  
  # Crear un data frame vacío con estas columnas
  odata <- data.frame(matrix(ncol = length(cnames), nrow = 0))
  colnames(odata) <- cnames
  
  while (tcurrent <= endDate) {
    
    # Actualizar el objeto ModelState (io) con los datos climáticos
    io$ETref <- wth %>%
      dplyr::filter(`Year-DOY` == !!tcurrent) %>%
      dplyr::pull(ETref)
      
    io$rain <- wth %>%
      dplyr::filter(`Year-DOY` == !!tcurrent) %>%
      dplyr::pull(Rain)
    
    io$wndsp <- if (exists("Wndsp", where = wth) && 
                    nrow(filter(wth, `Year-DOY` == tcurrent)) > 0) {
      wth %>%
        filter(`Year-DOY` == tcurrent) %>%
        pull(Wndsp)
    } else {
      NA  # 
    }

    if (is.na(io$wndsp)) {
      io$wndsp <- 2.0
    }
    
    io$rhmin <- if (exists("RHmin", where = wth) && 
                    nrow(filter(wth, `Year-DOY` == tcurrent)) > 0) {
      wth %>%
        filter(`Year-DOY` == tcurrent) %>%
        pull(RHmin)
    } else {
      NA  # 
    }
    
    if (is.na(io$rhmin)) {
      tmax <- wth %>%
        dplyr::filter(`Year-DOY` == !!tcurrent) %>%
        dplyr::pull(Tmax)
      tmin <- wth %>%
        dplyr::filter(`Year-DOY` == !!tcurrent) %>%
        dplyr::pull(Tmin)
      tdew <-    io$rhmin <- if (exists("Tdew", where = wth) && 
                                 nrow(filter(wth, `Year-DOY` == tcurrent)) > 0) {
        wth %>%
          filter(`Year-DOY` == tcurrent) %>%
          pull(Tdew)
      } else {
        tmin
      }
      # Cálculos según ASCE (2005) Eqs. 7 y 8
      emax <- 0.6108 * exp((17.27 * tmax) / (tmax + 237.3))
      ea <- 0.6108 * exp((17.27 * tdew) / (tdew + 237.3))
      io$rhmin <- ea / emax * 100
    }
    
    if (is.na(io$rhmin)) {
      io$rhmin <- 45
    }
    
    # Suponiendo que 'odata' es un data frame y 'tcurrent' es la fecha actual en R
    # Convertir 'tcurrent' a componentes de fecha
    
    year <- lubridate::year(as.Date(tcurrent, format = "%Y-%j"))
    doy <- lubridate::yday(as.Date(tcurrent, format = "%Y-%j"))  # Día del año
    dow <- lubridate::wday(as.Date(tcurrent, format = "%Y-%j"))  # Día de la semana
    dat <- lubridate::date(as.Date(tcurrent, format = "%Y-%j"))  # Fecha mm/dd/yy
    
    # Avanzar al dia siguiente de la siulacion
    io <- advance(io)
    
    # Crear una fila de datos
    data <- data.frame(Year = year, DOY = doy, DOW = dow, Date = dat, ETref = io$ETref, tKcb = io$tKcb, 
                       Kcb = io$Kcb, h = io$h, Kcmax = io$Kcmax, fc = io$fc, fw = io$fw, 
                       few = io$few, De = io$De, Kr = io$Kr, Ke = io$Ke, E = io$E, 
                       DPe = io$DPe, Kc = io$Kc, ETc = io$ETc, TAW = io$TAW, TAWrmax = io$TAWrmax, 
                       TAWb = io$TAWb, Zr = io$Zr, p = io$p, RAW = io$RAW, Ks = io$Ks, 
                       Kcadj = io$Kcadj, ETcadj = io$ETcadj, `T` = io$Transpiration, DP = io$DP, Dinc = io$Dinc, 
                       Dr = io$Dr, fDr = io$fDr, Drmax = io$Drmax, fDrmax = io$fDrmax, 
                       Db = io$Db, fDb = io$fDb, Rain = io$rain)
    
    # Añadir la fila al data frame 'odata'
    # 'mykey' puede ser una combinación de año y día del año para crear un identificador único
    mykey <- paste(year, doy, sep = "-")
    odata[mykey, ] <- data
    
    # Avanzar al siguiente día
    tcurrent <- lubridate::date(as.Date(tcurrent, format = "%Y-%j")) + lubridate::day(1)
    tcurrent <- paste0(lubridate::year(tcurrent), "-", lubridate::yday(tcurrent))
    io$i <- io$i + 1
    
  }
  
  return(odata)
}

  
advance <- function(io) {
  # Basal crop coefficient (Kcb) - De Tablas FAO-56 11 y 17
  s1 <- io$Lini
  s2 <- s1 + io$Ldev
  s3 <- s2 + io$Lmid
  s4 <- s3 + io$Lend
   if (0 <= io$i && io$i <= s1) {
    io$tKcb <- io$Kcbini
    io$Kcb <- io$Kcbini
  } else if (s1 < io$i && io$i <= s2) {
    io$tKcb <- io$tKcb + (io$Kcbmid - io$Kcbini) / (s2 - s1)
    io$Kcb <- io$Kcb + (io$Kcbmid - io$Kcbini) / (s2 - s1)
  } else if (s2 < io$i && io$i <= s3) {
    io$tKcb <- io$Kcbmid
    io$Kcb <- io$Kcbmid
  } else if (s3 < io$i && io$i <= s4) {
    io$tKcb <- io$tKcb + (io$Kcbmid - io$Kcbend) / (s3 - s4)
    io$Kcb <- io$Kcb + (io$Kcbmid - io$Kcbend) / (s3 - s4)
  } else if (s4 < io$i) {
    io$tKcb <- io$Kcbend
    io$Kcb <- io$Kcbend
  }

  # Altura de la planta (h, m)
  io$h <- max(io$hini + (io$hmax - io$hini) * (io$Kcb - io$Kcbini) / (io$Kcbmid - io$Kcbini), 0.001, io$h)
  
  # Profundidad de raíz (Zr, m) - FAO-56 página 279
  io$Zr <- max(io$Zrini + (io$Zrmax - io$Zrini) * (io$tKcb - io$Kcbini) / (io$Kcbmid - io$Kcbini), 0.001, io$Zr)
  
  # Coeficiente máximo del cultivo (Kcmax) - FAO-56 Ecuación 72
  u2 <- io$wndsp * (4.87 / log(67.8 * io$wndht - 5.42))
  u2 <- sort(c(1.0, u2, 6.0))[2]
  rhmin <- sort(c(20.0, io$rhmin, 80.0))[2]
  # Tipo de cultivo: alto o bajo
  if (io$rfcrp == 'S') {
    io$Kcmax <- max(1.2 + (0.04 * (u2 - 2.0) - 0.004 * (rhmin - 45.0)) * (io$h / 3.0)^0.3, io$Kcb + 0.05)
  } else if (io$rfcrp == 'T') {
    io$Kcmax <- max(1.0, io$Kcb + 0.05)
  }
  
  # Fracción de cobertura del dosel (fc, 0.0-0.99) - FAO-56 Ecuación 76
  io$fc <- min(max(0.0, ((io$Kcb - io$Kcbini) / (io$Kcmax - io$Kcbini)) ^ (1.0 + 0.5 * io$h)), 0.99)
  
  # Fracción de la superficie del suelo mojada (fw) - FAO-56 Tabla 20, página 149
  if (io$rain > 0.0) {
    # fw=fw input
  } else if (io$rain <= 0.0) {
    # fw=fw input
  } else if (io$rain >= 3.0) {
    io$fw <- 1.0
  } else {
    # fw = previous fw
  }
  
  # Fracción expuesta y mojada del suelo (few, 0.01-1.0) - FAO-56 Ecuación 75
  io$few <- min(max(0.01, min(1.0 - io$fc, io$fw)), 1.0)
  
  # Coeficiente de reducción de la evaporación (Kr, 0-1) - FAO-56 Ecuación 74
  io$Kr <- pmin(pmax(0.0, (io$TEW - io$De) / (io$TEW - io$REW)), 1.0)
  
  # Coeficiente de evaporación (Ke) - FAO-56 Ecuación 71
  io$Ke <- min(io$Kr * (io$Kcmax - io$Kcb), io$few * io$Kcmax)
  
  # Evaporación del agua del suelo (E, mm) - FAO-56 Ecuación 69
  io$E <- io$Ke * io$ETref
  
  # Percolación profunda bajo suelo expuesto (DPe, mm) - FAO-56 Ecuación 79
  runoff <- 0.0
  io$DPe <- max(io$rain - runoff - io$De, 0.0)
  
  # Profundidad acumulada de evaporación (De, mm) - FAO-56 Ecuaciones 77 y 78
  De <- io$De - (io$rain - runoff)   / io$fw + io$E / io$few + io$DPe
  io$De <- pmin(pmax(0.0, De), io$TEW)
  
  # Coeficiente del cultivo (Kc) - FAO-56 Ecuación 69
  io$Kc <- io$Ke + io$Kcb
  
  # Evapotranspiración del cultivo no estresado (ETc, mm) - FAO-56 Ecuación 69
  io$ETc <- io$Kc * io$ETref
  
  # Decidir el método de cálculo de Agua Disponible Total (TAW)
  if (io$solmthd == 'D') {
    # Agua Disponible Total (TAW, mm) - FAO-56 Ecuación 82
    io$TAW <- 1000.0 * (io$thetaFC - io$thetaWP) * io$Zr
  } else if (io$solmthd == 'L') {
    io$TAW <- 0.0
    # Iterar a través del perfil del suelo en incrementos de 1 mm
    for (dpthmm in 1:(max(io$lyr_dpths) * 10 + 1)) {
      # Encontrar el índice de la capa del suelo que contiene dpthmm
      lyr_idx <- which(dpthmm <= io$lyr_dpths * 10)[1]
      # Calcular el Agua Disponible Total (TAW, mm)
      if (dpthmm <= io$Zr * 1000) {
        diff <- (io$lyr_thFC[lyr_idx] - io$lyr_thWP[lyr_idx])
        io$TAW <- io$TAW + diff
      }
    }
    # Agua Disponible Total en la capa inferior (TAWb, mm)
    io$TAWb_prev <- io$TAWb
    io$TAWb <- io$TAWrmax - io$TAW
  }
  
  # Fracción de TAW agotada (p, 0.1-0.8) - FAO-56 p162 y Tabla 22
  if (io$cons_p) {
    io$p <- io$pbase
  } else {
    io$p <- min(max(0.1, io$pbase + 0.04 * (5.0 - io$ETc)), 0.8)
  }
  
  # Agua Disponible Fácilmente (RAW, mm) - FAO-56 Ecuación 83
  io$RAW <- io$p * io$TAW
  
  # Factor de reducción de la transpiración (Ks, 0.0-1.0) - FAO-56 Ecuación 84
  io$Ks <- pmin(pmax(0.0, (io$TAW - io$Dr) / (io$TAW - io$RAW)), 1.0)
  
  # Coeficiente de cultivo ajustado (Kcadj) - FAO-56 Ecuación 80
  io$Kcadj <- io$Ks * io$Kcb + io$Ke
  
  # Evapotranspiración del cultivo ajustada (ETcadj, mm) - FAO-56 Ecuación 80
  io$ETcadj <- io$Kcadj * io$ETref
  
  # Transpiración del cultivo ajustada (T, mm)
  io$Transpiration <- (io$Ks * io$Kcb) * io$ETref
  
  # Métodos de balance hídrico
  if (io$solmthd == 'D') {
    # Percolación profunda (DP, mm) - FAO-56 Ecuación 88
    # La capa límite se considera en la profundidad de la zona radicular (Zr)
    io$DP <- pmax(io$rain - runoff  - io$ETcadj - io$Dr, 0.0)
    
    # Agotamiento del agua del suelo en la zona radicular (Dr, mm) - FAO-56 Ecuaciones 85 y 86
    Dr <- io$Dr - (io$rain - runoff)  + io$ETcadj + io$DP
    io$Dr <- pmin(pmax(0.0, Dr), io$TAW)
    
    # Fracción de agotamiento del agua del suelo en la zona radicular (fDr, mm/mm)
    io$fDr <- 1.0 - ((io$TAW - io$Dr) / io$TAW)
    
    # FAO-56 no considera por defecto las siguientes variables
    io$Dinc <- -99.999
    io$Drmax <- -99.999
    io$fDrmax <- -99.999
    io$Db <- -99.999
    io$fDb <- -99.999
  } else if (io$solmthd == 'L') {
    # Percolación profunda (DP, mm)
    # La capa límite está en la profundidad máxima de las raíces (Zrmax)
    io$DP <- pmax(io$rain - runoff - io$ETcadj - io$Drmax, 0.0)
    
    # Incremento de agotamiento debido al crecimiento de las raíces (Dinc, mm)
    # Calculado a partir de Db basado en el cambio incremental en TAWb
    if (io$TAWb_prev > 0.0) {
      io$Dinc <- io$Db * (1.0 - (io$TAWb / io$TAWb_prev))
    } else {
      io$Dinc <- 0.0
    }
    
    # Agotamiento del agua del suelo en la zona radicular (Dr, mm)
    Dr <- io$Dr - (io$rain - runoff) + io$ETcadj + io$Dinc
    io$Dr <- pmin(pmax(0.0, Dr), io$TAW)
    
    # Fracción de agotamiento del agua del suelo en la zona radicular (fDr, mm/mm)
    io$fDr <- 1.0 - ((io$TAW - io$Dr) / io$TAW)
    
    # Agotamiento del agua del suelo en la profundidad máxima de raíz (Drmax, mm)
    Drmax <- io$Drmax - (io$rain - runoff) + io$ETcadj + io$DP
    io$Drmax <- pmin(pax(0.0, Drmax), io$TAWrmax)
    
    # Fracción de agotamiento del agua del suelo en Zrmax (fDrmax, mm/mm)
    io$fDrmax <- 1.0 - ((io$TAWrmax - io$Drmax) / io$TAWrmax)
    
    # Agotamiento del agua del suelo en la capa inferior (Db, mm)
    Db <- io$Drmax - io$Dr
    io$Db <- pmin(pmax(0.0, Db), io$TAWb)
    
    # Fracción de agotamiento del agua del suelo en la capa inferior (fDb, mm/mm)
    if (io$TAWb > 0.0) {
      io$fDb <- 1.0 - ((io$TAWb - io$Db) / io$TAWb)
    } else {
      io$fDb <- 0.0
    }
  }
  
  return(io)
}







