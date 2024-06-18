
#CURV("LIN", TB[1], TO1[1], TO2[1], TM[1], 11.402)

CURV <- function(CTYPE, XB, X1, X2, XM, X) {
  CURV <- 1.0
  if (CTYPE == 'NON' || CTYPE == 'non') {
    return(CURV)
  }
  
  if (CTYPE == 'LIN' || CTYPE == 'lin') {
    CURV <- 0.0
    if (X > XB && X < X1) CURV <- (X - XB) / (X1 - XB)
    if (X >= X1 && X <= X2) CURV <- 1.0
    if (X > X2 && X < XM) CURV <- 1.0 - (X - X2) / (XM - X2)
    CURV <- max(CURV, 0.0)
    CURV <- min(CURV, 1.0)
  }
  
  if (CTYPE == 'QDR' || CTYPE == 'qdr') {
    CURV <- 0.0
    if (X > XB && X < X1) CURV <- 1.0 - ((X1 - X) / (X1 - XB))^2
    if (X >= X1 && X <= X2) CURV <- 1.0
    if (X > X2 && X < XM) CURV <- 1.0 - ((X - X2) / (XM - X2))^2
    CURV <- max(CURV, 0.0)
    CURV <- min(CURV, 1.0)
  }
  
  if (CTYPE == 'INL' || CTYPE == 'inl') {
    CURV <- 1.0
    if (X > X1 && X < X2) CURV <- 1.0 - (1.0 - XM) * ((X - X1) / (X2 - X1))
    if (X >= X2) CURV <- XM
    CURV <- max(CURV, XM)
    CURV <- min(CURV, 1.0)
  }
  
  if (CTYPE == 'SHO' || CTYPE == 'sho') {
    if (X <= X1) {
      CURV <- 1.0
    } else if (X > X1 && X < X2) {
      CURV <- 1.0 - (1.0 - XM) * ((X - X1) / (X2 - X1))
    } else if (X >= X2) {
      CURV <- XM
    }
    CURV <- max(CURV, XM)
    CURV <- min(CURV, 1.0)
  }
  
  if (CTYPE == 'LON' || CTYPE == 'lon') {
    if (X < X2) {
      CURV <- XM
    } else if (X >= X2 && X < X1) {
      CURV <- 1.0 - (1.0 - XM) * ((X1 - X) / (X1 - X2))
    } else {
      CURV <- 1.0
    }
    CURV <- max(CURV, XM)
    CURV <- min(CURV, 1.0)
  }
  
  # SIN (Sinusoidal Curve): Se utiliza una curva sinusoidal que varía entre 0 y 1, 
  # adaptando los parámetros para ajustarse a las fases de la curva en función de X.
  if (CTYPE == 'SIN' || CTYPE == 'sin') {
    CURV <- 0.0
    if (X > XB && X < X1) {
      CURV <- 0.5 * (1 + cos(2 * pi * (X - X1) / (2 * (X1 - XB))))
    }
    if (X >= X1 && X <= X2) {
      CURV <- 1.0
    }
    if (X > X2 && X < XM) {
      CURV <- 0.5 * (1 + cos(2 * pi * (X2 - X) / (2 * (XM - X2))))
    }
    CURV <- max(CURV, 0.0)
    CURV <- min(CURV, 1.0)
  }
  
  # REV (Reversible Process): Representa un proceso reversible, donde la curva 
  # puede descender y luego ascender, con una normalización basada en los parámetros 
  # proporcionados.
  if (CTYPE == 'REV' || CTYPE == 'rev') {
    CURV <- 1.0
    if (X > XB && X < X1) {
      CURV <- 1.0 - ((X - XB) / (X1 - XB))
    }
    if (X >= X1 && X <= X2) {
      CURV <- 0.0 - ((X - X1) / (X2 - X1))
    }
    if (X > X2) {
      CURV <- -1.0
    }
    CURV <- max(CURV, -1.0)
    CURV <- min(CURV, 1.0)
    CURV <- CURV * XM
  }
  
  # DHD (Cold Dehardening): Representa el proceso de dehardening 
  # (pérdida de resistencia al frío) con una tasa que aumenta con la temperatura.
  if (CTYPE == 'DHD' || CTYPE == 'dhd') {
    CURV <- 0.0
    if (X > XB && X < X1) {
      CURV <- (X - XB) / (X1 - XB)
    }
    if (X >= X1 && X <= X2) {
      CURV <- 1.0
    }
    if (X > X2) {
      CURV <- 1.0
    }
    CURV <- max(CURV, 0.0)
    CURV <- min(CURV, 1.0)
    CURV <- CURV * XM
  }
  
  # DRD (Dormancy Related Decrease): Esta curva se utiliza para reducir las 
  # tasas de procesos a medida que avanza la dormancia. La tasa de reducción 
  # depende de la longitud del día, disminuyendo a medida que esta aumenta y 
  # tiene el máximo efecto en la dormancia completa (XM).
  if (CTYPE == 'DRD' || CTYPE == 'drd') {
    CURV <- X2
    if (X > XB && X < X1) {
      CURV <- X2 + (XM - X2) * (X - XB) / (X1 - XB)
    }
    if (X >= X1) {
      CURV <- XM
    }
    CURV <- max(CURV, X2)
    CURV <- min(CURV, XM)
  }
  
  # CDD (Curvilinear Dormancy Related Decrease): Similar a DRD, pero usa una 
  # función curvilínea en lugar de una función lineal para la reducción. 
  # La tasa de reducción también depende de la longitud del día.
  if (CTYPE == 'CDD' || CTYPE == 'cdd') {
    CURV <- X2
    if (X > XB && X < X1) {
      CURV <- XM - ((XM - X2) * ((X1 - X) / (X1 - XB))^2)
    }
    if (X >= X1) {
      CURV <- XM
    }
    CURV <- max(CURV, X2)
    CURV <- min(CURV, XM)
  }
  
  # EXK (Exponential Curve with k): Representa una función exponencial donde XB 
  # establece la amplitud de la curva, X1/XM determina la curvatura, y 
  # X2 desplaza la curva a lo largo del eje X.
  if (CTYPE == 'EXK' || CTYPE == 'exk') {
    CURV <- XB - exp(X1 * (X - X2) / XM)
  }
  
  # VOP (Variable Order Polynomial): Representa un polinomio de orden variable, 
  # donde X2 determina el orden del polinomio. La función escala a 0 por debajo 
  # de XB y por encima de XM.
  if (CTYPE == 'VOP' || CTYPE == 'vop') {
    CURV <- 0.0
    if (X > XB && X < XM) {
      CURV <- ((X - XB)^X2 * (XM - X)) / ((X1 - XB)^X2 * (XM - X1))
    }
    if (X >= XM) {
      CURV <- 0.0
    }
    CURV <- max(CURV, 0.0)
  }
  
  # Q10: Utiliza la función Q10, que es una forma de cuantificar cómo cambia 
  # la velocidad de una reacción química con un aumento de temperatura de 10 grados Celsius. 
  # XB es la temperatura de referencia, y X2 es el factor de incremento Q10.
  if (CTYPE == 'Q10' || CTYPE == 'q10') {
    CURV <- X1 * (X2^((X - XB) / 10))
  }
  
  # PWR (Power): Eleva X a la potencia especificada por X1, con XB como 
  # multiplicador para la función principal y X2 como un multiplicador de escala.
  if (CTYPE == 'PWR' || CTYPE == 'pwr') {
    if (X < 0.0) {
      CURV <- X2 * XB * (0^X1)
    } else {
      CURV <- X2 * XB * (X^X1)
    }
  }
  
  return(CURV)
}




soilt <- function(ALBEDO, B, CUMDPT, DOY, DP, HDAY, NLAYR, SRAD, TAMP, TAV = NA, TAVG, TMAX, WW = 0.20, PESW, DSMID, TMA = NA, ATOT = NA) {
  # Inicializar variables
  if (any(is.na(TMA))) {
    TMA <- numeric(5) # Asumir que TMA se pasa o se define de alguna manera anteriormente
  }
  
  if (is.na(TAV)) {
    TAV = 20
  }
  
  if (is.na(TAVG)) {
    TAVG = 20
  }
  
  if (is.na(ATOT)) {
    ATOT = TAV * 5
  }
  
  # Cálculos preliminares
  ALX <- (DOY - HDAY) * 0.0174
  ATOT <- ATOT - TMA[5]
  
  # Actualizar TMA
  for (K in 5:2) {
    TMA[K] <- TMA[K - 1]
  }
  
  TMA[1] <- (1 - ALBEDO) * (TAVG + (TMAX - TAVG) * sqrt(SRAD * 0.03)) + ALBEDO * TMA[1]
  TMA[1] <- round(TMA[1] * 10000) / 10000
  ATOT <- ATOT + TMA[1]
  
  # Calcular contenido de agua WC
  WC <- max(0.01, PESW) / (WW * CUMDPT) * 10
  
  # Más cálculos
  FX <- exp(B * ((1 - WC) / (1 + WC))^2)
  DD <- FX * DP
  TA <- TAV + TAMP * cos(ALX) / 2
  DT <- ATOT / 5 - TA
  
  # Calcular temperatura del suelo por capa
  ST <- numeric(NLAYR)
  for (L in 1:NLAYR) {
    ZD <- -DSMID[L] / DD
    ST[L] <- TAV + (TAMP / 2 * cos(ALX + ZD) + DT) * exp(ZD)
    ST[L] <- round(ST[L] * 1000) / 1000
  }
  
  # Calcular temperatura de la superficie del suelo
  SRFTEMP <- TAV + (TAMP / 2 * cos(ALX) + DT)
  
  # Retornar resultados
  return(list(ST = ST, SRFTEMP = SRFTEMP, TMA = TMA))
}

# Ejemplo de uso de la función
# result <- soilt(ALBEDO, B, CUMDPT, DOY, DP, HDAY, NLAYR, PESW, SRAD, TAMP, TAV = NA, TAVG, TMAX, WW, DSMID)
# print(result)


stemp <- function(SRAD, TAVG, TMAX, XLAT, TAV = NA, TMA = NA, ATOT = NA, TAMP, DOY,
                  NLAYR, DLAYR, BD, DUL, LL, DPLAYR, DSMID, MSALB ) {

  # Establecer HDAY basado en la latitud
  HDAY <- if (XLAT < 0) 20 else 200 # DOY (hottest) for southern/northern hemisphere
  ALX <- (as.numeric(DOY) - HDAY) * 0.0174

  TBD <- sum(BD * DLAYR)
  TDL <- sum(DUL * DLAYR)
  TLL <- sum(LL * DLAYR)

  ABD <- TBD / DS[NLAYR]
  FX <- exp(ABD / (ABD + 686 * exp(-5.63 * ABD)))
  DP <- 1000 + 2500 * FX
  WW <- 0.356 - 0.144 * ABD
  B <- log(500 / DP)
  ALBEDO <- MSALB
  
  PESW <- max(0, TDL - TLL)
  
  # Llamada a otra función que debería estar definida en R similar a SOILT en Fortran
  result_soilt <- soilt(ALBEDO = ALBEDO, B = B, CUMDPT = CUMDPT, DOY = DOY,
                        DP = DP, HDAY = HDAY, NLAYR = NLAYR, SRAD = SRAD,
                        TAMP = TAMP, TAV = TAV, TAVG = TAVG, TMAX = TMAX,
                        TMA = TMA, WW = WW, PESW = PESW, DSMID = DSMID)
  
  # Retorna los resultados
  return(result_soilt)
}





