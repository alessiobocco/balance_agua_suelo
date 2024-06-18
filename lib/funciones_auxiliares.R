# Función para crear una serie híbrida de tiempo
# 
# Esta función genera una serie híbrida de tiempo combinando datos climáticos
# de dos periodos: datos actuales y datos históricos. La serie se genera a partir
# de una fecha de inicio, una fecha de inflexión y una fecha de fin, utilizando
# datos históricos a partir del año de inflexión.
#
# Argumentos:
# id_estacion: Identificador de la estación climática.
# datos: Información climatica.
# fecha_inicio: Fecha de inicio del periodo actual (formato "YYYY-MM-DD").
# fecha_inflexion: Fecha de inflexión, punto de cambio entre datos actuales e históricos (formato "YYYY-MM-DD").
# fecha_fin: Fecha de fin del periodo actual (formato "YYYY-MM-DD").
# year_inflexion: Año de inicio para los datos históricos.
#
# Retorno:
# Una lista que contiene:
# - fechas: Secuencia completa de fechas desde el inicio hasta el final del periodo híbrido.
# - datos: Data frame combinado de datos actuales e históricos ordenados cronológicamente.

pr_crear_serie <- function(id_estacion, datos, fecha_inicio, fecha_inflexion, fecha_fin, year_inflexion) {
  # Convertir fechas a objetos Date
  fecha_inicio <- as.Date(fecha_inicio)
  fecha_inflexion <- as.Date(fecha_inflexion)
  fecha_fin <- as.Date(fecha_fin)
  
  # Extraer el mes y día de la fecha de inflexión
  monthday <- format(fecha_inflexion, "%m-%d")
  
  # Crear fecha híbrida de inicio utilizando el año de inflexión
  fecha_inicio_hibrida <- as.Date(paste(year_inflexion, format(fecha_inicio, "%m-%d"), sep="-"))
  
  # Ajustar la fecha de inicio de inflexión si es el 29 de febrero
  if(monthday == "02-29") {
    fecha_inicio_inflexion <- as.Date(paste(year_inflexion, "03-01", sep="-"))
  } else {
    fecha_inicio_inflexion <- as.Date(paste(year_inflexion, monthday, sep="-")) + 1
  }
  
  # Calcular la fecha de fin de inflexión
  fecha_fin_inflexion <- fecha_inicio_inflexion + as.numeric(difftime(fecha_fin, fecha_inflexion, units="days")) - 1
  
  # Filtrar los datos climáticos actuales
  datos_actuales <- subset(datos, station_id == id_estacion & date >= fecha_inicio & date <= fecha_inflexion)
  
  # Filtrar los datos climáticos históricos
  datos_historicos <- subset(datos, station_id == id_estacion & date >= fecha_inicio_inflexion & date <= fecha_fin_inflexion)
  
  # Añadir columna de orden para combinar y ordenar los datos
  datos_actuales$orden <- 1
  datos_historicos$orden <- 2
  datos_combinados <- rbind(datos_actuales, datos_historicos)
  
  # Ordenar los datos combinados por orden y fecha
  datos_combinados <- datos_combinados[order(datos_combinados$orden, datos_combinados$date),]
  
  # Generar secuencia completa de fechas desde la fecha híbrida de inicio hasta la fecha de fin de inflexión
  fechas_generadas <- seq(from = fecha_inicio_hibrida, to = fecha_fin_inflexion, by = "day")
  
  return(list("fechas" = fechas_generadas, "datos" = datos_combinados))
}
