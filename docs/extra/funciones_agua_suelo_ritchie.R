
ClasesTexturales <- function(arena, limo, arcilla, carbono, sistema_clasificacion = "USDA.TT") {

    # sistema_clasificacion: Single text string. Text code of the texture classification
    # system to be used for the classification of 'tri.data'.
    # Possible values are "none" (no classification plotted),
    # "USDA.TT" (USDA texture triangle),
    # "HYPRES.TT" (texture triangle of the European Soil Map),
    # "FR.AISNE.TT" (French texture triangle of the Aisne region soil survey),
    # "FR.GEPPA.TT" (French GEPPA texture triangle),
    # "DE.BK94.TT" (German texture triangle),
    # "UK.SSEW.TT" (Soil Survey of England and Wales),
    # "AU.TT" (Australian texture triangle),
    # "BE.TT" (Belgium texture triangle),
    # "CA.EN.TT" (Canadian texture triangle, with English class abbreviations)
    # "CA.FR.TT" (Canadian texture triangle, with French class abbreviations)

    # Creacion del data frame para la determinación de la textura correspondiente
    granulometria <- data.frame(
        "CLAY"  = arcilla,
        "SILT"  = limo,
        "SAND"  = arena,
        "OC"    = carbono
    )

    # Determinacion de la clase textural
    clase <- soiltexture::TT.points.in.classes(
        tri.data    = granulometria,
        class.sys   = sistema_clasificacion,
        PiC.type    = "t",
        tri.sum.tst = FALSE
    )

    # Devolver resultados
    return(clase)

}


FAOTextureClasses <- function(datos_suelo, profundidad_superficial = 30) {

    # Definicion de función para la clasificación de textura basada en la capa
    # superficial del suelo (0-30 cm) a partir del Mapa Global de Suelos
    # (FAO/Unesco, 1970-1980).
    # https://www.fao.org/3/aq361e/aq361e.pdf

    perfil_i <- datos_suelo %>%
        dplyr::select(perfil, capa, limite_superior, limite_inferior, arena, limo, arcilla) %>%
        dplyr::mutate(profunidad = limite_inferior - limite_superior)

    granulometria_media_superficial <- purrr::pmap_dfr(
        .l = perfil_i,
        .f = function(...) {
            row = tibble::tibble(...)
            cm_seq = seq(row$limite_inferior - row$profunidad, row$limite_inferior, 1)
            result = tibble::tibble(horizon_id = row$capa,
                                    horizon_depth = cm_seq) %>%
                dplyr::mutate(arena = row$arena,
                              limo = row$limo,
                              arcilla = row$arcilla)
        }
    ) %>%
        dplyr::filter(horizon_depth != 0) %>%
        dplyr::filter(horizon_depth <= config$params$profundidad$granulometria) %>%
        dplyr::group_by(horizon_id) %>%
        dplyr::summarise(arena = mean(arena),
                         limo = mean(limo),
                         arcilla = mean(arcilla),
                         cantidad_observaciones = n()) %>%
        dplyr::ungroup() %>%
        dplyr::mutate(peso_horizonte = cantidad_observaciones/sum(cantidad_observaciones)) %>%
        dplyr::summarise(arena = weighted.mean(x = arena, w = peso_horizonte),
                         limo = weighted.mean(x = limo, w = peso_horizonte),
                         arcilla = weighted.mean(x = arcilla, w = peso_horizonte))

    # Extraer promedio de capa superficial
    arena <- granulometria_media_superficial$arena
    limo <- granulometria_media_superficial$limo
    arcilla <- granulometria_media_superficial$arcilla


    # Clasificación textural de la FAO
    if (arcilla < 18 & arena >= 65) {
        clase <- 'gruesa'
    } else if ((arcilla < 35 & arena < 65) | (arcilla >= 18 & arena < 82)) {
        clase <- 'media'
    } else if (arcilla > 35) {
        clase <- 'fina'
    }

    # Devolver resultados
    return(clase)

}

#FAOTextureClasses(datos_suelo = datos_analiticos_perfiles,
#                               perfil_id = 'USDA 11')

FAOWaterHoldingCapacityClasses <- function(datos_suelo, profundidad_raices = 100, clases_completas = FALSE) {

    # Definicion de función para la clasificación de suelos en función
    # de su capacidad de almacenar agua. Las clases están basadas en
    # el Mapa GLobal de Suelos (FAO/Unesco, 1970-1980).
    # https://www.fao.org/3/aq361e/aq361e.pdf


    perfil_i <- datos_suelo %>%
        dplyr::select(perfil, capa, limite_superior, limite_inferior, punto_marchitez, capacidad_campo, saturacion) %>%
        dplyr::mutate(agua_disponible = capacidad_campo - punto_marchitez,
                      profunidad = limite_inferior - limite_superior)

    capacidad_alamacenaje <- purrr::pmap_dfr(
        .l = perfil_i,
        .f = function(...) {
            row = tibble::tibble(...)
            cm_seq = seq(row$limite_inferior - row$profunidad, row$limite_inferior, 1)
            result = tibble::tibble(horizon_id = row$capa,
                                    horizon_depth = cm_seq) %>%
                dplyr::mutate(punto_marchitez = row$punto_marchitez,
                              capacidad_campo = row$capacidad_campo,
                              saturacion = row$saturacion,
                              agua_disponible   = row$agua_disponible)
        }
    ) %>%
        # Eliminar la primer fila que produce profundidades de 0
        dplyr::filter(horizon_depth != 0) %>%
        dplyr::filter(horizon_depth <= profundidad_raices) %>%
        dplyr::pull(agua_disponible) %>%
        sum()*10

    if (clases_completas) {
        if (capacidad_alamacenaje >= 150) {
            clase <- '1'
        } else if (capacidad_alamacenaje >= 125 & capacidad_alamacenaje < 150) {
            clase <- '2'
        } else if (capacidad_alamacenaje >= 100 & capacidad_alamacenaje < 125) {
            clase <- '3'
        } else if (capacidad_alamacenaje >= 75 & capacidad_alamacenaje < 100) {
            clase <- '4'
        } else if (capacidad_alamacenaje >= 50 & capacidad_alamacenaje < 75) {
            clase <- '5'
        } else if (capacidad_alamacenaje >= 15 & capacidad_alamacenaje < 50) {
            clase <- '6'
        } else { # Almacenaje menor a 15 mm
            clase <- '7'
        }
    } else{
        if (capacidad_alamacenaje >= 150) {
            clase <- 'profundo'
        } else if (capacidad_alamacenaje >= 125 & capacidad_alamacenaje < 150) {
            clase <- 'media'
        } else {
            clase <- 'somero'
        }
    }

    # Devolver resultados
    return(clase)

}

#FAOWaterHoldingCapacityClasses(datos_suelo = datos_analiticos_perfiles,
#                               perfil_id = 'USDA 11')

FAOFertilityClasses <- function(datos_suelo, profundidad_superficial = 30) {

    # Definicion de función para la clasificación de suelos en función
    # de su capacidad de almacenar agua. Las clases están basadas en
    # el Mapa GLobal de Suelos (FAO/Unesco, 1970-1980).
    # https://www.fao.org/3/aq361e/aq361e.pdf


    perfil_i <- datos_suelo %>%
        dplyr::select(perfil, capa, limite_superior, limite_inferior,carbono) %>%
        dplyr::mutate(profunidad = limite_inferior - limite_superior)

    carbono_superficial <- purrr::pmap_dfr(
        .l = perfil_i,
        .f = function(...) {
            row = tibble::tibble(...)
            cm_seq = seq(row$limite_inferior - row$profunidad, row$limite_inferior, 1)
            result = tibble::tibble(horizon_id = row$capa,
                                    horizon_depth = cm_seq) %>%
                dplyr::mutate(carbono = row$carbono)
        }
    ) %>%
        # Eliminar la primer fila que produce profundidades de 0
        dplyr::filter(horizon_depth != 0) %>%
        dplyr::filter(horizon_depth <= config$params$profundidad$carbono) %>%
        dplyr::group_by(horizon_id) %>%
        dplyr::summarise(carbono = mean(carbono),
                         cantidad_observaciones = n()) %>%
        dplyr::ungroup() %>%
        dplyr::mutate(peso_horizonte = cantidad_observaciones/sum(cantidad_observaciones)) %>%
        dplyr::summarise(carbono = weighted.mean(x = carbono, w = peso_horizonte))

    if (carbono_superficial >= 1.2) {
        clase <- 'alta'
    } else if (carbono_superficial >= 0.7 & carbono_superficial <1.2) {
        clase <- 'media'
    } else { # carbono_superficial < 0.7
        clase <- 'baja'
    }
    # Devolver resultados
    return(clase)

}

#FAOFertilityClasses(datos_suelo = datos_analiticos_perfiles,
#                    perfil_id = 'USDA 11')


# Funciones
SoilRootGrowthFactor <- function(depth, clay) {

    # Factor de profundidad
    if (0 < depth & depth <= 30) {
        soil_depth_factor = 1
    } else if (30 < depth & depth <= 45) {
        soil_depth_factor = 0.8
    } else if (45 < depth & depth <= 60) {
        soil_depth_factor = 0.7
    } else if (60 < depth & depth <= 90) {
        soil_depth_factor = 0.5
    } else if (90 < depth & depth <= 120) {
        soil_depth_factor = 0.4
    } else if (120 < depth & depth <= 150) {
        soil_depth_factor = 0.3
    } else if (150 < depth & depth <= 180) {
        soil_depth_factor = 0.2
    } else if (180 < depth & depth <= 210) {
        soil_depth_factor = 0.1
    } else if (210 < depth) {
        soil_depth_factor = 0.05
    }

    # Factor de impedancia
    if (depth >= 30) {
        if (clay <= 32) {
            soil_impedance_factor = 1
        } else if (32 < clay & clay <= 40) {
            soil_impedance_factor = 0.4
        } else if (clay > 40) {
            soil_impedance_factor = 0.4
        }
    } else {
        soil_impedance_factor = 1
    }

    # Calcular factor de crecimiento
    soil_root_growth_factor = soil_depth_factor * soil_impedance_factor

    # Devolver resultados
    return(soil_root_growth_factor)
}


# Definicion de funcion para el cálculo del límite inferior
# de agua en el suelo según Ritchie et al. (1987,1989)

LowerLimitRitchie <- function(clay = NA, silt = NA, sand = NA) {

    # Check that clay and silt contents are present
    if (is.na(clay) | is.na(silt)) {
        stop('Clay and silt contents must be present in the dataset')
    }
    # Check that sum of all component yields 100%
    if (all(is.na(clay), is.na(silt), is.na(sand))) {
        total_sum = clay + silt + sand
        if (!all.equal(total_sum, 100)) {
            stop('The sum of all particle sizes must be equal to 100')
        }
    }
    # If sand content is not provided, estimate it as the difference
    # of clay and silt
    if (is.na(sand)) {
        sand = 100 - clay - silt
    }
    # Funcion para la estimacion del Punto de Marchitez Permanente
    # de la fraccion mineral basado en Ritchie and Crum (1989).
    # For sand higher than 75%
    if (sand > 75) {
        # Lower limit (LOL) for the mineral fraction
        LOLm = 18.8 - 0.168 * sand
        # For sand lower than 75%
    } else {
        # For sand lower than 75% and silt lower than 70%
        if (sand <= 75 & silt <= 70) {
            # Lower limit (LOL) for the mineral fraction
            LOLm = 3.62 + 0.444 * clay
        } else if (sand <= 75 & silt >= 70) {
            # Lower limit (LOL) for the mineral fraction
            LOLm = 5.0 + 0.0244 * clay ** 2
        }
    }

    return(LOLm)
}

#LOLm = LowerLimitRitchie(clay = clay, silt = silt)



# Approximating the influence of bulk density, organic matter and
# rock fragments on Extractable Soil Water Limits.

MineralBulkDensityRitchie <- function(clay = NA, silt = NA, sand = NA) {

    # Check that clay and silt contents are present
    if (is.na(clay) | is.na(silt)) {
        stop('Clay and silt contents must be present in the dataset')
    }
    # Check that sum of all component yields 100%
    if (all(is.na(clay), is.na(silt), is.na(sand))) {
        total_sum = clay + silt + sand
        if (!all.equal(total_sum, 100)) {
            stop('The sum of all particle sizes must be equal to 100')
        }
    }
    # If sand content is not provided, estimate it as the difference
    # of clay and silt
    if (is.na(sand)) {
        sand = 100 - clay - silt
    }

    # Estimation of Bulk Density for the mineral fraction (Dm)

    # For sand higher than 80%
    if (sand >= 80) {
        Dm = 1.709 - 0.01134 * clay
        # For sand between 20% and 80%
    } else if (20 <= sand & sand < 80) {
        Dm = 1.118 + 0.00816 * sand + clay * (0.00834 - 0.3606 / (100 - sand))
        # For sand lower than 20%
    } else if (sand < 20) {
        Dm = 1.453 - 0.004330 * sand
    }
    return(Dm)
}

#Dm = MineralBulkDensityRitchie(clay = clay, silt = silt)

FieldBulkDensityRitchie <- function(Dm = NA, sloc = NA) {

    # Estimation of Field Bulk Density by correcting the bulk density
    # from the mineral fraction by soil organic carbon
    Df = 100/(((100-sloc/0.58)/Dm) + (sloc/0.58/0.224))

    return(Df)
}

# Df = FieldBulkDensityRitchie(Dm = Dm, sloc = sloc)


DrainageUpperLimitRitchie <- function(clay = NA, silt = NA, sand = NA, sloc = NA, LOLm = NA, Dm = NA, Df = NA) {

    # Check that clay and silt contents are present
    if (is.na(clay) | is.na(silt)) {
        stop('Clay and silt contents must be present in the dataset')
    }
    # Check that sum of all component yields 100%
    if (all(is.na(clay), is.na(silt), is.na(sand))) {
        total_sum = clay + silt + sand
        if (!all.equal(total_sum, 100)) {
            stop('The sum of all particle sizes must be equal to 100')
        }
    }
    # If sand content is not provided, estimate it as the difference
    # of clay and silt
    if (is.na(sand)) {
        sand = 100 - clay - silt
    }

    if (sand >= 80) {
        DULc = LOLm + (0.01 * (42.3-0.381 * sand) * 100 - 3.5 * (Dm - Df) + 0.55*sloc/0.58)
    } else {
        DULc = LOLm + (0.01 * (10.79 + 0.05004 * silt) * 100 - 3.5 * (Dm - Df) + 0.55*sloc/0.58)
    }

    return(DULc)
}


# DrainageUpperLimitRitchie(clay = clay, silt = silt, sloc = sloc, LOLm = LOLm, Dm = Dm, Df = Df)

ksat <- function(sand, clay, soc, field_capacity, wilting_point,
                 saturation, bulk_density,
                 DF = 1, gravel = 0) {
    fcdf <- field_capacity
    wp <- wilting_point
    thetasdf = saturation
    mdens = bulk_density
    lambda <- (log(fcdf) - log(wp)) / (log(1500) - log(33))  # = 1/Beta
    theta_sdf_fcdf <- thetasdf - fcdf
    theta_sdf_fcdf <- ifelse(theta_sdf_fcdf < 0, 0, theta_sdf_fcdf) # FC ! > por.
    kbks <- (1 - gravel) / (1 - gravel * (1 - 1.5 * (mdens / 2.65)))
    ks <- 1930 * (theta_sdf_fcdf)^(3 - lambda) * kbks
    return(ks)
}
