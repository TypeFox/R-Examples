# Verified 1.3.18
# Version 3.0
read.MINAMB <-
function (file, state = NA, YSPLIT = 20) {
        con <- file(file, open='r', encoding='cp1252')
        cat("Note: Missing values will be coerced to NAs. Ignore coercion warnings.\n")
        cat("Note: Two digit date split is:", YSPLIT, "\n")
        meses.abb.pos <- function(x) {
                m = c("ENE", "FEB", "MAR", "ABR", "MAY", "JUN", "JUL", "AGO", "SEP", "OCT", "NOV", "DIC")
                return( ( 1 : 12 )[ m == x ] )
        }
        line = readLines(con, n=1)
        es.cnt = 0
        V.meta = list()
        V = list()
        while (length(line)) {
                if (str_detect(line, '^$')) { line = readLines(con, n = 1); next }
                if (str_detect(line, 'ESTACION')) {
                        # Solo interesa nombre de stacion:
                        Estacion = str_extract(line, ':.*$')
                        Estacion = str_trim(substring(Estacion, 2))
                        line = readLines(con, n = 1)
                        next
        
                }
                else if (str_detect(line, 'TIPO  SERIAL   DATO')) {
                        line = readLines(con, n = 1)
                        next
                }
                else if (str_detect(line, 'BASE DE DATOS')) {
                        tmp = sub('^.*\\*[ ]*', '', line)
                        tmp = unlist(strsplit(tmp, '[ ]+'))
                        Medicion = tmp[1]
                        Estacion = Estacion
                        Serial = tmp[2]
                        MSNM = as.integer(tmp[4])
                        Lat = as.integer(str_sub(tmp[5],0,1))
                        Lat = Lat + as.integer(str_sub(tmp[5],2,3))/60
                        Lat = Lat + as.integer(str_sub(tmp[5],4,5))/3600
                        Long = as.integer(str_sub(tmp[6],0,2))
                        Long = Long + as.integer(str_sub(tmp[6],3,4))/60
                        Long = Long + as.integer(str_sub(tmp[6],5,6))/3600
                        Long = -Long
                        Unidades = tmp[7]
                        Instalacion = c(as.numeric(str_sub(tmp[8],5,6)), meses.abb.pos(str_sub(tmp[8],0,3)))
                        if ( Instalacion[1] > YSPLIT ) {
                                Instalacion[1] = 1900 + Instalacion[1]
                        } else {
                                Instalacion[1] = 2000 + Instalacion[1]
                        }
                        if ( ! is.na(tmp[9]) ) {
                                if ( str_detect(tmp[9], "FUNC.") ) {
                                        Final = c(NA, NA)
                                } else {
                                        Final = c(as.numeric(str_sub(tmp[9],5,6)), meses.abb.pos(str_sub(tmp[9],0,3)))
                                        if ( Final[1] > YSPLIT ) {
                                                Final[1] = 1900 + Final[1]
                                        } else {
                                                Final[1] = 2000 + Final[1]
                                        }
                                }
                        } else {
                                Final = c(NA, NA)
                        }
                        es.cnt = es.cnt + 1
                        V.meta[[es.cnt]] = list( Estacion = unlist(Estacion),
                                               MSNM = unlist(MSNM),
                                               Lat = unlist(Lat),
                                               Long = unlist(Long),
                                               Medicion = unlist(Medicion),
                                               Unidades = unlist(Unidades),
                                               Instalacion = unlist(Instalacion),
                                               Final = unlist(Final),
                                               Serial = unlist(Serial)
                                                 )
                        line = readLines(con, n = 1)
                }
                else if (str_detect(line, 'ENE      FEB      MAR')) {
                        # Leer datos
                        aquired = F
                        tmp = ''
                        tmp.y = vector()
                        tmp.1 = tmp.2 = tmp.3 = tmp.4  = tmp.5  = tmp.6  = vector()
                        tmp.7 = tmp.8 = tmp.9 = tmp.10 = tmp.11 = tmp.12 = vector()
                        tmp.t = vector()
                        while (TRUE) {
                                line = readLines(con, n = 1)
                                if (!length(line)) { break }
                                if ( (str_detect(line, '^$')) & (!aquired) ) { next }
                                if ( (str_detect(line, '^$')) & (aquired) ) { break }
                                aquired = T
                                tmp    = unlist(strsplit(line, '[ ]+'))
                                tmp.y  = c(tmp.y, as.numeric(tmp[2]))
                                tmp.1  = c(tmp.1, as.numeric(tmp[3]))
                                tmp.2  = c(tmp.2, as.numeric(tmp[4]))
                                tmp.3  = c(tmp.3, as.numeric(tmp[5]))
                                tmp.4  = c(tmp.4, as.numeric(tmp[6]))
                                tmp.5  = c(tmp.5, as.numeric(tmp[7]))
                                tmp.6  = c(tmp.6, as.numeric(tmp[8]))
                                tmp.7  = c(tmp.7, as.numeric(tmp[9]))
                                tmp.8  = c(tmp.8, as.numeric(tmp[10]))
                                tmp.9  = c(tmp.9, as.numeric(tmp[11]))
                                tmp.10 = c(tmp.10, as.numeric(tmp[12]))
                                tmp.11 = c(tmp.11, as.numeric(tmp[13]))
                                tmp.12 = c(tmp.12, as.numeric(tmp[14]))
                                tmp.t  = c(tmp.t, as.numeric(tmp[15]))
                                V[[es.cnt]] = data.frame( Year      = tmp.y,
                                                          January   = tmp.1,
                                                          February  = tmp.2,
                                                          March     = tmp.3,
                                                          April     = tmp.4,
                                                          May       = tmp.5,
                                                          June      = tmp.6,
                                                          July      = tmp.7,
                                                          August    = tmp.8,
                                                          September = tmp.9,
                                                          October   = tmp.10,
                                                          November  = tmp.11,
                                                          December  = tmp.12,
                                                          Annual     = tmp.t )
                        }
                }
                # Ignore any other line: MED, etc...
                line = readLines(con, n = 1)
        }
        close(con)
        
        for ( es.cnt in 1 : length(V) ) {
                if ( is.na(V.meta[[es.cnt]]$Final[1]) ) {
                        V.meta[[es.cnt]]$Final = c(V[[es.cnt]]$Year[length(V[[es.cnt]]$Year)], 12)
                }
        }
        V.meta = t(simplify2array(V.meta))
        desc = data.frame(V.meta)
        desc = list( Nombre = unlist(desc$Estacion),
              MSNM = unlist(desc$MSNM),
              Lat = unlist(desc$Lat),
              Long = unlist(desc$Long),
              Variable.Codigo = unlist(desc$Medicion),
              Variable = unlist(desc$Unidades),
              Disponible = list(Instalacion=desc$Instalacion, Final=desc$Final),
              Instalacion = desc$Instalacion,
              Serial = unlist(desc$Serial))
        
        j = sapply(V, function(x) { x$Year[1] })
        desc$Inicio = cbind(matrix(j, length(V), 1), rep(1, length(V)))
        VV = list()
        V.ts = list()
        for ( i in 1:length(V) ) {
                VV[[i]] = list(
                        Name = desc$Nombre[[i]],
                        Altitude = desc$MSNM[[i]],
                        Latitude = desc$Lat[[i]],
                        Longitude = desc$Long[[i]],
                        Measure.code = desc$Variable.Codigo[[i]],
                        Measure.unit = desc$Variable[[i]],
                        Install = desc$Instalacion[[i]],
                        Start = desc$Instalacion[[i]],
                        State = ifelse(is.na(state), NA, state),
                        Avble.yrs = seq(desc$Inicio[i, 1], desc$Disponible$Final[[i]][1]),
                        Serial = desc$Serial[[i]]
                        )
                V.ts[[i]] = matrix(t(as.matrix(V[[i]][, 2:13])), nrow=1)
                V.ts[[i]] = ts(as.numeric(V.ts[[i]]), start=V[[i]]$Year[1], frequency=12)
        }
        col = list(data = V.ts, catalog = VV)
        class(col) <- "Catalog"
        return(col)
        
}
