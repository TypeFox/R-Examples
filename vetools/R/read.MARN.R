# Verified 
# Version 3.0
read.MARN <-
function (file) {
        con <- file(file, open='r', encoding='cp1252')
        cat("Note: Missing values will be coerced to NAs. Ignore coercion warnings.\n")
        line = readLines(con, n=1)
        ignorar_estacion = FALSE
        ignorar_variable = FALSE
        ene = mar = may = jul = ago = oct = dic = array(dim=31)
        abr = jun = sep = nov = array(dim=30)
        feb = array(dim=28)
        SERIE = numeric(0)
        anios.disponibles = numeric(0)
        while (length(line)) {
                if ((!ignorar_variable) & (str_detect(line, 'DIVISION DE HIDROLOGIA,'))) {
                        line = readLines(con, n = 1)
                        Variable.tmp = unlist(strsplit(line, '\\.'))
                        Variable = str_trim(Variable.tmp[1])
                        Variable.Codigo = as.numeric(str_trim(Variable.tmp[2]))
                        cat("Variable :", Variable, "\n")
                        cat("Code     :", Variable.Codigo, "\n")
                        ignorar_variable = TRUE
                }
                # Non empty line:
                if (str_detect(line, 'ESTACION')) {
                        an.dis = unlist(strsplit(line, ':'))
                        anios.disponibles = c(anios.disponibles, as.numeric(an.dis[length(an.dis)]))
                        n = length(anios.disponibles)
                        if ( n > 1) {
                                for (y in (anios.disponibles[n-1]+1):anios.disponibles[n]){
                                        if (y == anios.disponibles[n]) { break }
                                        if ( lubridate::leap_year(y) ) {
                                                cat("Note: Padding year: ", y, " as leap year\n")
                                                SERIE = c(SERIE, rep(NA, 366))
                                        } else {
                                                cat("Note: Padding year: ", y, " as regular year\n")
                                                SERIE = c(SERIE, rep(NA, 365))
                                        }
                                }
                        }
                }
                if ((!ignorar_estacion) & (str_detect(line, 'ESTACION'))) {
                        tmp = unlist(strsplit(line, "SERIAL:"))
                        Estacion = unlist(strsplit(tmp[1], ':'))
                        Estacion = str_trim(Estacion[2])
                        tmp = unlist(strsplit(line, "SERIAL: "))
                        Serial = unlist(strsplit(tmp[2], ' '))
                        Serial = str_trim(Serial[1])
                        tmp = unlist(strsplit(line, "EDO: "))
                        Estado = unlist(strsplit(tmp[2], ' '))
                        Estado = str_trim(Estado[1])
                        tmp = unlist(strsplit(line, "EDO: .*A.*O: "))
                        Inicio = as.numeric(tmp[2])
                        cat("Station name: ", Estacion, " (",Serial, "). State: ", Estado, "\n", sep="")
                        line = readLines(con, n = 1)
                        if (!str_detect(line, 'LATITUD')) { stop("Parse error, LATITUD line expected but not found.") }
                        tmp = sub('^.*\\*[ ]*', '', line)
                        tmp = unlist(strsplit(tmp, '[ ]+'))
                        # LONGITUD
                        ch = "\u00a7"
                        tmp.proc = unlist(strsplit(tmp[5], ch))
                        Long = as.numeric(tmp.proc[1])
                        tmp.proc = unlist(strsplit(tmp.proc[2], "'"))
                        Long = Long + as.numeric(tmp.proc[1]) / 60
                        tmp.proc = unlist(strsplit(tmp.proc[2], "\""))
                        Long = Long + as.numeric(tmp.proc[1]) / 3600
                        # LATITUD
                        # tmp.proc = unlist(strsplit(tmp[3], "\xc2\xa7"))
                        tmp.proc = unlist(strsplit(tmp[3], ch))
                        Lat = as.numeric(tmp.proc[1])
                        tmp.proc = unlist(strsplit(tmp.proc[2], "'"))
                        Lat = Lat + as.numeric(tmp.proc[1]) / 60
                        tmp.proc = unlist(strsplit(tmp.proc[2], "\""))
                        Lat = Lat + as.numeric(tmp.proc[1]) / 3600
                        # ALTITUD
                        Altitud = as.numeric(tmp[7])
                        cat("Long: ", Long, ", Lat: ", Lat, ", Altitude (MSNM): ", Altitud, "\n", sep="")
                        Instalacion = as.Date(paste("01",tmp[10],sep="/"), "%d/%m/%Y")
                        cat("Instalation date: ", format(Instalacion), ", Initial year of operation: ", Inicio, "\n", sep="")
                        ignorar_estacion = TRUE
                        next
                }
                if (str_detect(line, 'DIA   ENE   FEB')) {
                        line = readLines(con, n = 1)
                        for ( d in 1:28 ) {
                                line = readLines(con, n = 1)
                                if (str_detect(line, '^ $')) { 
                                        line = readLines(con, n = 1);
                                }
                                tmp = unlist(strsplit(line, '[ ]+'))
                                tmp[tmp=="*"] = -9999
                                tmp = as.numeric(tmp[3:14])
                                ene[d] = tmp[1]
                                feb[d] = tmp[2]
                                mar[d] = tmp[3]
                                abr[d] = tmp[4]
                                may[d] = tmp[5]
                                jun[d] = tmp[6]
                                jul[d] = tmp[7]
                                ago[d] = tmp[8]
                                sep[d] = tmp[9]
                                oct[d] = tmp[10]
                                nov[d] = tmp[11]
                                dic[d] = tmp[12]
                        }
                        line = readLines(con, n = 1)
                        tmp = unlist(strsplit(line, '[ ]+')) # day 29
                        d = 29
                        if (length(tmp) == 12+3) { # leap year
                                y.leap = TRUE
                                tmp[tmp=="*"] = -9999
                                tmp = as.numeric(tmp[3:14])
                                ene[d] = tmp[1]
                                feb.bis = tmp[2]
                                mar[d] = tmp[3]
                                abr[d] = tmp[4]
                                may[d] = tmp[5]
                                jun[d] = tmp[6]
                                jul[d] = tmp[7]
                                ago[d] = tmp[8]
                                sep[d] = tmp[9]
                                oct[d] = tmp[10]
                                nov[d] = tmp[11]
                                dic[d] = tmp[12]
                        } else { # regular year
                                y.leap = FALSE
                                tmp[tmp=="*"] = -9999
                                tmp = as.numeric(tmp[3:13])
                                ene[d] = tmp[1]
                                # NO HAY
                                feb.bis = -5555
                                mar[d] = tmp[2]
                                abr[d] = tmp[3]
                                may[d] = tmp[4]
                                jun[d] = tmp[5]
                                jul[d] = tmp[6]
                                ago[d] = tmp[7]
                                sep[d] = tmp[8]
                                oct[d] = tmp[9]
                                nov[d] = tmp[10]
                                dic[d] = tmp[11]
                        }
                        d = 30
                        line = readLines(con, n = 1)
                        tmp = unlist(strsplit(line, '[ ]+')) # DIA 30
                        tmp[tmp=="*"] = -9999
                        tmp = as.numeric(tmp[3:13])
                        ene[d] = tmp[1]
                        # Skip feb
                        mar[d] = tmp[2]
                        abr[d] = tmp[3]
                        may[d] = tmp[4]
                        jun[d] = tmp[5]
                        jul[d] = tmp[6]
                        ago[d] = tmp[7]
                        sep[d] = tmp[8]
                        oct[d] = tmp[9]
                        nov[d] = tmp[10]
                        dic[d] = tmp[11]
                        d = 31
                        line = readLines(con, n = 1)
                        tmp = unlist(strsplit(line, '[ ]+')) # day 31
                        tmp[tmp=="*"] = -9999
                        tmp = as.numeric(tmp[3:9])
                        ene[d] = tmp[1]
                        # Skip feb
                        mar[d] = tmp[2]
                        # Skip
                        may[d] = tmp[3]
                        # Skip
                        jul[d] = tmp[4]
                        ago[d] = tmp[5]
                        # Skip
                        oct[d] = tmp[6]
                        # Skip
                        dic[d] = tmp[7]
                        # SERIES
                        if ( y.leap == FALSE ) {
                                SERIE = c(SERIE, ene, feb, mar, abr, may, jun, jul, ago, sep, oct, nov, dic)        
                        } else {
                                SERIE = c(SERIE, ene, feb, feb.bis, mar, abr, may, jun, jul, ago, sep, oct, nov, dic)
                        }
                }
                
                # Ignore any other line: MED, etc...
                line = readLines(con, n = 1)
        }
        close(con)
        cat("Years with available measurements: ", anios.disponibles, "\n")
        desc = list(Name = Estacion,
                    Serial = Serial,
                    State = Estado,
                    Install = Instalacion,
                    Start = Inicio,
                    Avble.yrs = anios.disponibles,
                    Longitude = Long,
                    Latitude = Lat,
                    Altitude = Altitud,
                    Measure.unit = Variable,
                    Measure.code = Variable.Codigo)
        
        pr = ts(SERIE, frequency=365.25, start=Inicio)
        col = list(data = pr, catalog = desc)
        class(col) <- "Catalog"
        return(col)
}
