# Verified 1.3.18
# Version 3.0
read.HIDROX <-
function(file, state=NA, altitudes=NA, serial=NA, unit=NA) {
        a <- read.csv(file, header=T)
        Nombres <- levels(a[,6])
        datos.m = list()
        catalogo = list()
        k = 0
        k.dis = 0
        m = length(Nombres)
        for ( n in Nombres ) {
                k = k + 1
                k.dis = k.dis + 1
                cat("Processing station: [",k.dis,"/",m,"] ", n, "\n")
                idx = factor(a$Estacion) == n        
                b = a[idx, ]
                mk = matrix(unique(cbind(b$Longitud, b$Latitud)), ncol=2)
                colnames(mk)<-c("Longitud", "Latitud")
                print(mk)
                if ( length(unique(b$Latitud)) > 1 ) {
                        cat("Note: Station appears to have more than one Lat/Long. Could there be more than one station with same name?\n")
                }
                bb = b
                k = k - 1
                for ( i in 1 : nrow(mk) ) {
                        k = k + 1
                        b = bb[  ( (bb$Longitud == mk[i, "Longitud"]) & (bb$Latitud == mk[i, "Latitud"])  ) ,  ]
                        b$Fecha <- as.Date(b$Fecha)
                        b.xts = xts::xts(b$Valor, b$Fecha)
                        cat("   + Processing station: [",k.dis,"/",m,"] ", n, ": sub estation ", i, "of", nrow(mk), "\n")
                        datos.m[[k]] = xts2ts(b.xts)
                        catalogo[[k]] = list(
                                Name = gsub(" *$", "", n),
                                Serial = ifelse(is.na(serial), NA, serial[k]),
                                Ss = nrow(mk),
                                S = i,
                                Altitude = ifelse(is.na(altitudes), NA, altitudes[k]),
                                Latitude = b$Latitud[1],
                                Longitude = b$Longitud[1],
                                Measure.code = as.character(a$Unidad.de.medida[1]),
                                Measure.unit = ifelse(is.na(unit), as.character(a$Unidad.de.medida[1]), unit),
                                Install = time2ym(start(b.xts)),
                                Start = time2ym(start(b.xts)),
                                State = ifelse(is.na(state), NA, state),
                                Avble.yrs = sort(unique(as.numeric(format(b$Fecha, "%Y"))))
                        )
                }
        }
        col = list(data = datos.m, catalog = catalogo)
        class(col) <- "Catalog"
        return(col)
}
