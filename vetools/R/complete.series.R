# Verified 1.3.18
# ASUMPTIONS:
# input: collection = c(data, catalog)
# k.ubic = k.m [KMEANS] Modelo  de Conglomedaro Jerarquico kmeans()
# modelo = fit [LISTA] Modelos de efectos mensuales LISTA, se extrae modelos$fitted.values
# NOTE:
# It is assumed that "modelo", does not have any intermediary NAs, like using the command: 
#    fit[[k]] = lm(datos.m.sqrt[[k]] ~ ZZ - 1, singular.ok=T, na.action=na.omit)
complete.series <-
function( collection, model, k.ubic=NA, centers=3, nstart=3, weps=0.05, MAX.ITER=100, AEM.debug=T ) {
        indices = function(x, c) {
                if ( is.logical(c) ) {
                        return( (1 : length(x))[c] )
                } else {
                        return( (1 : length(x))[x == c] )
                }
        }
        
        datos = collection$data
        modelo = model
        
        if ( is.na(k.ubic) ) {
                cat("Building cluster from model...\n")
                Mat = lapply(modelo, function(x) { x$coefficients })
                Mat = t(matrix(unlist(Mat), ncol=length(modelo)))
                cat("k.ubic = kmean(modelo, centers =", centers, "nstart = ", nstart, ")\n")
                k.ubic = kmeans(x=Mat, centers=centers, nstart=nstart)
                cat("Groups found:", k.ubic$cluster, "\n")
        }

        # Sanity checks
        cl = length(unique(k.ubic$cluster))
        for ( i in 1 : cl ) {
                if (sum( k.ubic$cluster == i ) == 1) { stop("Each group must have at least two individuals... Decrease \"centers\" a/o \"nstart\".")}
        }

        # Convert data to a list with length of the longest vector.
        # Before migrating to list, move all NAs withing vector to the end. Does this not scamble seasonality (???)
        M = unlist(lapply(datos, function(x) { length(x) }))
        m = max(M)
        y = lapply(modelo, function(x) { as.numeric(x$residuals) })
        modelo.res = lapply(y, function(x, m) { c(x, rep(NA, m - length(x))) }, m=m)
        modelo.res = matrix(unlist(modelo.res), ncol=length(modelo.res))

        Y.compl = list()
        prop.falt = c()

        for ( kkk in 1 : cl ) {
                if ( AEM.debug ) { cat("Cluster group ", kkk, "of", cl, "\n") }
                res = modelo.res[ , k.ubic$cluster == kkk]
                cat("Cluster group of stations ", (1:length(datos))[k.ubic$cluster == kkk], "\n")
                # START
                n = dim(res)[1]
                k = dim(res)[2]
                prop.falt[kkk] = sum(is.na(res)) / (n * k)
                res.aux = res  # Initial guess
                sd.j = apply(res, 2, sd, na.rm=T)
                for (i in 1 : dim(res)[1]) {
                        for (j in 1 : dim(res)[2]) {
                                if ( is.na(res[i, j]) ) { res.aux[i, j] = rnorm(1, mean=0, sd=sd.j[j]) }
                        }
                }
                mu.ini = apply(res.aux, 2, mean)
                V.ini = cov(res.aux)
                res.compl = res
                cambio = 1
                ult = cambio
                iter = 0
                alto = FALSE
                while ( alto == FALSE ) { # AEM start here with table data 'kkk' geo group.
                        iter = iter + 1
                        if (iter > MAX.ITER) { stop("MAX.ITER reached. Increase value and turn on debugging with AEM.debug=T.") }
                        com.mu = mu.ini
                        com.V = V.ini
                        # Expectation
                        x.suma = rep(0, k)
                        V.suma = matrix(0, k, k)
                        V.aux = matrix(NA, k, k)
                        for (i in 1 : n) {
                                x.i = res[i, ]
                                falta = is.na(x.i)
                                cuantos = sum(falta)
                                if ( (cuantos > 0) & (cuantos < k) ) {
                                        falta.ord = sort(seq(k) * falta)
                                        falta = falta.ord[seq(k - cuantos + 1, k)]
                                        V.12 = V.ini[falta, -falta];
                                        if ( cuantos == 1 ) {
                                                V.12 = matrix(c(V.12), nrow=1, ncol=k-1)
                                        }
                                        V.22 = as.matrix(V.ini[-falta, -falta]);
                                        V.22 = solve(V.22)
                                        x.i[falta] = mu.ini[falta] + ( V.12 %*% V.22 ) %*% ( x.i[-falta] - mu.ini[-falta] )
                                        x.suma = x.suma + x.i
                                        V.aux[falta, falta] = as.matrix( (V.ini[falta, falta] - V.12 %*% V.22 %*% t(V.12) ) + x.i[falta] %*% t(x.i[falta]) )
                                        V.aux[-falta, -falta] = as.matrix( x.i[-falta] %*% t(x.i[-falta]) )
                                        V.aux.21 = x.i[-falta] %*% t(x.i[falta])
                                        if ( k - cuantos == 1 ) {
                                                V.aux.21 = matrix(c(V.aux.21), nrow=1, ncol=cuantos)
                                        }
                                        V.aux[-falta, falta] = V.aux.21
                                        V.aux[falta, -falta] = t(V.aux.21)
                                        V.suma = V.suma + V.aux
                                } else if ( cuantos == k ) {
                                        x.i = mu.ini
                                        x.suma = x.suma + x.i
                                        V.suma = V.suma + ( V.ini - mu.ini %*% t(mu.ini) )
                                } else {
                                        x.suma = x.suma + x.i
                                        V.suma = V.suma + x.i %*% t(x.i)
                                }
                                if (cambio <= weps) {
                                        res.compl[i, ] = x.i;
                                        alto = TRUE
                                } # Vector is filled
                        }
                        # Maximize
                        mu.ini = (1/n) * x.suma
                        V.ini= (1/n) * V.suma - ( mu.ini %*% t(mu.ini) )
                        ultimo = cambio
                        cambio = max( max( com.mu - mu.ini ), max( V.ini - com.V ) )
                        if ( AEM.debug ) { cat("Cluster:", kkk, "iter:", iter, "change update:", cambio, "\n") }
                        try(is.na(cambio)) # fail-safe
                        if ( ( iter > 1 ) & ( ultimo < cambio ) ) { stop("AEM stoped converging.") } # Puede deberse a un desafortunado agrupamiento (?)
                } # AEM finishes here
                Y.compl[[kkk]] = res.compl
        }

        # Back to ts class. UGLY (!!)
        datos.completados = datos
        for ( cl.in in 1 : cl ) {
                A = Y.compl[[cl.in]]
                for ( l in indices(k.ubic$cluster, cl.in) ) {
                        idx.data = indices(datos[[l]], !is.na(datos[[l]]))
                        idx.na = indices(datos[[l]], is.na(datos[[l]]))
                        # Con NAs
                        ret = datos[[l]]
                        if ( AEM.debug ) { plot(ret, typ="l", lwd=4, main=paste("Station", l, "(completed with AEM)"), xlab="", ylab="[mm]") }
                        # Monthly adjustment, NEED TO SINCHRONIZE BY MONTH:
                        # SERIES MAY NOT BEGIN AT JANUARY!
                        mo = modelo[[l]]$fitted.values[apply(matrix(idx.na, ncol=1), 1, function(x, ini) { a=((ini - 1 + x) %% 12);if (a==0) {return(12)} else {return(a)} }, ini=start(ret)[2])]
                        ret[idx.na] = mo + A[length(idx.data)+(1:length(idx.na)), 1]
                        if ( AEM.debug ) { lines(ret, col="blue", lwd=2) }
                        datos.completados[[l]] = ret
                        if ( AEM.debug ) { scan() }
                }
        }
        col = list(catalog = collection$catalog, data = datos.completados)
        class(col) <- "Catalog"
        return( col )
}
