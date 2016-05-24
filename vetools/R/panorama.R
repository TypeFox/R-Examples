# Verified 1.3.18
# Version 4.0
panorama <-
function(collection, main, cut, 
                     ylab.push.factor = 10,
                     cut.col = "darkred",
                     cut.lty = 1,
                     cut.lwd = 2,
                     col = "RoyalBlue",
                     col.ramp = c("red", "pink","blue"),
                     col.line = "gray30",
                     mar = c(5, 4+ylab.push.factor, 3, 2),
                     cex.axis = 0.8,
                     cex.yaxis = 0.7,
                     xlab = "Year",
                     color.by.data = FALSE,
                     ...) {
        
        if ( length(names(collection$catalog)) != 0 ) {
                cat = collection$catalog; catalogo = list(); catalogo[[1]] = cat; rm("cat")
                dat = collection$data; datos = list(); datos[[1]] = dat; rm("dat")
        } else {
                catalogo = collection$catalog
                datos = collection$data
        }
        
        if ( length(catalogo) != length(datos) ) { stop("Collection: catalog and data lengths differ.") }
        
        disponibles <- function(x) { return( start(x)[1] : end(x)[1] ) }

        n = length(catalogo)

        xcol = col
        colf = function(x) { colorRamp(col.ramp)(x) }
        if ( ! missing(datos) ) { #&& ( color.by.data == FALSE ) ) {
                dpa = list()
                xcol = list()
                for ( k in 1 : n ) {
                        s = start(datos[[k]])
                        e = end(datos[[k]])
                        f = frequency(datos[[k]]) # = 12
                        kk = 1
                        an = array()
                        for ( a in s[1] : e[1] ) {
                                if  ( color.by.data == FALSE ) {
                                        an[[kk]] = sum(!is.na(window(datos[[k]], start=c(a,1), end=c(a, f), extend=T))) / f
                                } else {
                                        an[[kk]] = mean((window(datos[[k]], start=c(a,1), end=c(a, f), extend=T)), na.rm=TRUE)
                                }
                                dpa[[k]] = an
                                kk = kk + 1
                        }
                        xcol[[k]] = rgb(colf(an)/255)
                }
        } 
#         else if ( ( ! missing(datos) ) && ( color.by.data == TRUE ) ) {
#                 xcol = list()
#                 for ( k in 1 : n ) {
#                         xcol[[k]] = rgb(colf(datos[[k]])/255)
#                 }
#         }

        dis = unlist(lapply(datos, function(x) { c(start(x)[1], end(x)[1]) } ))
        ylabs = unlist(lapply(catalogo, function(x) { x$Name } ))
        xdat = range(dis)
        xlim = xdat + c(0, 1)
        xlim.names = c(0, 0)
        ylim.names = c(0, 4)
        ylim = c(0.5, n)
        old.par <- par(no.readonly = TRUE)
        on.exit(par(old.par))
        layout(1)
        par(bty="n", mar = mar, ...)
        plot(axes=F, xdat, ylim+ylim.names,type="n", xlab=xlab,ylab=NA, xlim=xlim+xlim.names, ylim=ylim+ylim.names)
        abline(h=seq(1,n,2), col=col.line, lty=3, lwd=1)
        if (!missing(cut)) { abline(v=cut, col=cut.col, lty=cut.lty, lwd=cut.lwd) }
        points(disponibles(datos[[1]]), rep(1, length(disponibles(datos[[1]]))), pch=22, bg=xcol[[1]])

        text(xdat[1]+5, ylim[2]-1.5+ylim.names[2],labels="Available data", pos=3, cex=0.85)
        points(xdat[1], ylim[2]-2+ylim.names[2], pch=22, bg=rgb(colf(1)/255))
        text(xdat[1], ylim[2]-2+ylim.names[2],labels="100%", pos=4, cex=0.85)
        points(xdat[1]+5, ylim[2]-2+ylim.names[2], pch=22, bg=rgb(colf(0.5)/255))
        text(xdat[1]+5, ylim[2]-2+ylim.names[2],labels="50%", pos=4, cex=0.85)
        points(xdat[1]+10, ylim[2]-2+ylim.names[2], pch=22, bg=rgb(colf(0.0)/255))
        text(xdat[1]+10, ylim[2]-2+ylim.names[2],labels="0%", pos=4, cex=0.85)

        if ( n > 1 ) {
                for ( f in 2:n ) {
                        if (missing(datos)) {
                                points(disponibles(datos[[f]]), rep(f, length(disponibles(datos[[f]]))), pch=22, bg=xcol, type="p")
                        } else {
                                points(disponibles(datos[[f]]), rep(f, length(disponibles(datos[[f]]))), pch=22, bg=xcol[[f]], type="p")
                        }
                }
        }
        axis(1)
        axis(2, 1:n, ylabs, hadj=1, las=1,tick=F, cex.axis=cex.yaxis)
        axis(4, 1:n, 1:n, las=1,tick=F, hadj=0.5, col.axis=col.line, cex.axis=cex.yaxis)
        if (missing(main)) {
                main="Longevity of stations"
        }
        title(main=main)
        invisible()
}
