# file OrdinalLogisticBiplot/R/BiplotDensity.R
# copyright (C) 2012-2013 J.C. Hernandez and J.L. Vicente-Villardon
#
#  This program is free software; you can redistribute it and/or modify
#  it under the terms of the GNU General Public License as published by
#  the Free Software Foundation; either version 2 or 3 of the License
#  (at your option).
#
#  This program is distributed in the hope that it will be useful,
#  but WITHOUT ANY WARRANTY; without even the implied warranty of
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#  GNU General Public License for more details.
#
#  A copy of the GNU General Public License is available at
#  http://www.r-project.org/Licenses/
#


BiplotDensity <- function(X, y = NULL, nlevels = max(y), grouplabels = 1:nlevels, ncontours = 6,
                           groupcols = 1:nlevels, img = TRUE, separate = FALSE, ncolors=20, ColorType=4,
                           xliml=-1,xlimu=1,yliml=-1,ylimu=1,plotInd = FALSE) {
        if (is.null(y))
                y = matrix(1, dim(X)[1], 1)
        nlevels = max(y)

        switch(ColorType, "1" = {colores = rainbow(ncolors)},
                          "2" = {colores = heat.colors(ncolors)},
                          "3" = {colores = terrain.colors(ncolors)},
                          "4" = {colores = topo.colors(ncolors)},
                          "5" = {colores = cm.colors(ncolors)})
        #print(colores)
        
        if(img){
            if(separate){
                ncols = round(sqrt(nlevels))
                nrows = ceiling(nlevels/ncols)
                par(mfrow = c(nrows, ncols))
                par(mar = c(0, 0, 0, 0), oma = c(4, 4, 0.5, 0.5))
                par(tcl = -0.25)
                par(mgp = c(2, 0.6, 0))
            }else{
                f1 <- kde2d(X[,1], X[,2], n = 400, lims = c(xliml, xlimu,yliml, ylimu))
                image(f1, zlim = c(0, 1), col=colores, add=TRUE)
            }
        }

        for (i in 1:nlevels) {
                x = X[which(y == i), ]
                if(!is.null(nrow(x))){
                    f1 <- kde2d(x[, 1], x[, 2], n = 400, lims = c(xliml, xlimu,yliml, ylimu))
                    contours = round(seq(0, max(f1$z), length.out = ncontours), digits = 2)
                    #print(contours)
                    centre = apply(x, 2, mean)
                }else{
                    centre = x
                }
                if (separate){                        
                    if(img){
                        plot(X[, 1], X[, 2], pch = y, cex = 0, asp = 1, axes = FALSE,xlim=c(xliml,xlimu),ylim=c(yliml,ylimu))
                        if(!is.null(nrow(x))){
                          image(f1, zlim = c(0, 1), col = colores, add = TRUE)
                        }
                        if(plotInd){
                          points(x[, 1], x[, 2], pch = y, cex = 0.1, asp = 1, col=groupcols[i])
                        }
                        box(col = "grey60")
                    }else{
                        if(plotInd){
                            plot(X[, 1], X[, 2], pch = y, cex = 0.1, asp = 1, axes = FALSE,xlim=c(xliml,xlimu),ylim=c(yliml,ylimu))
                        }
                    }
                }
                if(!is.null(nrow(x))){
                    contour(f1, levels = contours, col = groupcols[i], add = TRUE)
                }
                box(col = "grey60")
                
                points(centre[1], centre[2])
                text(centre[1], centre[2], labels = grouplabels[i],pos=2,offset=0.1)
        }
}

