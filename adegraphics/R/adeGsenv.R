## At the loading of the package, creation of an environment .ADEgEnv to store:
## - the list of graphical parameters   
## - the theme adeg                     
## - the last plotted graphics          

.ADEgEnv <- new.env()

.onLoad <- function(libname, pkgname) {
    assign("padegraphic",
           list(p1d = list(horizontal = TRUE, reverse = FALSE, rug = list(draw = TRUE, tck = 0.5, margin = 0.07, line = TRUE)),
                
                parrows = list(angle = 15, ends = "last", length = 0.1),
                
                paxes = list(aspectratio = "iso", draw = FALSE, x = list(draw = TRUE), y = list(draw = TRUE)),
                
                pbackground = list(col = "white", box = TRUE),
                
                pellipses = list(alpha = 0.5, axes = list(draw = TRUE, col = "black", lty = 4, lwd = 1), border = "black", col = "transparent", lty = 1, lwd = 1),
                
                pgrid = list(col = "grey", draw = TRUE, lty = 1, lwd = 1, nint = 5, text = list(cex = 1, col = "black", pos = "topright")),
                
                plabels = list(alpha = 1, cex = 1, col = "black", srt = "horizontal", optim = FALSE, 
                    boxes = list(alpha = 1, border = "black", col = "white", draw = TRUE, lwd = 1, lty = 1)),
                
                plegend = list(drawKey = TRUE, drawColorKey = FALSE, size = 1), 
                
                plines = list(col = "black", lty = 1, lwd = 1),
                
                pnb = list(edge = list(col = "black", lwd = 1, lty = 1), node = list(pch = 20, cex = 1, col = "black", alpha = 1)),
                
                porigin = list(alpha = 1, col = "black", draw = TRUE, include = TRUE, lty = 1, lwd = 1, origin = c(0, 0)),
                
                ppalette = list(quanti = colorRampPalette(c("white", "black")),
                    quali = function(n, name = "Set1") {
                        if(n > 9)
                            return(rainbow(n))
                        else if(n > 2)
                            return(brewer.pal(n, name))
                        else
                            return(brewer.pal(n + 2, name)[1:n])
                    }),  ## see http://colorbrewer2.org/
                
                ppoints = list(alpha = 1, cex = 1, col = "black", pch = 20, fill = "black"),
                
                ppolygons = list(border = "black", col = "grey", lty = 1, lwd = 1, alpha = 0.4),
                
                pSp = list(col = "grey", border = "black", lwd = 1, lty = 1, alpha = 1, cex = 3, pch = 20),
                
                psub = list(cex = 1, col = "black", position = "bottomleft", text = ""),
                
                ptable = list(x = list(srt = 0, pos = "top", tck = 5, adj = NA),
                              y = list(srt = 90, pos = "right", tck = 5, adj = NA),
                              margin = list(bottom = 5, left = 5, top = 5, right = 5))
                
                ),
           envir = .ADEgEnv)
    
    assign("adegtheme",
           list(layout.heights = list(
                    top.padding = 0, main.key.padding = 0,
                    key.axis.padding = 0,	axis.xlab.padding = 0,
                    xlab.key.padding = 0,	key.sub.padding = 0,
                    bottom.padding = 0),
                layout.widths = list(left.padding = 0, key.ylab.padding = 0, ylab.axis.padding = 0, axis.key.padding = 0, right.padding = 0),
                background = list(col = "transparent", alpha = 1),
                plot.polygon = list(col = "#F2F2F2"),
                plot.line = list(col = "#000000"),
                add.line = list(col =  "#000000", lty = 2),    
                ## clipping allows drawing to go outside panel (i.e : drawings) limits
                as.table = TRUE
                ), envir = .ADEgEnv
           )
    
    changelatticetheme(get("adegtheme", envir = .ADEgEnv))
    assign("currentadeg", list(), envir = .ADEgEnv)
}
