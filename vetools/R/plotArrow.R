# Verified 1.3.18
plotArrow <- function(
        shape="",
        pos = 1,
        offset.arrow = c(0, 0),
        north.lwd = par()$lwd+2,
        north.col = par()$col,
        ...)
        {
                position = pos
                x = c(0,0,0.5)
                y = c(0,3.5,2)
                nx = c(0,0,1,1)-0.5
                ny = c(0,1,0,1)+0.4
                dots = list(...)
                cin = par()$cin
                cex = par()$cex
                col = par()$col
                lwd = par()$lwd
                if ( is.element("cex", names(dots)) == TRUE ) { cex = dots$cex}
                if ( is.element("col", names(dots)) == TRUE ) { col = dots$col}
                if ( is.element("lwd", names(dots)) == TRUE ) { lwd = dots$lwd}
                
                north.arrow = cbind(x, y) * cex * cin[1]
                north.text =  cbind(nx, ny) * cex * cin[1]
                xy.arrow = c(0, 0)
                
                if ( ! is.character(shape) ) {
                        xwidth = range(north.arrow[,1])[2]
                        yheight = range(north.arrow[,2])[2]
                        pos = matrix(get.shape.range(shape), ncol = 2)
                        # 4-----3
                        # |     |
                        # |     |
                        # 1-----2
                        x = range(pos[,1])
                        y = range(pos[,2])
                        xy.arrow = c(x[1], y[1])
                        if ( position == 2 ) { # South-East
                                xy.arrow = c(x[2]-xwidth, y[1])
                        } else if ( position == 3 ) { # Nort-East
                                xy.arrow = c(x[2]-xwidth, y[2]-yheight)
                                
                        } else if ( position == 4 ) { # North-West
                                xy.arrow = c(x[1], y[2]-yheight)
                        }  # position == 1  South-West 
                        # Do nothing!
                }

        lines(north.arrow[, 1] + xy.arrow[1] + offset.arrow[1], 
              north.arrow[, 2] + xy.arrow[2] + offset.arrow[2], 
              lwd = lwd, 
              col = col)
        lines(north.text[, 1] + xy.arrow[1] + offset.arrow[1], 
              north.text[, 2] + xy.arrow[2] + offset.arrow[2], 
              lwd = north.lwd, 
              col = north.col)
}