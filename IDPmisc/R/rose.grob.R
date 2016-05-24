## rose.grob.R

rose.grob <- function(rho,
                      cyclVar,
                      circle,
                      transf,
                      general,
                      grid,
                      title,
                      gdat)
    ## Creating the grob rose without legend
    ## No testing is made for any argument!
    ## Author: Rene Locher
    ## Version: 2009-03-16

{

    ## for drawing data
    if (general$type=="s") {## creating segments of a circle
        nude.rose <- segments.circle(
             rho = transf(rho)-transf(grid$ray$lim[1]),
             shift = general$shift,
             circle = circle,
             ncp = general$ncp,
             gp = if (general$stacked) ## colored areas
                      gpar(col = if(general$lwd>=0) general$col else
                                                    "black",
                           fill = general$col,
                           lwd = general$lwd) else ## colored lines
                      gpar(col = general$col,
                           lwd = general$lwd,
                           lty = general$lty)
             )
    } else { ## drawing lines between observations
        x.dat <-
            as.vector(sweep(transf(rho)-transf(grid$ray$lim[1]),
              MARGIN=1,sin(2*pi*(cyclVar+general$shift)/circle),"*"))
        y.dat <-
            as.vector(sweep(transf(rho)-transf(grid$ray$lim[1]),
              MARGIN=1,cos(2*pi*(cyclVar+general$shift)/circle),"*"))
        id.dat <- rep(1:ncol(rho), rep(nrow(rho),ncol(rho)))

        nude.rose <-
            polygonGrob(name = "data",
                        x = x.dat,
                        y = y.dat,
                        id = id.dat,
                        default.units = "native",
                        gp = if (general$stacked) ## colored areas
                        gpar(col = "black", fill = general$col,
                             lwd = 1) else ## colored lines
                        gpar(col = general$col,
                             lwd = general$lwd,
                             lty = general$lty)
                        )
    }

    rose <-
        gTree(name = "rose",
              children = gList(
              nude.rose,

              if (grid$circ$sub$plot)
              circleGrob(name = "subcircles",
                         x = 0, y = 0,
                         r = grid$circ$sub$r[grid$circ$sub$r>0],
                         default.units = "native",
                         gp=gpar(col = grid$circ$sub$col,
                                 fill="transparent",
                                 lwd = grid$circ$sub$lwd)),

              circleGrob(name = "circles",
                         x = 0, y = 0,
                         r = grid$circ$r[grid$circ$r>0],
                         default.units = "native",
                         gp = gpar(col = grid$circ$col,
                                   fill="transparent",
                                   lwd = grid$circ$lwd)),

              textGrob(name = "circ.lab",
                       grid$circ$value,
                       x = gdat$circ$lab$x,
                       y = gdat$circ$lab$y,
                       just = c("center","center"),
                       default.units = "native",
                       gp = gpar(cex=grid$circ$cex,
                                 fill="transparent")),

              segmentsGrob(name = "rays",
                           x0 = 0, y0 = 0,
                           x1 = grid$circ$r[length(grid$circ$r)] *
                           sin(2*pi/grid$ray$n*(1:grid$ray$n)),
                           y1 = grid$circ$r[length(grid$circ$r)] *
                           cos(2*pi/grid$ray$n*(1:grid$ray$n)),
                           default.units = "native",
                           gp=gpar(col=grid$circ$col,
                                   fill="transparent",
                                   lwd=grid$circ$lwd)),

              textGrob(name = "cyclVar.lab",
                       label = grid$cyclVar$lab,
                       x = gdat$cyclVar$lab$x,
                       y = gdat$cyclVar$lab$y,
                       default.units = "native",
                       just = c("center","center"),
                       gp=gpar(cex = grid$cyclVar$cex,
                               fill="transparent")),

              if (!is.null(title$text))
              textGrob(name = "title",
                       label = title$text,
                       x = 0,
                       y = gdat$title$y,
                       default.units = "native",
                       just = c("center","bottom"),
                       gp = gpar(cex=title$cex,
                                 fill="transparent"))))
    return(rose)
} ## rose.grob
