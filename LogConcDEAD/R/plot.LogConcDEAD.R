plot.LogConcDEAD <- function (x, uselog = FALSE, type = "ic", addp =
                              TRUE,
                              drawlabels = TRUE,
                              gridlen = 100,
                              g, marg, g.marg, main, xlab, ylab,...) 
{
  d <- ncol(x$x)

  ## one dimensional data
  if (d == 1) {
    y <- x$x[, 1]
    o <- order(y)
    type <- "l"
    if( missing( xlab ) ) xlab <- "X"
    if (uselog) {
      if( missing( ylab ) ) ylab <- "log density estimate"
      plot(y[o], x$logMLE[o], type = type, ylab = ylab, xlab = xlab, ...)
      if (addp == TRUE) points(y[o], x$logMLE[o])
    } else {
      if( missing( ylab ) ) ylab <- "density estimate"
      ysample <- seq(min(y),max(y),length.out = 500)
      logMLEsample <- dlcd(ysample,x,uselog=TRUE)
      plot(ysample, exp(logMLEsample), ylab = ylab,  xlab = xlab, type = type, ...)
      if (addp == TRUE) points(y[o], exp(x$logMLE[o]))
    }
  }



  ## marginals
  else if ( !missing(marg) ) {

    ##If we have a valid marginal
    if (is.element(marg, 1:d)) {

      ##Check whether already calculated and, if not, calculate
      if( missing( g.marg) ) g.marg <- interpmarglcd(x,marg=marg)
      if( missing( main ) ) main <- paste( "Marginal for X_",marg, sep="" )
      if( missing( xlab ) ) xlab <- paste("X", marg )
      if(uselog) {
        if( missing( ylab ) ) ylab <- "log estimated marginal density"
        plot(g.marg$xo, log(g.marg$marg), type = "l", xlab = xlab, ylab = ylab, main=main,...) }
      else {
        if( missing( ylab ) ) ylab <- "estimated marginal density"
        plot(g.marg$xo, g.marg$marg, type = "l", xlab = xlab, 
             ylab = ylab, main=main,...)
      }
    }    
    else stop(cat("Marginal should be one of 1, ...",d,"\n"))
  }
  else {

    if (d > 2) {
      stop("It is not possible to plot for d>2. Marginal estimates may be plotted by setting the parameter marg")
    }
    ##A rather ugly solution: hardwiring the color map from the package "colorspace"; may be produced using the command heat_hcl(128) from that package
    mycolors <- c("#E2E6BD","#E2E794","#E2E78C","#E2E686","#E2E581","#E2E47D","#E3E279","#E3E176","#E3E072","#E3DF6F","#E4DE6C","#E4DD69","#E4DC67","#E4DA64","#E5D961","#E5D85F","#E5D75C","#E5D65A","#E5D458","#E6D355","#E6D253","#E6D151","#E6CF4F","#E7CE4C","#E7CD4A","#E7CC48","#E7CB46","#E7C944","#E8C842","#E8C741","#E8C63F","#E8C43D","#E8C33B","#E8C23A","#E8C038","#E9BF37","#E9BE35","#E9BD34","#E9BB32","#E9BA31","#E9B930","#E9B82F","#E9B62E","#E9B52D","#E9B42C","#EAB22B","#EAB12A","#EAB02A","#EAAF29","#EAAD29","#EAAC28","#EAAB28","#EAA928","#EAA828","#EAA728","#EAA628","#EAA428","#E9A328","#E9A229","#E9A029","#E99F2A","#E99E2A","#E99C2B","#E99B2C","#E99A2C","#E9982D","#E9972E","#E9962F","#E89430","#E89331","#E89231","#E89132","#E88F33","#E88E34","#E78D35","#E78B37","#E78A38","#E78939","#E7873A","#E6863B","#E6853C","#E6833D","#E6823E","#E5803F","#E57F40","#E57E41","#E47C42","#E47B43","#E47A45","#E47846","#E37747","#E37648","#E37449","#E2734A","#E2714B","#E1704C","#E16F4D","#E16D4E","#E06C4F","#E06A50","#E06951","#DF6852","#DF6653","#DE6554","#DE6355","#DE6256","#DD6057","#DD5F58","#DC5E59","#DC5C5A","#DB5B5B","#DB595C","#DA585D","#DA565E","#D9555F","#D95360","#D85161","#D85061","#D74E62","#D74D63","#D64B64","#D64A65","#D54866","#D54667","#D44567","#D44368","#D34169","#D33F6A")
    y <- x$x
      z <- x$logMLE
      if ( missing( g) ) g <- interplcd(x, gridlen = gridlen)
      
      if(uselog) {
         if( missing( main ) ) main <- "Log density estimate"
      } else {
        g$z <- exp(g$z)
        g$z[is.na(g$z)] <- 0
        z <- exp(z)
        if( missing( main ) ) main <- "Density estimate"
      }
      if (type == "p") 
        persp(g, zlab = main, xlab = "X_1", 
                ylab = "X_2", ...)
      else if (type == "i") 
        image(g, col = mycolors, main=main,...)
      else if (type == "c") 
        contour(g, main=main, drawlabels=drawlabels, ...)
      else if (type == "ic") {
        image(g, col = mycolors, main=main, ...)
        contour(g, add = TRUE, drawlabels=drawlabels, ...)
      }
      else if (type == "r") {
        if (!require("rgl", quietly = TRUE)) 
          stop("you need to install the rgl package")
         zlim <- range(g$z[!is.na(g$z)])

        zcolors <- mycolors[(g$z - min(g$z,na.rm=TRUE))*128/diff(zlim) + 1]
        open3d()
        par3d(cex=0.8)
        if (!uselog) {
          surface3d(g$x, g$y, g$z, color = zcolors, back = "lines")
          decorate3d(xlim = range(y[, 1]), ylim = range(y[, 
                                             2]), zlim = zlim, zlab = "Density", ...)
          zscale <- diff(range(c(g$x,g$y)))/diff(zlim)*0.8
          par3d(scale=c(1,1,zscale))
        }
        else {
          decorate3d(xlim = range(y[, 1]), ylim = range(y[, 
                                             2]), zlim = zlim, zlab = "Log density", 
                     ...)
          rgl.surface(g$x, g$y, g$z, coords = c(1, 3, 2), 
                      color = zcolors, back = "lines")
          zscale <- diff(range(c(g$x,g$y)))/diff(zlim)*0.8
          par3d(scale=c(1,1,zscale))
        }
        if (addp) {
          plot3d(x$x[, 1], x$x[, 2], z, pch = 4, size = 2, add = TRUE, col = "black")
        }
         
      }
      else {
        stop("type should be one of r, p, i, c, or ic")
      }
      if (addp && (type == "i" || type == "c" || type == 
                   "ic")) 
        points(y, ...)
    }
  }
