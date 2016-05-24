########## R-function: set.plot.dflts ##########

# Sets plotting defaults (except for ylim)

# Last changed: 06 JAN 2005

set.plot.dflts <- function(object,plot.params,num.curves)
{
   # First treat curve plotting parameters.

   lin.present <- !is.null(object$info$lin)
   pen.present <- !is.null(object$info$pen)
   krige.present <- !is.null(object$info$krige)

   if (lin.present|pen.present)
   {
      x.vals <- NULL
      if (lin.present)
         x.vals <- cbind(x.vals,object$info$lin$x)
      if (pen.present)
         x.vals <- cbind(x.vals,object$info$pen$x)

      if (is.null(plot.params$plot.it))
         plot.params$plot.it <- rep(TRUE,num.curves)

      if (is.null(plot.params$bty))
         plot.params$bty <- rep("l",num.curves)

      if (is.null(plot.params$main))
         plot.params$main <- rep("",num.curves)

      if (is.null(plot.params$xlab))
         plot.params$xlab <- c(object$info$lin$name,object$info$pen$name)

      if (is.null(plot.params$ylab))
         plot.params$ylab <- rep("",num.curves)

      if (is.null(plot.params$xlim))
      {
         plot.params$xlim$lower <- apply(x.vals,2,min)
         plot.params$xlim$upper <- apply(x.vals,2,max)
      }

      if (is.null(plot.params$grid.size))
         plot.params$grid.size <- rep(101,num.curves)

      if (is.null(plot.params$lty))
         plot.params$lty <- rep(1,num.curves)

      if (is.null(plot.params$lwd))
         plot.params$lwd <- rep(2,num.curves)

      if (is.null(plot.params$col))
         plot.params$col <- rep("black",num.curves)

      if (is.null(plot.params$se.lty))
         plot.params$se.lty <- rep(2,num.curves)

      if (is.null(plot.params$se.lwd))
         plot.params$se.lwd <- rep(2,num.curves)

      if (is.null(plot.params$se.col))
         plot.params$se.col <- rep("black",num.curves)

      if (is.null(plot.params$shade.col))
         plot.params$shade.col <- rep("grey70",num.curves)

      if (is.null(plot.params$rug.col))
         plot.params$rug.col <- rep("black",num.curves)

      if (is.null(plot.params$jitter.rug))
         plot.params$jitter.rug <- FALSE

      if (is.null(plot.params$zero.line))
         plot.params$zero.line <- TRUE

   }

   # Second treat image plotting parameters.

   if (krige.present&plot.params$plot.image)
   {
      xvals <- object$info$krige$x

      if (is.null(plot.params$image.bg))
         plot.params$image.bg <- "white"

      if (is.null(plot.params$image.bty))
         plot.params$image.bty <- "l"

      if (is.null(plot.params$image.main))
         plot.params$image.main <- ""

      if (is.null(plot.params$image.xlab))
         plot.params$image.xlab <- object$info$krige$name[1]

      if (is.null(plot.params$image.ylab))
         plot.params$image.ylab <- object$info$krige$name[2]

      if (is.null(plot.params$image.zlab))
         plot.params$image.zlab <- ""

      if (is.null(plot.params$image.xlim))
         plot.params$image.xlim <- range(xvals[,1])

      if (is.null(plot.params$image.ylim))
         plot.params$image.ylim <- range(xvals[,2])

      if (is.null(plot.params$image.grid.size))
         plot.params$image.grid.size <- c(64,64)


      if (is.null(plot.params$bdry))
      {
         plot.params$bdry <- default.bdry(xvals[,1],xvals[,2])

         cat("\n\n Would you like to save this boundary information in a file? (y/n): ")

         ans <- readline()

         if (ans=="y")
         {
            cat("\n\n Enter filename for storage of boundary information: ")

            bdry.filename <- readline()

            write(t(plot.params$bdry),file=bdry.filename,ncolumns=2)
         }
      }

      if (is.null(plot.params$leg.loc))
      {
         image.xlim <- plot.params$image.xlim
         image.ylim <- plot.params$image.ylim
         plot.params$leg.loc <- c(0.45*image.xlim[1]+0.55*image.xlim[2],
                                   0.8*image.ylim[1]+0.2*image.ylim[2])
      }

      if (is.null(plot.params$leg.dim))
      {
         image.xlim <- plot.params$image.xlim
         image.ylim <- plot.params$image.ylim

         plot.params$leg.dim <- c(0.3*(image.xlim[2]-image.xlim[1]),
                                   0.1*(image.ylim[2]-image.ylim[1]))
      }

      if (is.null(plot.params$image.zlab.col))
      {
         plot.params$image.zlab.col <- "black"
      }
   }

   return(plot.params)
}

######### End of set.plot.dflts ########
