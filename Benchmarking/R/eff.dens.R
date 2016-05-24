# $Id: eff.dens.R 140 2015-05-15 21:48:02Z B002961 $

# Plot af taethed for efficiencer, bruger spejlingsprcincip

eff.dens <- function(eff, bw="nrd0")  {
   if ( class(eff) == "Farrell" )
      E <- eff$eff
   else
      E <- eff

   if ( max(E) <= 1 )
      orient <- "in"
   else if ( min(E) >= 1 )
      orient <- "out"
   else
      stop("Efficiencies should be on one side of 1 only, not on both sides")

   # Reflection around the boundary 1
   refl <- c(E,2-E)
   sr <- sort(refl)
   if ( orient=="in" ) {
      den <- density(sr, bw=bw, from=0, to=1, na.rm=TRUE)
   } else {
      den <- density(sr, bw=bw, from=1, na.rm=TRUE)
   }
   x <- den$x
   y <- 2*den$y
   eden <- list(x=x, y=y, bw=den$bw, n=den$n, data.name=den$data.name)
   class(eden) <- "density"
   return (eden)
} # eff.dens



eff.dens.plot  <- function(obj, bw="nrd0", ..., xlim, ylim, xlab, ylab)  {
   if ( class(obj) != "list" )  {
      o_ <- eff.dens(obj, bw=bw)
      obj <- o_
   }
   if ( missing(xlim) )  {
      if ( min(obj$x) < 1 ) {
         xl <- min(obj$x)
         xlim <- c(xl,1)
      } else {
         xr <- max(obj$x)
         xlim <- c(1,xr)
      }
   } 
   if ( missing(ylim) )  {
      ylim <- c(0, max(obj$y))
   }
   if ( missing(xlab) ) xlab <- "Efficiency"
   if ( missing(ylab) ) ylab <- "Density"
   plot(obj$x, obj$y, type="l",xlim=xlim,ylim=ylim,ylab=ylab,xlab=xlab,
               frame=FALSE,...)
} # eff.dens.plot

