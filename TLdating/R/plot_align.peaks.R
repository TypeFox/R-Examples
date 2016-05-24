#' Plots mod_alignPeaks results 
#' 
#' This function plots the results obtained by mod_alignPeaks.
#' 
#' @param temperatures
#'  \link{numeric}: Vector containing the temperature step
#' @param old.TL
#'  \link{numeric}: Matrix containing the luminescence signal before the peak alignment.
#' @param new.TL
#'  \link{numeric}: Matrix containing the luminescence signal after the peak alignment.
#' @param ref.TL
#'  \link{numeric}: Matrix containing the luminescence signal used as reference to define the peak position.
#' @param pos.peak
#'  \link{numeric}: Average peak position.
#' @param plotting.parameters
#'  \link{list} (with default): list containing the plotting parameters. See details.
#'  
#'@details
#'  \bold{Plotting parameters} \cr
#'  The plotting parameters are:  \cr
#'  \describe{
#'  \item{\code{plot.Tmin}}{
#'    \link{logical}: Minimum temperature which is plotted.}
#'  \item{\code{plot.Tmax}}{
#'    \link{logical}: Maximum temperature which is plotted.}
#' }
#' 
#' @seealso   
#'  \link{mod_align.peaks}
#'  
#' @author David Strebler
#' 
## @export plot_align.peaks

plot_align.peaks <- function(

  temperatures,
  old.TL,
  new.TL,
  ref.TL,
  pos.peak,
  
  plotting.parameters=list(plot.Tmin=0,
                           plot.Tmax=NA)
  
){
  # ------------------------------------------------------------------------------  
  # Integrity Check
  # ------------------------------------------------------------------------------  
  if (missing(temperatures)){
    stop("[plot_align.peaks] Error: Input 'temperatures' is missing.") 
  }else if (!is.numeric(temperatures)){
    stop("[plot_align.peaks] Error: Input 'temperatures' is not of type 'numeric'.")
  } 
  
  if (missing(old.TL)){
    stop("[plot_align.peaks] Error: Input 'old.TL' is missing.") 
  }else if (!is.numeric(old.TL)){
    stop("[plot_align.peaks] Error: Input 'old.TL' is not of type 'numeric'.")
  } 

  if (missing(new.TL)){
    stop("[plot_align.peaks] Error: Input 'new.TL' is missing.") 
  }else if (!is.numeric(new.TL)){
    stop("[plot_align.peaks] Error: Input 'new.TL' is not of type 'numeric'.")
  } 
  
  if (missing(ref.TL)){
    stop("[plot_align.peaks] Error: Input 'ref.TL' is missing.") 
  }else if (!is.numeric(ref.TL)){
    stop("[plot_align.peaks] Error: Input 'ref.TL' is not of type 'numeric'.")
  }
  
  if (missing(pos.peak)){
    stop("[plot_align.peaks] Error: Input 'pos.peak' is missing.") 
  }else if (!is.numeric(pos.peak)){
    stop("[plot_align.peaks] Error: Input 'pos.peak' is not of type 'numeric'.")
  }
  
  if(!is.list(plotting.parameters)){
    stop("[plot_align.peaks] Error: Input 'plotting.parameters' is not of type 'list'.")
  }
  # ------------------------------------------------------------------------------
  
  Tmax <- max(temperatures)
  nPoints <- length(temperatures)

  plot.Tmin <- plotting.parameters$plot.Tmin
  plot.Tmax <- plotting.parameters$plot.Tmax
  
  # Check Values -------------------
  
  # Plotting parameters 
  if(!is.numeric(plot.Tmin)){
    if(!is.finite(plot.Tmin)  || is.null(plot.Tmin)){
      plot.Tmin <- 0
    }else{
      stop("[plot_align.peaks] Error: plot.Tmin is not numeric.")    
    }
  }
  
  if(!is.numeric(plot.Tmax)){
    if(!is.finite(plot.Tmax) || is.null(plot.Tmax)){
      plot.Tmax <- Tmax
    }else{
      stop("[plot_align.peaks] Error: plot.Tmax is not numeric.")      
    }
  }
  
  if(plot.Tmin > plot.Tmax){
    stop("[plot_align.peaks] Error: plot.Tmin > plot.Tmax")
  }
  
  if(plot.Tmin < 0){
    plot.Tmin <- 0
  }
  
  if(plot.Tmax > Tmax){
    plot.Tmax <- Tmax
  }
  # -------------------------------
  Tstep <- Tmax/nPoints
  
  plot.min <- ceiling(plot.Tmin/Tstep)
  plot.max <-floor(plot.Tmax/Tstep)
  
  #----------------------------------------------------------------------------------------------
  #Plot results
  #----------------------------------------------------------------------------------------------

  #Layout
  old.par <- par( no.readonly = TRUE )
  par( oma = c(0.5, 0, 3, 0 ) )
  layout(matrix(c(1,2,3,3), 2, 2, byrow = TRUE))
  
  #Plot not aligned
  #Boundary
  
  plot.TL.max <- max(old.TL,na.rm = TRUE)
  
  #color
  colors <- 1:ncol(old.TL)
  
  for(i in 1 : ncol(old.TL)){  
    temp.TL <- old.TL[,i]
    temp.color <- colors[i]
    
    if(i == 1) {
      plot(x=temperatures, 
           y=temp.TL, 
           xlim=c(0,Tmax),
           ylim=c(0,plot.TL.max), 
           xlab="Temperature (\u00b0C)",
           ylab = "Luminescence signal",
           main="TL before peaks alignement",
           type="l", 
           col=temp.color)
      
      par(new = TRUE)
      
    }else{
      lines(x=temperatures, 
            y=temp.TL, 
            col=temp.color, 
            xlim=c(0,Tmax),
            ylim=c(0,plot.TL.max)
      )
    }
  }
  par(new = FALSE)
  
  #Plot Reference TL (testdose)
  #Boundary
  plot.TL.max <- max(ref.TL,na.rm = TRUE)
  
  #color
  colors <- 1:ncol(ref.TL)
  
  for(i in 1 : ncol(ref.TL)){  
    temp.TL <- ref.TL[,i]
    temp.color <- colors[i]
    
    if(i == 1) {
      plot(x=temperatures, 
           y=temp.TL, 
           xlim=c(0,Tmax),
           ylim=c(0,plot.TL.max), 
           xlab="Temperature (\u00b0C)",
           ylab = "Luminescence signal (Tx)",
           main="Peak position",
           type="l", 
           col=temp.color)
      
      par(new = TRUE)
      
    }else{
      lines(x=temperatures, 
            y=temp.TL, 
            col=temp.color, 
            xlim=c(0,Tmax),
            ylim=c(0,plot.TL.max)
      )
    }
  }    
  abline(v=pos.peak,col=2,lty=3)
  par(new = FALSE)
    
  
  #Plot aligned
  #Boundary
  plot.TL.max <- max(new.TL[plot.min:plot.max,],na.rm = TRUE)
  
  #color
  colors <- 1:ncol(new.TL)
  
  for(i in 1 : ncol(new.TL)){  
    temp.TL <- new.TL[,i]
    temp.color <- colors[i]
    
    if(i == 1) {
      plot(x=temperatures, 
           y=temp.TL, 
           xlim=c(plot.Tmin,plot.Tmax),
           ylim=c(0,plot.TL.max), 
           xlab="Temperature (\u00b0C)",
           ylab = "Luminescence signal (peaks shifted)",
           main="TL after peak alignment",
           type="l", 
           col=temp.color)
      par(new = TRUE)
      
    }else{
      lines(x=temperatures, 
            y=temp.TL, 
            col=temp.color, 
            xlim=c(plot.Tmin,plot.Tmax),
            ylim=c(0,plot.TL.max)
      )
    }
  }  
  par(new = FALSE)
  
  #Page title ---------------------------------------------------------
  page.title <- paste("Peak Alignment")
  
  mtext(page.title, outer=TRUE,font = 2)
  
  #clean layout...
  layout(1)
  par(old.par)
  
}
