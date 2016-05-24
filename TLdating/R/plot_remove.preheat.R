#' Plotting of the preheat and TL curves
#'
#' This functions plots the results  obtained by mod_remove.preheat
#'
#' @param PH.signal
#'  \link{numeric}: matrix containing the preheat curves.
#' @param PH.temperatures
#'  \link{numeric}: matrix containing the temperature steps for each preheat curve.
#' @param PH.times
#'  \link{numeric}: matrix containing the time steps for each preheat curve.
#' @param TL.signal
#'  \link{numeric}: matrix containing the TL curves.
#' @param TL.temperatures
#'  \link{numeric}: matrix containing the temperature steps for each TL curve.
#'
#' @seealso
#'  \link{mod_remove.preheat}
#'
#' @author David Strebler
#'
## @export plot_remove.preheat

plot_remove.preheat <- function(

  PH.signal,
  PH.temperatures,
  PH.times,
  TL.signal,
  TL.temperatures
){
  # ------------------------------------------------------------------------------
  # Integrity Check
  # ------------------------------------------------------------------------------
  if (missing(PH.signal)){
    stop("[plot_remove.preheat] Error: Input 'PH.signal' is missing.")
  }else if (!is.numeric(PH.signal)){
    stop("[plot_remove.preheat] Error: Input 'PH.signal' is not of type 'numeric'.")
  }

  if (missing(PH.temperatures)){
    stop("[plot_remove.preheat] Error: Input 'PH.temperatures' is missing.")
  }else if (!is.numeric(PH.temperatures)){
    stop("[plot_remove.preheat] Error: Input 'PH.temperatures' is not of type 'numeric'.")
  }

  if (missing(PH.times)){
    stop("[plot_remove.preheat] Error: Input 'PH.times' is missing.")
  }else if (!is.numeric(PH.times)){
    stop("[plot_remove.preheat] Error: Input 'PH.times' is not of type 'numeric'.")
  }

  if (missing(TL.signal)){
    stop("[plot_remove.preheat] Error: Input 'TL.signal' is missing.")
  }else if (!is.numeric(TL.signal)){
    stop("[plot_remove.preheat] Error: Input 'TL.signal' is not of type 'numeric'.")
  }

  if (missing(TL.temperatures)){
    stop("[plot_remove.preheat] Error: Input 'TL.temperatures' is missing.")
  }else if (!is.numeric(TL.temperatures)){
    stop("[plot_remove.preheat] Error: Input 'TL.temperatures' is not of type 'numeric'.")
  }
  # ------------------------------------------------------------------------------

  # ------------------------------------------------------------------------------
  # Value check
  if(length(PH.signal) != length(PH.temperatures)){
    stop("[plot_remove.preheat] Error: PH.signal and PH.temperatures have a different size.")
  }else if(length(PH.signal) != length(PH.times)){
    stop("[plot_remove.preheat] Error: PH.signal and PH.times have a different size.")
  }
  if(length(TL.signal) != length(TL.temperatures)){
    stop("[plot_remove.preheat] Error: TL.signal and TL.temperatures have a different size.")
  }
  # ------------------------------------------------------------------------------

  old.par <- par( no.readonly = TRUE )
  par( oma = c(0.5, 0, 3, 0 ) )

  #Layout
  layout(matrix(c(1,1,2,2), 2, 2, byrow = TRUE))

  #Plot preheat

  if(length(PH.signal) > 0)
  {
    #Boundary
    plot.PH.Tmax <- max(PH.temperatures)
    plot.PH.Smax <- max(PH.times)
    plot.PH.Lmax <- max(PH.signal)

    #color
    colors <- 1:ncol(PH.signal)

    for(i in 1 : ncol(PH.signal)){
      temp.temperatures <- PH.temperatures[,i]
      temp.times <- PH.times[,i]
      temp.PH <- PH.signal[,i]
      temp.color <- colors[i]

      if(i == 1) {
        par(mar = c(5,5,4,5) )
        #Temperature
        plot(x=temp.times,
             y=temp.temperatures,
             xlim=c(0,plot.PH.Smax),
             ylim=c(0,plot.PH.Tmax),
             yaxt = "n",
             xaxt = "n",
             xlab = "",
             ylab = "",
             type="l",
             lty=2
        )
        axis(4)
        mtext(side = 4,
              text = "Temperature (\u00b0C)",
              line = 2.5,
              cex = 0.8
        )

        par(new = TRUE)

        plot(main= "Preheat signal",
             x=temp.times,
             y=temp.PH,
             xlim=c(0,plot.PH.Smax),
             ylim=c(0,plot.PH.Lmax),
             xlab="Time (s)",
             ylab = "Luminescence signal (PH)",
             type="l",
             col=temp.color
        )
        par(new = TRUE)

      }else{
        lines(x=temp.times,
              y=temp.PH,
              xlim=c(0,plot.PH.Smax),
              ylim=c(0,plot.PH.Lmax),
              col=temp.color
        )
      }
    }
    par(new = FALSE)
  }else{
    textplot(" ")
    title("Preheat signal")
  }

  #Plot TL
  if(length(TL.signal) > 0)
  {
    #Boundary
    plot.TL.Tmax <- max(TL.temperatures)
    plot.TL.Lmax <- max(TL.signal)

    #color
    colors <- 1:ncol(TL.signal)

    for(i in 1 : ncol(TL.signal)){
      temp.temperatures <- TL.temperatures[,i]
      temp.TL <- TL.signal[,i]
      temp.color <- colors[i]

      if(i == 1) {
        plot(main= "Thermoluminescence signal",
             x=temp.temperatures,
             y=temp.TL,
             xlim=c(0,plot.TL.Tmax),
             ylim=c(0,plot.TL.Lmax),
             xlab="Temperature (\u00b0C)",
             ylab = "Luminescence signal (TL)",
             type="l",
             col=temp.color
        )
        par(new = TRUE)

      }else{
        lines(x=temp.temperatures,
              y=temp.TL,
              xlim=c(0,plot.TL.Tmax),
              ylim=c(0,plot.TL.Lmax),
              col=temp.color
        )
      }
    }
    par(new = FALSE)
  }else{
    textplot(" ")
    title("Thermoluminescence signal")
  }

  #Page title ---------------------------------------------------------
  page.title <- paste("Preheat Removal")

  mtext(page.title, outer=TRUE,font = 2)

  #clean layout...
  layout(matrix(c(1), 1, 1, byrow = TRUE))
  par(old.par)

}
