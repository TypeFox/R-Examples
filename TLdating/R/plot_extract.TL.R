#' plot the TL curves
#'
#' This function plots the results obtained by mod_extract.TL.
#'
#'
#' @param temperatures
#'  \link{numeric}: matrix containing the temperature steps for each TL curve.
#' @param TL
#'  \link{numeric}: Matrix containing the luminescence signal for the TL curves.
#'
#' @seealso
#'  \link{mod_extract.TL}
#'
#' @author David Strebler
#'
## @export plot_extract.TL

plot_extract.TL <- function(

  temperatures,

  TL

){
  # ------------------------------------------------------------------------------
  # Integrity Check
  # ------------------------------------------------------------------------------
  if (missing(temperatures)){
    stop("[plot_extract.TL] Error: Input 'temperatures' is missing.")
  }else if (!is.list(temperatures)){
    stop("[plot_extract.TL] Error: Input 'temperatures' is not of type 'list'.")
  }

  if (missing(TL)){
    stop("[plot_extract.TL] Error: Input 'TL' is missing.")
  }else if (!is.list(TL)){
    stop("[plot_extract.TL] Error: Input 'TL' is not of type 'list'.")
  }
  # ------------------------------------------------------------------------------

  # ------------------------------------------------------------------------------
  # Value check
  if(length(TL) != length(temperatures)){
    stop("[plot_extract.TL] Error: TL and temperatures have a different size.")
  }
  # ------------------------------------------------------------------------------

  old.par <- par( no.readonly = TRUE )
  par( oma = c(0.5, 0, 3, 0 ) )

  #Layout
  layout(matrix(c(1), 1, 1, byrow = TRUE))

  #Plot TL
  plot.TL.Tmax <- max(unlist(temperatures), na.rm=TRUE)
  plot.TL.Lmax <- max(unlist(TL), na.rm=TRUE)

  for(i in 1 : length(TL)){
    temp.temperatures <- temperatures[[i]]
    temp.TL <- TL[[i]]

    if(i == 1) {
      plot(main= "Thermoluminescence signals",
           x=temp.temperatures,
           y=temp.TL,
           xlim=c(0,plot.TL.Tmax),
           ylim=c(0,plot.TL.Lmax),
           xlab="Temperature (\u00b0C)",
           ylab = "Luminescence signal",
           type="l",
           col=i
      )
      par(new = TRUE)

    }else{
      lines(x=temp.temperatures,
            y=temp.TL,
            xlim=c(0,plot.TL.Tmax),
            ylim=c(0,plot.TL.Lmax),
            col=i
      )
    }
  }
  par(new = FALSE)

  #Page title ---------------------------------------------------------
  page.title <- paste("TL curves")

  mtext(page.title, outer=TRUE,font = 2)

  #clean layout...
  layout(matrix(c(1), 1, 1, byrow = TRUE))
  par(old.par)
}
