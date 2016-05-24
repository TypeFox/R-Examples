#' Plotting function for mod_substract.background.
#'
#' This function plots the results of the \link{mod_substract.background} function.
#' It plots the TL curves, the background (BG) curves and the background substracted curves.
#'
#' @param old.TL
#'  \link{numeric}: Matrix containing the luminescence signal before background subtraction.
#' @param BG
#'  \link{numeric}: Matrix containing the luminescence signal from the background curves.
#' @param new.TL
#'  \link{numeric}: Matrix containing the luminescence signal after background subtraction.
#' @param temperatures
#'  \link{numeric}: Vector containing the temperature step
#'
#' @seealso
#'  \link{mod_substract.background}
#'
#' @author David Strebler
#'
## @export plot_substract.background

plot_substract.background <- function(

  old.TL,
  BG,
  new.TL,
  temperatures

){
  # ------------------------------------------------------------------------------
  # Integrity Check
  # ------------------------------------------------------------------------------
  if (missing(old.TL)){
    stop("[plot_substract.background] Error: Input 'old.TL' is missing.")
  }else if (!is.numeric(old.TL)){
    stop("[plot_substract.background] Error: Input 'old.TL' is not of type 'numeric'.")
  }

  if (missing(BG)){
    stop("[plot_substract.background] Error: Input 'BG' is missing.")
  }else if (!is.numeric(BG)){
    stop("[plot_substract.background] Error: Input 'BG' is not of type 'numeric'.")
  }

  if (missing(new.TL)){
    stop("[plot_substract.background] Error: Input 'new.TL' is missing.")
  }else if (!is.numeric(new.TL)){
    stop("[plot_substract.background] Error: Input 'new.TL' is not of type 'numeric'.")
  }

  if (missing(temperatures)){
    stop("[plot_substract.background] Error: Input 'temperatures' is missing.")
  }else if (!is.numeric(temperatures)){
    stop("[plot_substract.background] Error: Input 'temperatures' is not of type 'numeric'.")
  }
  # ------------------------------------------------------------------------------

  # ------------------------------------------------------------------------------
  # Value check
  if(length(old.TL) != length(BG)){
    stop("[plot_remove.preheat] Error: old.TL and BG have a different size.")
  }else if(length(old.TL) != length(new.TL)){
    stop("[plot_remove.preheat] Error: old.TL and new.TL have a different size.")
  }

  if(length(temperatures) != nrow(old.TL)){
    stop("[plot_remove.preheat] Error: temperatures and old.TL have a different size.")
  }
  # ------------------------------------------------------------------------------

  #Layout
  old.par <- par( no.readonly = TRUE )
  par( oma = c(0.5, 0, 3, 0 ) )
  layout(matrix(c(1,2,3,3), 2, 2, byrow = TRUE))

  #Plot TL (old.TL)
  #Boundary
  plot.Tmax <- max(temperatures)
  plot.TL.max <- max(old.TL)

  #color
  colors <- 1:ncol(old.TL)

  for(i in 1 : ncol(old.TL)){
    temp.TL <- old.TL[,i]
    temp.color <- colors[i]

    if(i == 1) {
      plot(main="TL before background substraction",
           x=temperatures,
           y=temp.TL,
           xlim=c(0,plot.Tmax),
           ylim=c(0,plot.TL.max),
           xlab="Temperature (C)",
           ylab="Luminescence signal (TL)",
           type="l",
           col=temp.color)

      par(new = TRUE)

    }else{
      lines(x=temperatures,
            y=temp.TL,
            xlim=c(0,plot.Tmax),
            ylim=c(0,plot.TL.max),
            col=temp.color)
    }
  }
  par(new = FALSE)

  #Plot BG
  #Boundary
  plot.Tmax <- max(temperatures)
  plot.TL.max <- max(old.TL)

  #color
  colors <- 1:ncol(BG)

  for(i in 1 : ncol(BG)){
    temp.TL <- BG[,i]
    temp.color <- colors[i]

    if(i == 1) {
      plot(main="Background",
           x=temperatures,
           y=temp.TL,
           xlim=c(0,plot.Tmax),
           ylim=c(0,plot.TL.max),
           xlab="Temperature (C)",
           ylab="Luminescence signal (BG)",
           type="l",
           col=temp.color)

      par(new = TRUE)

    }else{
      lines(x=temperatures,
            y=temp.TL,
            xlim=c(0,plot.Tmax),
            ylim=c(0,plot.TL.max),
            col=temp.color)
    }
  }
  par(new = FALSE)

  # Plot TL-BG (new.TL)
  #Boundary
  plot.Tmax <- max(temperatures)
  plot.TL.max <- max(new.TL)

  #color
  colors <- 1:ncol(new.TL)

  for(i in 1 : ncol(new.TL)){
    temp.TL <- new.TL[,i]
    temp.color <- colors[i]

    if(i == 1) {
      plot(main="TL after background substraction",
           x=temperatures,
           y=temp.TL,
           xlim=c(0,plot.Tmax),
           ylim=c(0,plot.TL.max),
           xlab="Temperature (C)",
           ylab="Luminescence signal (TL-BG)",
           type="l",
           col=temp.color)

      par(new = TRUE)

    }else{
      lines(x=temperatures,
            y=temp.TL,
            xlim=c(0,plot.Tmax),
            ylim=c(0,plot.TL.max),
            col=temp.color)
    }
  }
  par(new = FALSE)

  #Page title ---------------------------------------------------------
  page.title <- paste("Background substraction")

  mtext(page.title, outer=TRUE,font = 2)

  #clean layout...
  layout(matrix(c(1), 1, 1, byrow = TRUE))
  par(old.par)
}
