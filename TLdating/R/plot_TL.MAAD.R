#' plot MAAD result
#' 
#' This function plots the results for \link{analyse_TL.MAAD}.
#' The first page regroups all the information about the additive curves (names, doses, intensity vs. temperature and plateau test for Lx, Tx and Lx/Tx).
#' The second page regroups all the information about the regenerative curves (names, doses, intensity vs. temperature and plateau test for Lx, Tx and Lx/Tx).
#' The third page regroups all the information about the equivalent dose (dose plateau for the palaeodose and the sublinearity correction, growth curves, rejection criteria,...).
#' 
#' @param sample.name
#'  \link{character} (\bold{required}): Sample name.
#' @param temperatures
#'  \link{numeric} (\bold{required}): temperature vector 
#' @param eval.Tmin
#'  \link{integer} (\bold{required}): Temperature (°C) of the lower boundary for the signal integration.
#' @param eval.Tmax 
#'  \link{integer} (\bold{required}): Temperature (°C) of the upper boundary for the signal integration.
#' @param aNames
#'  \link{character} (\bold{required}): Name vector for the additive curves.
#' @param aDoses
#'  \link{numeric} (\bold{required}): Dose vector for the additive curves.
#' @param aLx
#'  \link{numeric} (\bold{required}): Lx matrix for the additive curves.
#' @param aTx
#'  \link{numeric} (\bold{required}): Tx matrix for the additive curves.
#' @param aLxTx
#'  \link{numeric} (\bold{required}): Lx/Tx matrix for the additive curves.
#' @param aLx.plateau
#'  \link{numeric} (\bold{required}): Ln/Lx matrix for the additive curves.
#' @param aTx.plateau
#'  \link{numeric} (\bold{required}): Ln/Tx matrix for the additive curves.
#' @param aLxTx.plateau
#'  \link{numeric} (\bold{required}): (Ln/Tn)/(Lx/Tx) matrix for the additive curves.
#' @param rNames
#'  \link{character} (\bold{required}): Name vector for the regenerative curves.
#' @param rDoses
#'  \link{numeric} (\bold{required}): Dose vector for the regenerative curves.
#' @param rLx
#'  \link{numeric} (\bold{required}): Lx matrix for the regenerative curves.
#' @param rTx
#'  \link{numeric} (\bold{required}): Tx matrix for the regenerative curves.
#' @param rLxTx
#'  \link{numeric} (\bold{required}): Lx/Tx matrix for the regenerative curves.
#' @param rLx.plateau
#'  \link{numeric} (\bold{required}): Ln/Lx matrix for the regenerative curves.
#' @param rTx.plateau
#'  \link{numeric} (\bold{required}): Tn/Tx matrix for the regenerative curves.
#' @param rLxTx.plateau
#'  \link{numeric} (\bold{required}): (Ln/Tn)/(Lx/Tx) matrix for the regenerative curves.
#' @param DP.Q.line
#'  \link{numeric} (\bold{required}): Vector containing the estimation of Q for each T° step.
#' @param DP.Q.line.error
#'  \link{numeric} (\bold{required}): Vector containing the uncertainty on the estimation of Q for each T° step.
#' @param GC.Q.line 
#'  \link{numeric} (\bold{required}): growth curve for Q
#' @param GC.Q.slope
#'  \link{numeric} (\bold{required}): growth curve parameters for Q
#' @param GC.Q.LxTx
#'  \link{numeric} (\bold{required}): Lx/Tx vector used for Q estimation using the growth curve approach.
#' @param GC.Q.LxTx.error
#'  \link{numeric} (\bold{required}): Error on the Lx/tx vector used for Q estimation using the growth curve approach.
#' @param DP.I.line
#'  \link{numeric} (\bold{required}): Vector containing I for each temperature step.
#' @param DP.I.line.error
#'  \link{numeric} (\bold{required}): Vector containing the uncertainty on I for each temperature step.
#' @param GC.I.line
#'  \link{numeric} (\bold{required}): growth curve for I
#' @param GC.I.slope
#'  \link{numeric} (\bold{required}): growth curve parameters for I.
#' @param GC.I.LxTx
#'  \link{numeric} (\bold{required}): Lx/tx vector used for I estimation using the growth curve approach.
#' @param GC.I.LxTx.error
#'  \link{numeric} (\bold{required}): Error on the Lx/tx vector used for I estimation using the growth curve approach.
#' @param Q.DP
#'  \link{numeric} (\bold{required}): Q estimation using the dose plateau approach
#' @param Q.DP.error
#'  \link{numeric} (\bold{required}): Uncertainty on the Q estimation using the dose plateau approach
#' @param Q.GC
#'  \link{numeric} (\bold{required}): Q estimation using the growth curve approach
#' @param Q.GC.error
#'  \link{numeric} (\bold{required}): Uncertainty on the Q estimation using the growth curve approach
#' @param I.DP
#'  \link{numeric} (\bold{required}): I estimation using the dose plateau approach
#' @param I.DP.error
#'  \link{numeric} (\bold{required}): Uncertainty on the I estimation using the dose plateau approach
#' @param I.GC
#'  \link{numeric} (\bold{required}): I estimation using the growth curve approach
#' @param I.GC.error
#'  \link{numeric} (\bold{required}): Uncertainty on the I estimation using the growth curve approach
#' @param De.GC,
#'  \link{numeric} (\bold{required}): ED (Q+I) estimation using the growth curve approach
#' @param De.GC.error,
#'  \link{numeric} (\bold{required}): Uncertainty on the ED (Q+I) estimation using the growth curve approach
#' @param De.DP,
#'  \link{numeric} (\bold{required}): ED (Q+I) estimation using the dose plateau approach
#' @param De.DP.error
#'  \link{numeric} (\bold{required}): Uncertainty on the ED (Q+I) estimation using the dose plateau approach
#' @param rejection.values
#'  \link{list} (\bold{required}): result of the rejection tests.
#' @param fitting.parameters
#'  \link{list} (with default): list containing the fitting parameters. See details.
#' @param plotting.parameters
#'  \link{list} (with default): list containing the plotting parameters. See details.
#'  
#' @details 
#' 
#' \bold{Fitting parameters} \cr
#' The fitting parameters are:  \cr
#' \describe{
#'  \item{\code{method}}{
#'    \link{character}: Fitting method (\code{LIN}, \code{EXP}, \code{EXP+LIN} or \code{EXP+EXP}).}
#'  \item{\code{fit.weighted}}{
#'    \link{logical}: If the fitting is weighted or not.}
#'  \item{\code{fit.use.slope}}{
#'    \link{logical}: If the slope of the Q growth curve is reused for the sublinearity correction.}
#'  \item{\code{fit.aDoses.min}}{
#'    \link{numeric}: Lowest additive dose used for the fitting.}
#'  \item{\code{fit.aDoses.max}}{
#'    \link{numeric}: Highest additive dose used for the fitting.}
#'  \item{\code{fit.rDoses.min}}{
#'    \link{numeric}: Lowest regenerative dose used for the fitting.}
#'  \item{\code{fit.rDoses.max}}{
#'    \link{numeric}: Highest regenerative dose used for the fitting.}
#' }
#' See also \link{analyse_TL.MAAD}, \link{calc_TL.MAAD.fit.Q} and \link{calc_TL.MAAD.fit.I}. \cr
#' 
#' \bold{Plotting parameters} \cr
#' The plotting parameters are:  \cr
#' \describe{
#'  \item{\code{plot.Tmin}}{
#'    \link{numeric}: Lower temperature plotted.} 
#'  \item{\code{plot.Tmax}}{
#'    \link{numeric}: Higher temperature plotted.}
#'  \item{\code{no.plot}}{
#'    \link{logical}: If \code{TRUE}, the results will not be plotted.}
#' }
#' See also \link{analyse_TL.MAAD}. \cr
#'
#' @seealso 
#'  \link{analyse_TL.MAAD},
#'  \link{calc_TL.MAAD.fit.Q},
#'  \link{calc_TL.MAAD.fit.I}.
#'  
#' @author David Strebler
#' 
## @export plot_TL.MAAD

plot_TL.MAAD <- function(

  sample.name,
  temperatures,
  eval.Tmin,
  eval.Tmax,
  aNames,
  aDoses,
  aLx,
  aTx,
  aLxTx,
  aLx.plateau,
  aTx.plateau,
  aLxTx.plateau,
  rNames,
  rDoses,
  rLx,
  rTx,
  rLxTx,
  rLx.plateau,
  rTx.plateau,
  rLxTx.plateau,
  DP.Q.line,
  DP.Q.line.error,
  GC.Q.line,
  GC.Q.slope,
  GC.Q.LxTx,
  GC.Q.LxTx.error,
  DP.I.line,
  DP.I.line.error,
  GC.I.line,
  GC.I.slope,
  GC.I.LxTx,
  GC.I.LxTx.error,
  Q.DP,
  Q.DP.error,
  Q.GC,
  Q.GC.error,
  I.DP,
  I.DP.error,
  I.GC,
  I.GC.error,
  De.GC,
  De.GC.error,
  De.DP,
  De.DP.error,
  rejection.values,
  fitting.parameters,
  plotting.parameters=list(plot.Tmin=0,
                           plot.Tmax=NA)
  
  ){ 
  # ------------------------------------------------------------------------------  
  # Integrity Check
  # ------------------------------------------------------------------------------  
  if (missing(sample.name)){
    stop("[plot_TL.MAAD] Error: Input 'sample.name' is missing.") 
  }else if (!is.character(sample.name)){
    stop("[plot_TL.MAAD] Error: Input 'sample.name' is not of type 'character'.")
  } 
  
  # ...
  
  # ------------------------------------------------------------------------------  
  
  Tmax <- max(temperatures)
  nPoints <- length(temperatures)
  Tstep <- Tmax/nPoints
  eval.min <- ceiling(eval.Tmin/Tstep)
  eval.max <-floor(eval.Tmax/Tstep)
  
  uDoses <- unique(c(aDoses,rDoses))
  
  fit.method <- fitting.parameters$fit.method
  fit.weighted <- fitting.parameters$fit.weighted
  fit.rDoses.min <- fitting.parameters$fit.rDoses.min
  fit.rDoses.max <- fitting.parameters$fit.rDoses.max
  
  plot.Tmin <- plotting.parameters$plot.Tmin
  plot.Tmax <- plotting.parameters$plot.Tmax
  
  # ------------------------------------------------------------------------------  
  # Values check
  
  # Plotting parameters 
  if(!is.numeric(plot.Tmin)){
    if(!is.finite(plot.Tmin)  || is.null(plot.Tmin)){
      plot.Tmin <- 0
    }else{
      stop("[plot_TL.MAAD] Error: plot.Tmin is not numeric.")    
    }
  }
  
  if(!is.numeric(plot.Tmax)){
    if(!is.finite(plot.Tmax) || is.null(plot.Tmax) ){
      plot.Tmax <- Tmax
    }else{
      stop("[plot_TL.MAAD] Error: plot.Tmax is not numeric.")      
    }
  }
  
  if(plot.Tmin > plot.Tmax){
    stop("[plot_TL.MAAD] Error: plot.Tmin > plot.Tmax")
  }
  
  if(plot.Tmin < 0){
    plot.Tmin <- 0
  }
  
  if(plot.Tmax > Tmax){
    plot.Tmax <- Tmax
  }
  
  # -------------------------------
  
  #----------------------------------------------------------------------------------------------------------------
  #Plot results
  #----------------------------------------------------------------------------------------------------------------
  old.par <- par( no.readonly = TRUE )
  par( oma = c(0.5, 0, 3, 0 ) )
  
  page <- 0
  
  #---------------------------------------------------------------------------
  #Page 1
  #---------------------------------------------------------------------------
  if(length(aLxTx) > 0){
    page <- page+1
    #Layout
    layout(matrix(c(1,2,3,4,5,6,7,7), 4, 2, byrow = TRUE),heights = c(2,2,2,1))
    
    #color
    ref_colors <- rainbow(n=length(aNames)-1)
    colors <- seq(length(aNames))
    names(colors) <- aNames
    colors[names(colors)=="N"] <- 1
    colors[names(colors)!="N"] <- ref_colors
    
    #Lx (additive)
    plot.aLx.max <- max(aLx,na.rm = TRUE)*1.1
    
    for(i in 1 : length(aDoses)){
      temp.curve <- aLx[,i]
      temp.name <- colnames(aLx)[i]
      temp.color <- colors[temp.name]
      
      if(i == 1) {
        plot(x=temperatures, 
             y=temp.curve, 
             type="l", 
             col=temp.color, 
             xlim=c(plot.Tmin,plot.Tmax),
             ylim=c(0,plot.aLx.max), 
             main = "Lx",
             xlab = "Temperature (\u00b0C)", 
             ylab = "Luminescence signal")
        par(new = TRUE)
        
      }else{
        lines(x=temperatures, 
              y=temp.curve, 
              col=temp.color,
              xlim=c(plot.Tmin,plot.Tmax),
              ylim=c(0,plot.aLx.max) 
        )
      }
    }
    abline(v=eval.Tmin,col=3,lty=3)
    abline(v=eval.Tmax,col=2,lty=3)
    
    par(new = FALSE)
    
    #Lx.plateau (additive)
    plot.aLx.plateau.max <- max(aLx.plateau[eval.min:eval.max,],na.rm = TRUE)*1.2
    
    for(i in 1 : ncol(aLx.plateau)){  
      temp.curve <- aLx.plateau[,i]
      temp.name <- colnames(aLx.plateau)[i]
      temp.color <- colors[temp.name]
      
      if(i == 1) {
        plot(x=temperatures, 
             y=temp.curve, 
             type="l", 
             col=temp.color, 
             xlim=c(plot.Tmin,plot.Tmax),
             ylim=c(0,plot.aLx.plateau.max), 
             main = "Plateau test (Lx)",
             xlab = "Temperature (\u00b0C)", 
             ylab = "Luminescence signal"
        )
        par(new = TRUE)
      }else{
        lines(x=temperatures, 
              y=temp.curve, 
              col=temp.color, 
              xlim=c(plot.Tmin,plot.Tmax),
              ylim=c(0,plot.aLx.plateau.max) 
        )
      }
    }
    abline(v=eval.Tmin,col=3,lty=3)
    abline(v=eval.Tmax,col=2,lty=3)
    
    par(new = FALSE)
    
    #Tx (additive)
    plot.aTx.max <- max(aTx,na.rm = TRUE)*1.1
    
    for(i in 1 : length(aDoses)){
      temp.curve <- aTx[,i]
      temp.name <- colnames(aTx)[i]
      temp.color <- colors[temp.name]
      
      if(i == 1) {
        plot(x=temperatures, 
             y=temp.curve, 
             type="l", 
             col=temp.color, 
             xlim=c(plot.Tmin,plot.Tmax),
             ylim=c(0,plot.aTx.max), 
             main = "Tx",
             xlab = "Temperature (\u00b0C)", 
             ylab = "Luminescence signal"
        )
        par(new = TRUE)
        
      }else{
        lines(x=temperatures, 
              y=temp.curve, 
              col=temp.color, 
              xlim=c(plot.Tmin,plot.Tmax),
              ylim=c(0,plot.aTx.max) 
        )
      }
    }
    abline(v=eval.Tmin,col=3,lty=3)
    abline(v=eval.Tmax,col=2,lty=3)
    
    par(new = FALSE)
    
    #Tx.plateau (additive)
    plot.aTx.plateau.max <- max(aTx.plateau[eval.min:eval.max,],na.rm = TRUE)*1.2
    
    for(i in 1 : ncol(aTx.plateau)){   
      temp.curve <- aTx.plateau[,i]
      temp.name <- colnames(aTx.plateau)[i]
      temp.color <- colors[temp.name]
      
      if(i == 1) {
        plot(x=temperatures, 
             y=temp.curve, 
             type="l", 
             col=temp.color, 
             xlim=c(plot.Tmin,plot.Tmax),
             ylim=c(0,plot.aTx.plateau.max), 
             main = "Plateau test (Tx)",
             xlab = "Temperature (\u00b0C)", 
             ylab = "Luminescence signal"
        )
        par(new = TRUE)
        
      }else{
        lines(x=temperatures, 
              y=temp.curve, 
              col=temp.color, 
              xlim=c(plot.Tmin,plot.Tmax),
              ylim=c(0,plot.aTx.plateau.max) 
        )
      }
    }
    abline(v=eval.Tmin,col=3,lty=3)
    abline(v=eval.Tmax,col=2,lty=3)
    
    par(new = FALSE)
    
    #LxTx (additive)
    plot.aLxTx.max <- max(aLxTx[eval.min:eval.max,],na.rm = TRUE)*1.1
    
    for(i in 1 : length(aDoses)){
      temp.curve <- aLxTx[,i]
      temp.name <- colnames(aLxTx)[i]
      temp.color <- colors[temp.name]
      
      if(i == 1) {
        plot(x=temperatures, 
             y=temp.curve, 
             type="l", 
             col=temp.color, 
             xlim=c(plot.Tmin,plot.Tmax),
             ylim=c(0,plot.aLxTx.max), 
             main = "Lx/Tx",
             xlab = "Temperature (\u00b0C)", 
             ylab = "Signal"
        )
        par(new = TRUE)
        
      }else{
        lines(x=temperatures, 
              y=temp.curve, 
              col=temp.color, 
              xlim=c(plot.Tmin,plot.Tmax),
              ylim=c(0,plot.aLxTx.max) 
        )
      }
    }
    abline(v=eval.Tmin,col=3,lty=3)
    abline(v=eval.Tmax,col=2,lty=3)
    
    par(new = FALSE)
    
    #LxTx.plateau (additive)
    plot.aLxTx.plateau.max <- max(aLxTx.plateau[eval.min:eval.max,],na.rm = TRUE)*1.2
    
    for(i in 1 : ncol(aLxTx.plateau)){
      temp.curve <- aLxTx.plateau[,i]
      temp.name <- colnames(aLxTx.plateau)[i]
      temp.color <- colors[temp.name]
      
      
      if(i == 1) {
        plot(x=temperatures, 
             y=temp.curve, 
             type="l", 
             col=temp.color, 
             xlim=c(plot.Tmin, plot.Tmax),
             ylim=c(0, plot.aLxTx.plateau.max), 
             main = "Plateau test (Lx/Tx)",
             xlab = "Temperature (\u00b0C)", 
             ylab = "Signal"
        )
        par(new = TRUE)
        
      }else{
        lines(x=temperatures, 
              y=temp.curve, 
              col=temp.color, 
              xlim=c(plot.Tmin,plot.Tmax),
              ylim=c(0,plot.aLxTx.plateau.max) 
        )
      }
    }
    abline(v=eval.Tmin,col=3,lty=3)
    abline(v=eval.Tmax,col=2,lty=3)
    
    par(new = FALSE)
    
    #Legend
    legend.size <- length(aDoses)
    legend.text <- c(aNames, aDoses)
    legend.color <- matrix(data = colors,
                           nrow = 2,
                           ncol = legend.size,
                           byrow=TRUE
    )
    legend <- matrix(data = legend.text, 
                     nrow = 2,
                     ncol = legend.size,
                     byrow=TRUE,
                     dimnames = list(c("Names", "Doses (s)"),
                                     vector(mode = "character",length = legend.size)
                     )
    )
    
    textplot(object= legend, 
             col.data = legend.color,
             cex = 1.5,
             halign = "center",
             valign="center",
             show.colnames= FALSE,
             show.rownames= TRUE
    )
    
    #Page title
    page.title <- paste("MAAD: ",
                        sample.name,
                        " - page ",
                        page,
                        ": Additive doses",
                        sep = "")
    mtext(page.title, outer=TRUE,font = 2)
  }
  
  #---------------------------------------------------------------------------
  #Page 2
  #---------------------------------------------------------------------------
  
  if(length(rLxTx) > 0){
    
    page <- page+1
    
    #Layout
    layout(matrix(c(1,2,3,4,5,6,7,7), 4, 2, byrow = TRUE),heights = c(2,2,2,1))
    
    #color
    ref_colors <- rainbow(n=length(rNames)-1)
    colors <- seq(length(rNames))
    names(colors) <- rNames
    colors[names(colors)=="R0"] <- 1
    colors[names(colors)!="R0"] <- ref_colors
    
    #Lx (regenerative)
    plot.rLx.max <- max(rLx,na.rm = TRUE)*1.1
    
    for(i in 1 : length(rDoses)){
      temp.curve <- rLx[,i]
      temp.name <- colnames(rLx)[i]
      temp.color <- colors[temp.name]
      
      if(i == 1) {
        plot(x=temperatures, 
             y=temp.curve, 
             type="l", 
             col=temp.color, 
             xlim=c(plot.Tmin,plot.Tmax),
             ylim=c(0,plot.rLx.max), 
             main = "Lx (regenerative curve)",
             xlab = "Temperature (\u00b0C)", 
             ylab = "Luminescence signal"
        )
        par(new = TRUE)
        
      }else{
        lines(x=temperatures, 
              y=temp.curve, 
              col=temp.color,
              xlim=c(plot.Tmin,plot.Tmax),
              ylim=c(0,plot.rLx.max) 
        )
      }
    }
    abline(v=eval.Tmin,col=3,lty=3)
    abline(v=eval.Tmax,col=2,lty=3)
    
    par(new = FALSE)
    
    #Lx.plateau (regenerative)
    plot.rLx.plateau.max <- max(rLx.plateau[eval.min:eval.max,],na.rm = TRUE)*1.2
    
    for(i in 1 : ncol(rLx.plateau)){
      temp.curve <- rLx.plateau[,i]
      temp.name <- colnames(rLx.plateau)[i]
      temp.color <- colors[temp.name]
      
      if(i == 1) {
        plot(x=temperatures, 
             y=temp.curve, 
             type="l", 
             col=temp.color, 
             xlim=c(plot.Tmin, plot.Tmax),
             ylim=c(0, plot.rLx.plateau.max), 
             main = "Plateau test (Lx)",
             xlab = "Temperature (\u00b0C)", 
             ylab = "Signal"
        )
        par(new = TRUE)
        
      }else{
        lines(x=temperatures, 
              y=temp.curve, 
              col=temp.color, 
              xlim=c(plot.Tmin,plot.Tmax),
              ylim=c(0,plot.rLx.plateau.max) 
        )
      }
    }
    abline(v=eval.Tmin,col=3,lty=3)
    abline(v=eval.Tmax,col=2,lty=3)
    
    par(new = FALSE)
    
    #Tx (regenerative)
    plot.rTx.max <- max(rTx,na.rm = TRUE)*1.1
    
    for(i in 1 : length(rDoses)){        
      temp.curve <- rTx[,i]
      temp.name <- colnames(rTx)[i]
      temp.color <- colors[temp.name]
      
      if(i == 1) {
        plot(x=temperatures, 
             y=temp.curve, 
             type="l", 
             col=temp.color, 
             xlim=c(plot.Tmin,plot.Tmax),
             ylim=c(0,plot.rTx.max), 
             main = "Tx (regenerative curve)",
             xlab = "Temperature (\u00b0C)", 
             ylab = "Luminescence signal"
        )
        par(new = TRUE)
        
      }else{
        lines(x=temperatures, 
              y=temp.curve, 
              col=temp.color, 
              xlim=c(plot.Tmin,plot.Tmax),
              ylim=c(0,plot.rTx.max) 
        )
      }
    }
    abline(v=eval.Tmin,col=3,lty=3)
    abline(v=eval.Tmax,col=2,lty=3)
    
    par(new = FALSE)
    
    #Tx.plateau (regenerative)
    plot.rTx.plateau.max <- max(rTx.plateau[eval.min:eval.max,],na.rm = TRUE)*1.2
    
    for(i in 1 : ncol(rTx.plateau)){   
      temp.curve <- rTx.plateau[,i]
      temp.name <- colnames(rTx.plateau)[i]
      temp.color <- colors[temp.name]
      
      if(i == 1) {
        plot(x=temperatures, 
             y=temp.curve, 
             type="l", 
             col=temp.color, 
             xlim=c(plot.Tmin, plot.Tmax),
             ylim=c(0, plot.rTx.plateau.max), 
             main = "Plateau test (Lx)",
             xlab = "Temperature (\u00b0C)", 
             ylab = "Signal"
        )
        par(new = TRUE)
        
      }else{
        lines(x=temperatures, 
              y=temp.curve, 
              col=temp.color, 
              xlim=c(plot.Tmin,plot.Tmax),
              ylim=c(0,plot.rTx.plateau.max) 
        )
      }
    }
    abline(v=eval.Tmin,col=3,lty=3)
    abline(v=eval.Tmax,col=2,lty=3)
    
    par(new = FALSE)
    
    #LxTx (Regenerative)
    plot.rLxTx.max <- max(rLxTx[eval.min:eval.max,],na.rm = TRUE)*1.1
    
    for(i in 1 : length(rDoses)){  
      temp.curve <- rLxTx[,i]
      temp.name <- colnames(rLxTx)[i]
      temp.color <- colors[temp.name]
      
      if(i == 1) {
        plot(x=temperatures, 
             y=temp.curve, 
             type="l", 
             col=temp.color, 
             xlim=c(plot.Tmin,plot.Tmax),
             ylim=c(0,plot.rLxTx.max), 
             main = "Lx/Tx (regenerative curve)",
             xlab = "Temperature (\u00b0C)", 
             ylab = "Signal"
        )
        par(new = TRUE)
        
      }else{
        lines(x=temperatures, 
              y=temp.curve, 
              col=temp.color, 
              xlim=c(plot.Tmin,plot.Tmax),
              ylim=c(0,plot.rLxTx.max) 
        )
      }
    }
    abline(v=eval.Tmin,col=3,lty=3)
    abline(v=eval.Tmax,col=2,lty=3)
    
    par(new = FALSE)
    
    #LxTx.plateau (regenerative)
    plot.rLxTx.plateau.max <- max(rLxTx.plateau[eval.min:eval.max,],na.rm = TRUE)*1.2
    
    for(i in 1 : ncol(rLxTx.plateau)){ 
      temp.curve <- rLxTx.plateau[,i]
      temp.name <- colnames(rLxTx.plateau)[i]
      temp.color <- colors[temp.name]
      
      if(i == 1) {
        plot(x=temperatures, 
             y=temp.curve, 
             type="l", 
             col=temp.color, 
             xlim=c(plot.Tmin, plot.Tmax),
             ylim=c(0, plot.rLxTx.plateau.max), 
             main = "Plateau test (Lx)",
             xlab = "Temperature (\u00b0C)", 
             ylab = "Signal"
        )
        par(new = TRUE)
        
      }else{
        lines(x=temperatures, 
              y=temp.curve, 
              col=temp.color, 
              xlim=c(plot.Tmin,plot.Tmax),
              ylim=c(0,plot.rLxTx.plateau.max) 
        )
      }
    }
    abline(v=eval.Tmin,col=3,lty=3)
    abline(v=eval.Tmax,col=2,lty=3)
    
    par(new = FALSE)
    
    #Legend
    legend.size <- length(rDoses)
    legend.text <- c(rNames, rDoses)
    legend.color <- matrix(data = colors,
                           nrow = 2,
                           ncol = legend.size,
                           byrow=TRUE
    )
    legend <- matrix(data = legend.text, 
                     nrow = 2,
                     ncol = legend.size,
                     byrow=TRUE,
                     dimnames = list(c("Names", "Doses (s)"),
                                     vector(mode = "character",length = legend.size)
                     )
    )
    
    textplot(object= legend, 
             col.data = legend.color,
             cex = 1.5,
             halign = "center",
             valign="center",
             show.colnames= FALSE,
             show.rownames= TRUE
    )
    
    #Page title
    page.title <- paste("MAAD: ",
                        sample.name,
                        " - page ",
                        page,
                        ": Regenerative doses",
                        sep = "")
    mtext(page.title, outer=TRUE,font = 2)
  }
  
  #---------------------------------------------------------------------------
  #Page 3
  #---------------------------------------------------------------------------
  
  page <- page+1
  
  #Layout
  layout(matrix(c(1,1,2,1,1,4,3,3,5,3,3,6), 4, 3, byrow = TRUE))
  
  # Plotting  Palaeodose (Q) ----------------------------------------
  
  if(length(aLxTx) > 0){
    
    plot.DP.Q.line.max <- max(DP.Q.line[eval.min:eval.max],na.rm = TRUE)*1.5
    
    plot(x=temperatures, 
         y=DP.Q.line, 
         xlim=c(plot.Tmin, plot.Tmax),
         ylim=c(0, plot.DP.Q.line.max), 
         xlab = "Temperature (\u00b0C)", 
         ylab = "Dose (s)",
         main = "D\u2091 plateau - Palaeodose (Q)",
         sub = paste("Q =", 
                     round(Q.DP, digits = 2), "\u00b1", round(Q.DP.error, digits = 2),
                     paste( "(", round(Q.DP.error/Q.DP*100, digits = 2), "%)",sep = "")
         ),
         type="b",
         lty=2,
         pch=18,
         col=6)
    
    par(new = TRUE)
    
    abline(v=eval.Tmin,col=3,lty=3)
    abline(v=eval.Tmax,col=2,lty=3)
    
    arrows(temperatures, 
           DP.Q.line-DP.Q.line.error, 
           temperatures, 
           DP.Q.line+DP.Q.line.error, 
           length=0.05, 
           angle=90, 
           code=3)
    
    par(new = FALSE)
    
  }else if(length(rLxTx) > 0){
    
    plot.DP.I.line.max <- max(DP.I.line[eval.min:eval.max],na.rm = TRUE)*2
    
    plot(x=temperatures, 
         y=DP.I.line, 
         xlim=c(plot.Tmin, plot.Tmax),
         ylim=c(0, plot.DP.I.line.max), 
         main = "D\u2091 plateau - Sublinearity corr. (I)",
         sub = paste("I =", 
                     round(I.DP, digits = 2), "\u00b1", round(I.DP.error, digits = 2),
                     paste( "(", round(I.DP.error/I.DP*100, digits = 2), "%)",sep = "")
         ),
         xlab = "Temperature (\u00b0C)", 
         ylab = "Dose (s)",
         type="b",
         lty=2,
         pch=18,
         col=4)
    
    par(new = TRUE)
    
    abline(v=eval.Tmin,col=3,lty=3)
    abline(v=eval.Tmax,col=2,lty=3)
    
    arrows(temperatures, 
           DP.I.line-DP.I.line.error, 
           temperatures, 
           DP.I.line+DP.I.line.error, 
           length=0.05, 
           angle=90, 
           code=3)
    
    par(new = FALSE)
    
  }else{
    #Empty space
    textplot(" ")
    title("Palaeodose (Q)")
  }
  
  # Plotting sublinerarity correction (I) -----------------------------------------
  
  if(length(rLxTx) > 0 && length(aLxTx) > 0){
    
    plot.DP.I.line.max <- max(DP.I.line[eval.min:eval.max],na.rm = TRUE)*2
    
    plot(x=temperatures, 
         y=DP.I.line,
         xlim=c(plot.Tmin, plot.Tmax),
         ylim=c(0, plot.DP.I.line.max), 
         main = "D\u2091 plateau - Sublinearity corr. (I)",
         sub = paste("I =", 
                     round(I.DP, digits = 2), "\u00b1", round(I.DP.error, digits = 2),
                     paste( "(", round(I.DP.error/I.DP*100, digits = 2), "%)",sep = "")
         ),
         xlab = "Temperature (\u00b0C)", 
         ylab = "Dose (s)",
         type="l",
         pch=18,
         lty=1,
         col=4)
    
    par(new = TRUE)
    
    abline(v=eval.Tmin,col=3,lty=3)
    abline(v=eval.Tmax,col=2,lty=3)
    
    par(new = FALSE)
    
  }else if(length(aLxTx) > 0){
    #Empty space
    textplot(" ")
    title("D\u2091 plateau - Sublinearity corr. (I)")
  }else{
    #Empty space
    textplot(" ")
  }
  
  # Plotting  GC (Q&I) ----------------------------------------
  
  # ylim & xlim 
  
  if(length(aLx) == 0 && length(rLx) == 0){
    plot.GC.ymax <- 0
  }else if(length(aLx) == 0){
    plot.GC.ymax <- max(GC.Q.LxTx+GC.Q.LxTx.error,na.rm = TRUE)
  }else if(length(rLx) == 0){
    plot.GC.ymax <- max(GC.I.LxTx+GC.I.LxTx.error,na.rm = TRUE)
  }else{
    plot.GC.ymax <- max(GC.I.LxTx+GC.I.LxTx.error, GC.Q.LxTx+GC.Q.LxTx.error, na.rm = TRUE)
  }
  
  plot.GC.xmin <- -1.1*(De.DP+De.DP.error)
  plot.GC.xmax <- max(uDoses,na.rm = TRUE)    
  
  # Additive curve
  
  plot(x=NA, 
       y=NA,
       xlim = c(plot.GC.xmin,plot.GC.xmax),
       ylim = c(0,plot.GC.ymax), 
       main = "Growth curves",
       sub = paste(if(length(aLx) > 0){paste("Q (GC) =",
                                               paste(round(Q.GC, digits = 2), 
                                                     "\u00b1", 
                                                     round(Q.GC.error, digits = 2)),
                                               paste("(",
                                                     round(Q.GC.error/Q.GC*100, digits = 2), 
                                                     "%)",
                                                     sep = ""))},
                   if(length(aLx) > 0 && length(rLx) > 0){"|"},
                   if(length(rLx) > 0){paste("I (GC) = ",
                                               paste(round(I.GC, digits = 2), 
                                                     "\u00b1", 
                                                     round(I.GC.error, digits = 2)),
                                               paste("(", 
                                                     round(I.GC.error/I.GC*100, digits = 2), 
                                                     "%)",
                                                     sep = ""))}),
       xlab = "Dose (s)",
       ylab = "Signal (Lx/Tx)", 
       type = "p",
       pch = 18,
       col = 1)
  
  par(new = TRUE)
  
  if(length(GC.I.LxTx)>0){
    
    lines(x = aDoses[aNames!="N"], 
          y = GC.I.LxTx[aNames!="N"],
          type="p",
          pch=18,
          col=1)
    
    par(new = TRUE)
    
    # Natural
    points(x = aDoses[aNames=="N"],
           y = GC.I.LxTx[aNames=="N"],
           pch = 5,
           col = 1)
    
    # error on aLxTx
    arrows(aDoses, 
           GC.I.LxTx-GC.I.LxTx.error, 
           aDoses, 
           GC.I.LxTx+GC.I.LxTx.error, 
           length=0, #0.05, 
           angle=90, 
           code=3)
    
    # linear regression
    if(length(GC.Q.line)>0){
      abline(GC.Q.line)
      
      # Q.GC
      points(x=-Q.GC,
             y=0,
             pch=18,
             col=3)
      # error on Q.GC
      arrows(-Q.GC-Q.GC.error, 
             0, 
             -Q.GC+Q.GC.error, 
             0, 
             length=0.05, 
             angle=90, 
             code=3)
      
      # Q.DP
      points(x=-Q.DP,
             y=0,
             pch=18,
             col=2)
    }
    
  }else{
    plot(x=NA, 
         y=NA,
         xlim = c(plot.GC.xmin,plot.GC.xmax),
         ylim = c(0,plot.GC.ymax), 
         main = "Palaeodose (s)",
         xlab = "Dose (s)",
         ylab = "Signal (Lx/Tx)", 
         type = "p",
         pch = 18,
         col = 1)
    
    par(new = TRUE)
  }
  
  # Regenerative curve
  if(length(GC.Q.LxTx)>0){
    # rLxTx    
    lines(x = rDoses[rDoses < fit.rDoses.min | rDoses > fit.rDoses.max],
          y = GC.Q.LxTx[rDoses < fit.rDoses.min | rDoses > fit.rDoses.max],
          type = "p",
          pch = 5,
          col = 4) 
    
    lines(x = rDoses[rDoses >= fit.rDoses.min & rDoses <= fit.rDoses.max],
          y = GC.Q.LxTx[rDoses >= fit.rDoses.min & rDoses <= fit.rDoses.max],
          type = "p",
          pch = 18,
          col = 4) 
    
    lines(x = rDoses,
          y = GC.Q.LxTx,
          type = "l",
          pch = 5,
          col = 4)
    
    #error on rLxTx
    arrows(rDoses, 
           GC.Q.LxTx-GC.Q.LxTx.error, 
           rDoses, 
           GC.Q.LxTx+GC.Q.LxTx.error, 
           length=0, #0.05, 
           angle=90, 
           code=3)
    
    # linear regression
    if(length(GC.I.line)>0){
      abline(GC.I.line)
      
      # I.GC
      points(x = I.GC,
             y = 0,
             pch = 18,
             col = 3)  
      # error on I.GC
      arrows(I.GC-I.GC.error, 
             0, 
             I.GC+I.GC.error, 
             0, 
             length = 0.05, 
             angle = 90, 
             code = 3)
      
      # I
      points(x = I.DP,
             y = 0,
             pch = 18,
             col = 2)  
    }
  }
  
  # Legend
  legend(x = "topleft", 
         legend = c(if(length(GC.I.LxTx)>0){c("Natural", "Natural + \u03b2")}, 
                    if(length(GC.Q.LxTx)>0){c("REG points (not used)", "REG points (used)")}, 
                    if(length(GC.I.LxTx)>0){c("Q (DP)", "Q (GC)")}, 
                    if(length(GC.Q.LxTx)>0){c("I (DP)", "I (GC)")}),
         pch = c(if(length(GC.I.LxTx)>0){c(5, 18)},
                 if(length(GC.Q.LxTx)>0){c(5, 18)},
                 if(length(GC.I.LxTx)>0){c(18, 18)},
                 if(length(GC.Q.LxTx)>0){c(18, 18)}),
         col = c(if(length(GC.I.LxTx)>0){c(1,1)},
                 if(length(GC.Q.LxTx)>0){c(4,4)},
                 if(length(GC.I.LxTx)>0){c(2,3)},
                 if(length(GC.Q.LxTx)>0){c(2,3)}),
         bty = "n")
  
  par(new = FALSE)
  
  
  #Rejection criteria ----------------------------------------
  rejection.title <- "Rejection criteria"
  
  rejection.text <- c("Q:",
                      "",
                      "Lx error (max):",
                      paste(round(rejection.values$aLx.error.max*100,digits = 2),"%",sep = ""),                                  
                      "Tx error (max):",
                      paste(round(rejection.values$aTx.error.max*100,digits = 2),"%",sep = ""),
                      "I:",
                      "",
                      "Lx error (max):",
                      paste(round(rejection.values$rLx.error.max*100,digits = 2),"%",sep = ""),                                  
                      "Tx error (max):",
                      paste(round(rejection.values$rTx.error.max*100,digits = 2),"%",sep = ""))
  
  rejection <- matrix(data = rejection.text,
                      nrow = 6,
                      ncol = 2,
                      byrow=TRUE)
  
  rejection.color <- matrix(data = c(6,6,
                                     if(rejection.values$test.aLx.error){1}else{2}, if(rejection.values$test.aLx.error){1}else{2}, 
                                     if(rejection.values$test.aTx.error){1}else{2}, if(rejection.values$test.aTx.error){1}else{2},
                                     4,4,
                                     if(rejection.values$test.rLx.error){1}else{2}, if(rejection.values$test.rLx.error){1}else{2}, 
                                     if(rejection.values$test.rTx.error){1}else{2}, if(rejection.values$test.rTx.error){1}else{2}),
                            nrow = 6,
                            ncol = 2,
                            byrow=TRUE)
  
  textplot(object=rejection, 
           col.data=rejection.color,
           cex=1.2,
           halign="center",
           valign="top",
           show.colnames= FALSE,
           show.rownames= FALSE)
  
  title(rejection.title)
  
  # Curve fitting -----------------------------------------------------
  
  if(fit.method == "LIN"){
    fitting.title <- paste("Curve fitting (GC):",
                           "Linear",
                           if(fit.weighted){"(weighted)"},
                           "\n",
                           "y = a + bx")
    
    if(length(GC.Q.line)>0){
      fitting.text <- c("a (Q) =",
                        paste(format(GC.Q.slope$a, digits = 3, scientific = TRUE), "\u00b1", format(GC.Q.slope$a.error, digits = 3, scientific = TRUE)),
                        "b (Q) =",
                        paste(format(GC.Q.slope$b, digits = 3, scientific = TRUE), "\u00b1", format(GC.Q.slope$b.error, digits = 3, scientific = TRUE))
      )
      fitting <- matrix(data = fitting.text,
                        nrow = 2,
                        ncol = 2,
                        byrow=TRUE)
      
      fitting.color <- matrix(data = c(1,1,
                                       1,1),
                              nrow = 2,
                              ncol = 2,
                              byrow=TRUE) 
      
    }else if(length(GC.I.line)>0){
      fitting.text <- c("a (I) =",
                        paste(format(GC.I.slope$a, digits = 3, scientific = TRUE), "\u00b1", format(GC.I.slope$a.error, digits = 3, scientific = TRUE)),
                        "b (I) =",
                        paste(format(GC.I.slope$b, digits = 3, scientific = TRUE), "\u00b1", format(GC.I.slope$b.error, digits = 3, scientific = TRUE))
      )
      fitting <- matrix(data = fitting.text,
                        nrow = 2,
                        ncol = 2,
                        byrow=TRUE)
      
      fitting.color <- matrix(data = c(1,1,
                                       1,1),
                              nrow = 2,
                              ncol = 2,
                              byrow=TRUE) 
    }
  }
  
  textplot(object= fitting, 
           col.data = fitting.color,
           cex = 1.2,
           halign = "center",
           valign="top",
           show.colnames= FALSE,
           show.rownames= FALSE)
  
  title(main = fitting.title)
  
  #Results ------------------------------------------------------------
  results.title <- "Results"
  results.subtitle <-"D\u2091 = Q + I"
  
  results.text <- c("D\u2091 (DP):", 
                    paste(round(De.DP, digits = 2), "\u00b1", round(De.DP.error, digits = 2)),
                    " ",
                    paste( "(", round(De.DP.error/De.DP*100, digits = 2), "%)",sep = ""),
                    "D\u2091 (GC):", 
                    paste(round(De.GC, digits = 2), "\u00b1", round(De.GC.error, digits = 2)),
                    " ",
                    paste( "(", round(De.GC.error/De.GC*100, digits = 2), "%)",sep = ""))
  
  results <- matrix(data = results.text,
                    nrow = 4,
                    ncol = 2,
                    byrow=TRUE)
  
  results.color <- matrix(data = c(2,2,
                                   2,2,
                                   3,3,
                                   3,3),
                          nrow = 4,
                          ncol = 2,
                          byrow=TRUE)  
  
  textplot(object= results, 
           col.data = results.color,
           cex = 1.2,
           halign = "center",
           valign="center",
           show.colnames= FALSE,
           show.rownames= FALSE)
  
  title(main = results.title)
  
  #Page title ---------------------------------------------------------
  page.title <- paste("MAAD: ",
                      sample.name,
                      " - page ",
                      page,
                      ": Palaeodose estimation | fit: ",
                      fit.method,
                      if(fit.weighted){" (weighted)"},
                      sep = "")
  
  mtext(page.title, outer=TRUE,font = 2)
  
  #clean layout...
  layout(1)
  par(old.par)
}