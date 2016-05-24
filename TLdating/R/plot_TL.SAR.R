#' plots MAAD results
#' 
#' This function plots the results obtained by the analyse_TL.MAAD function.
#' 
#' @param sample.name
#'  \link{character} (\bold{required}): Sample name.
#' @param sample.position
#'  \link{integer} (\bold{required}): aliquot position.
#' @param plotting.parameters
#'  \link{list} (with default): list containing the plotting parameters. See details.
#' @param eval.Tmin
#'  \link{integer} (\bold{required}): Temperature (°C) of the lower boundary for the signal integration.
#' @param eval.Tmax
#'  \link{integer} (\bold{required}): Temperature (°C) of the upper boundary for the signal integration.
#' @param fitting.parameters
#'  \link{list} (with default): list containing the fitting parameters. See details.
#' @param temperatures
#'  \link{numeric} (\bold{required}): temperature vector 
#' @param names
#'  \link{character} (\bold{required}): Name vector for the regenerative curves.
#' @param names.duplicated
#'  \link{character} (\bold{required}): Name vector for the duplicated doses.
#' @param doses
#'  \link{numeric} (\bold{required}): Dose vector for the regenerative curves.
#' @param Lx
#'  \link{numeric} (\bold{required}): Lx matrix for the regenerative curves.
#' @param Tx
#'  \link{numeric} (\bold{required}): Tx matrix for the regenerative curves.
#' @param LxTx
#'  \link{numeric} (\bold{required}): Lx/Tx matrix for the regenerative curves.
#' @param Lx.plateau
#'  \link{numeric} (\bold{required}): Ln/Lx matrix for the regenerative curves.
#' @param Tx.plateau
#'  \link{numeric} (\bold{required}): Tn/Tx matrix for the regenerative curves.
#' @param LxTx.plateau
#'  \link{numeric} (\bold{required}): (Ln/Tn)/(Lx/Tx) matrix for the regenerative curves.
#' @param DP.Q.line
#'  \link{numeric} (\bold{required}): Vector containing the estimation of Q for each T° step.
#' @param DP.Q.line.error
#'  \link{numeric} (\bold{required}): Vector containing the uncertainty on the estimation of Q for each T° step.
#' @param Q.DP
#'  \link{numeric} (\bold{required}): Q estimation using the dose plateau approach
#' @param Q.DP.error
#'  \link{numeric} (\bold{required}): Uncertainty on the Q estimation using the dose plateau approach
#' @param GC.Q.line
#'  \link{numeric} (\bold{required}): growth curve for Q
#' @param GC.Q.LxTx
#'  \link{numeric} (\bold{required}): Lx/Tx vector used for Q estimation using the growth curve approach.
#' @param GC.Q.LxTx.error
#'  \link{numeric} (\bold{required}): Error on the Lx/tx vector used for Q estimation using the growth curve approach.
#' @param GC.Q.slope
#'  \link{numeric} (\bold{required}): growth curve parameters for Q
#' @param Q.GC
#'  \link{numeric} (\bold{required}): Q estimation using the growth curve approach
#' @param Q.GC.error
#'  \link{numeric} (\bold{required}): Uncertainty on the Q estimation using the growth curve approach
#' @param TxTn
#'  \link{numeric} (\bold{required}): average Tx/Tn value for the regenerative curves.
#' @param rejection.values
#'  \link{list} (\bold{required}): result of the rejection tests.

#' @details 
#' 
#' \bold{Fitting parameters} \cr
#' The fitting parameters are:  \cr
#' \describe{
#'  \item{\code{method}}{
#'    \link{character}: Fitting method (\code{LIN}, \code{EXP}, \code{EXP+LIN} or \code{EXP+EXP}).}
#'  \item{\code{fit.weighted}}{
#'    \link{logical}: If the fitting is weighted or not.}
#'  \item{\code{fit.rDoses.min}}{
#'    \link{numeric}: Lower regenerative dose used for the fitting.}
#'  \item{\code{fit.rDoses.max}}{
#'    \link{numeric}: Higher regenerative dose used for the fitting.}
#' }
#' See also \link{calc_TL.SAR.fit}. \cr
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
#' See also \link{plot_TL.SAR}. \cr
#' 
#' @author David Strebler
#'  
## @export plot_TL.SAR
 
plot_TL.SAR <- function(
  
  sample.name,
  sample.position,
  
  fitting.parameters=list(fit.method="LIN",
                          fit.weighted=FALSE,
                          fit.rDoses.min=NA,
                          fit.rDoses.max=NA),
  
  eval.Tmin,
  eval.Tmax,

  temperatures,
  
  names,
  names.duplicated,
  doses,

  Lx,
  Tx,
  LxTx,
  Lx.plateau,
  Tx.plateau,
  LxTx.plateau,

  DP.Q.line,
  DP.Q.line.error,
  GC.Q.line,
  GC.Q.LxTx,
  GC.Q.LxTx.error,
  GC.Q.slope,
  
  Q.DP,
  Q.DP.error,
  Q.GC,
  Q.GC.error,

  TxTn,
  
  rejection.values,
  
  plotting.parameters=list(plot.Tmin=0,
                           plot.Tmax=NA)
  
  ){ 
  # ------------------------------------------------------------------------------  
  # Integrity Check
  # ------------------------------------------------------------------------------  
  if (missing(sample.name)){
    stop("[plot_TL.SAR] Error: Input 'sample.name' is missing.") 
  }else if (!is.character(sample.name)){
    stop("[plot_TL.SAR] Error: Input 'sample.name' is not of type 'character'.")
  } 
  
  # ...
  
  # ------------------------------------------------------------------------------  
  
  Tmax <- max(temperatures)
  nPoints <- length(temperatures)
  Tstep <- Tmax/nPoints
  eval.min <- ceiling(eval.Tmin/Tstep)
  eval.max <-floor(eval.Tmax/Tstep)
  
  
  recycling.ratio <- rejection.values$recycling.ratio
  recuperation.rate <- rejection.values$recuperation.rate
  Lx.error.max <- rejection.values$Lx.error.max
  Tx.error.max <- rejection.values$Tx.error.max
  test.recycling <- rejection.values$test.recycling
  test.recuperation <- rejection.values$test.recuperation
  test.Lx.error <- rejection.values$test.Lx.error
  test.Tx.error <- rejection.values$test.Tx.error
  
  plot.Tmin <- plotting.parameters$plot.Tmin
  plot.Tmax <- plotting.parameters$plot.Tmax
  
  fit.method <- fitting.parameters$fit.method
  fit.weighted <- fitting.parameters$fit.weighted
  
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
  
  #----------------------------------------------------------------------------------------------------------------
  #Plot results
  #----------------------------------------------------------------------------------------------------------------
  
  
  old.par <- par( no.readonly = TRUE )
  par( oma = c( 0, 0, 3, 0 ) )
  
  # --------------------------------------------------------------------
  #page 1
  # --------------------------------------------------------------------
  
  #Layout
  layout(matrix(c(1,2,3,4,5,6,7,7), 4, 2, byrow = TRUE),heights = c(2,2,2,1))
  
  #color
  ref_colors <- rainbow(n=length(names)-1)
  colors <- seq(length(names))
  names(colors) <- names
  colors[names(colors)=="N"] <- 1
  colors[names(colors)!="N"] <- ref_colors
  
  #Lx
  plot.Lx.max <- max(Lx)
  
  rLx <- Lx[,colnames(Lx) != "N"]
  nLx <- Lx[,colnames(Lx) == "N"]
  
  # rLx
  for(i in 1 : ncol(rLx)){
    temp.curve <- rLx[,i]
    temp.name <- colnames(rLx)[i]
    temp.color <- colors[temp.name]
    
    if(i == 1) {
      plot(x=temperatures, 
           y=temp.curve, 
           type="l", 
           col=temp.color, 
           xlim=c(plot.Tmin,plot.Tmax),
           ylim=c(0,plot.Lx.max), 
           xlab = "Temperature (\u00b0C)", 
           ylab = "Lx"
      )
      par(new = TRUE)
      
    }else{
      lines(x=temperatures, 
            y=temp.curve, 
            col=temp.color,
            xlim=c(plot.Tmin,plot.Tmax),
            ylim=c(0,plot.Lx.max))
    }
  }
  
  # nLx
  lines(x=temperatures, 
        y=nLx, 
        col=colors["N"],
        lwd=2,
        xlim=c(plot.Tmin,plot.Tmax),
        ylim=c(0,plot.Lx.max))
  
  # T boundaries
  abline(v=eval.Tmin,col=3,lty=3)
  abline(v=eval.Tmax,col=2,lty=3)
  
  par(new = FALSE)
  
  #Lx.plateau
  plot.Lx.plateau.max <- max(Lx.plateau[eval.min:eval.max,])*1.5
  
  for(i in 1 : ncol(Lx.plateau)){  
    temp.curve <- Lx.plateau[,i]
    temp.name <- colnames(Lx.plateau)[i]
    temp.color <- colors[temp.name]
    
    if(i == 1) {
      plot(x=temperatures, 
           y=temp.curve, 
           type="l", 
           col=temp.color, 
           xlim=c(plot.Tmin,plot.Tmax),
           ylim=c(0,plot.Lx.plateau.max), 
           xlab = "Temperature (\u00b0C)", 
           ylab = "Lx Plateau"
      )
      par(new = TRUE)
      
    }else{
      lines(x=temperatures, 
            y=temp.curve, 
            col=temp.color, 
            xlim=c(plot.Tmin,plot.Tmax),
            ylim=c(0,plot.Lx.plateau.max) 
      )
    }
  }  
  abline(v=eval.Tmin,col=3,lty=3)
  abline(v=eval.Tmax,col=2,lty=3)
  
  par(new = FALSE)
  
  #Tx
  plot.Tx.max <- max(Tx)
  
  rTx <- Tx[,colnames(Tx) != "N"]
  nTx <- Tx[,colnames(Tx) == "N"]
  
  # rTx
  for(i in 1 : ncol(rTx)){
    temp.curve <- rTx[,i]
    temp.name <- colnames(rTx)[i]
    temp.color <- colors[temp.name]
    
    if(i == 1) {
      plot(x=temperatures, 
           y=temp.curve, 
           type="l", 
           col=temp.color, 
           xlim=c(plot.Tmin,plot.Tmax),
           ylim=c(0,plot.Tx.max), 
           xlab = "Temperature (\u00b0C)", 
           ylab = "Tx"
      )
      par(new = TRUE)
      
    }else{
      lines(x=temperatures, 
            y=temp.curve, 
            col=temp.color, 
            xlim=c(plot.Tmin,plot.Tmax),
            ylim=c(0,plot.Tx.max) 
      )
    }
  }  
  
  # nTx
  lines(x=temperatures, 
        y=nTx, 
        col=colors["N"],
        lwd=2,
        xlim=c(plot.Tmin,plot.Tmax),
        ylim=c(0,plot.Tx.max))
  
  # T boundaries
  abline(v=eval.Tmin,col=3,lty=3)
  abline(v=eval.Tmax,col=2,lty=3)
  
  par(new = FALSE)
  
  #Tx.plateau
  plot.Tx.plateau.max <- max(Tx.plateau[eval.min:eval.max,])*1.5
  
  for(i in 1 : ncol(Tx.plateau)){    
    temp.curve <- Tx.plateau[,i]
    temp.name <- colnames(Tx.plateau)[i]
    temp.color <- colors[temp.name]
    
    if(i == 1) {
      plot(x=temperatures, 
           y=temp.curve, 
           type="l", 
           col=temp.color, 
           xlim=c(plot.Tmin,plot.Tmax),
           ylim=c(0,plot.Tx.plateau.max), 
           xlab = "Temperature (\u00b0C)", 
           ylab = "Tx Plateau"
      )
      par(new = TRUE)
      
    }else{
      lines(x=temperatures, 
            y=temp.curve, 
            col=temp.color, 
            xlim=c(plot.Tmin,plot.Tmax),
            ylim=c(0,plot.Tx.plateau.max) 
      )
    }
  }  
  
  abline(v=eval.Tmin,col=3,lty=3)
  abline(v=eval.Tmax,col=2,lty=3)
  
  par(new = FALSE)
  
  #LxTx
  plot.LxTx.max <- max(LxTx[eval.min:eval.max,])
  
  rLxTx <- LxTx[,colnames(LxTx) != "N"]
  nLxTx <- LxTx[,colnames(LxTx) == "N"]
  
  # rLxTx
  for(i in 1 : ncol(rLxTx)){    
    temp.curve <- rLxTx[,i]
    temp.name <- colnames(rLxTx)[i]
    temp.color <- colors[temp.name]
    
    if(i == 1) {
      plot(x=temperatures, 
           y=temp.curve, 
           type="l", 
           col=temp.color, 
           xlim=c(plot.Tmin,plot.Tmax),
           ylim=c(0,plot.LxTx.max), 
           xlab = "Temperature (\u00b0C)", 
           ylab = "Lx/Tx"
      )
      par(new = TRUE)
      
    }else{
      lines(x=temperatures, 
            y=temp.curve, 
            col=temp.color, 
            xlim=c(plot.Tmin,plot.Tmax),
            ylim=c(0,plot.LxTx.max) 
      )
    }
  }  
  
  #nLxTx
  lines(x=temperatures, 
        y=nLxTx, 
        col=colors["N"],
        lwd=2,
        xlim=c(plot.Tmin,plot.Tmax),
        ylim=c(0,plot.LxTx.max))
  
  #T boundaries
  abline(v=eval.Tmin,col=3,lty=3)
  abline(v=eval.Tmax,col=2,lty=3)
  
  par(new = FALSE)
  
  #LxTx.plateau
  plot.LxTx.plateau.max <- max(LxTx.plateau[eval.min:eval.max,])*1.5
  
  for(i in 1 : ncol(LxTx.plateau)){      
    temp.curve <- LxTx.plateau[,i]
    temp.name <- colnames(LxTx.plateau)[i]
    temp.color <- colors[temp.name]
    
    if(i == 1) {
      plot(x=temperatures, 
           y=temp.curve, 
           type="l", 
           col=temp.color, 
           xlim=c(plot.Tmin, plot.Tmax),
           ylim=c(0, plot.LxTx.plateau.max), 
           xlab = "Temperature (\u00b0C)", 
           ylab = "Lx/Tx Plateau"
      )
      par(new = TRUE)
      
    }else{
      lines(x=temperatures, 
            y=temp.curve, 
            col=temp.color, 
            xlim=c(plot.Tmin,plot.Tmax),
            ylim=c(0,plot.LxTx.plateau.max) 
      )
    }
  }  
  abline(v=eval.Tmin,col=3,lty=3)
  abline(v=eval.Tmax,col=2,lty=3)
  
  par(new = FALSE)
  
  #Legend
  legend.title <- "Legend"
  legend.size <- length(doses)
  legend.text <- c(names, doses)
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
  title(legend.title)
  
  #Page title
  page.title <- paste("SAR: ",
                      sample.name,
                      ", Disk ", 
                      sample.position, 
                      " - page 1: Regenerative doses",
                      sep = "")
  mtext(page.title, outer=TRUE,font = 2)
  
  # --------------------------------------------------------------------
  #page 2
  # --------------------------------------------------------------------
  
  #Layout
  layout(matrix(c(1,1,3,1,1,4,2,2,5,2,2,6), 4, 3, byrow = TRUE))
  
  # Plotting  Palaeodose (DP.Q.line) ----------------------------------------
  
  plot.DP.Q.line.max <- max(DP.Q.line[eval.min:eval.max])*1.5
  
  par(mar=c(5,4,4,1))
  plot(x=temperatures, 
       y=DP.Q.line,
       xlim=c(plot.Tmin, plot.Tmax),
       ylim=c(0, plot.DP.Q.line.max), 
       main = "Palaeodose (DP)",
       sub = paste("D\u2091 =", 
                   round(Q.DP, digits = 2), "\u00b1", round(Q.DP.error, digits = 2),
                   paste( "(", round(Q.DP.error/Q.DP*100, digits = 2), "%)",sep = "")),
       xlab = "Temperature (\u00b0C)", 
       ylab = "Dose (s)",
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
  
  # Plotting  Growth curve (GC) ----------------------------------------
  doses.nat <- doses[names=="N"]
  LxTx.nat <- GC.Q.LxTx[names=="N"]
  LxTx.nat.error <- GC.Q.LxTx.error[names=="N"]
  
  doses.0 <- doses[names!="N" & doses==0]
  LxTx.0 <- GC.Q.LxTx[names!="N" & doses==0]
  LxTx.0.error <- GC.Q.LxTx.error[names!="N" & doses==0]
  
  doses.dupl <- doses[names %in% names.duplicated]
  LxTx.dupl <- GC.Q.LxTx[names %in% names.duplicated]
  LxTx.dupl.error <- GC.Q.LxTx.error[names %in% names.duplicated]
  
  doses.r <- doses[doses!=0 & !(names %in% names.duplicated)]
  LxTx.r <- GC.Q.LxTx[doses!=0 & !(names %in% names.duplicated)]
  LxTx.r.error <- GC.Q.LxTx.error[doses!=0 & !(names %in% names.duplicated)]
  
  #LxTx
  par(mar=c(5,4,4,1))
  plot(x=doses.r, 
       y=LxTx.r, 
       xlim=c(0, max(doses)),
       ylim=c(0, max(LxTx.r)), 
       main = "Growth curve",
       sub = paste("D\u2091 (GC) =", 
                   round(Q.GC, digits = 2), "\u00b1", round(Q.GC.error, digits = 2),
                   paste( "(", round(Q.GC.error/Q.GC*100, digits = 2), "%)",sep = "")),
       xlab = "Doses", 
       ylab = "Lx/Tx",
       type="p", 
       pch=18,
       col=1)
  par(new = TRUE)
  
  # error
  arrows(doses.r, 
         LxTx.r-LxTx.r.error, 
         doses.r, 
         LxTx.r+LxTx.r.error, 
         length=0, #0.05,
         angle=90, 
         code=3)
  
  # R0
  points(x=doses.0,
         y=LxTx.0,
         pch=5,
         col=1
  )
  # error
  arrows(doses.0, 
         LxTx.0-LxTx.0.error, 
         doses.0, 
         LxTx.0+LxTx.0.error, 
         length=0, #0.05,
         angle=90, 
         code=3)
  
  
  # Duplicated
  points(x=doses.dupl,
         y=LxTx.dupl,
         pch=2)
  # error
  arrows(doses.dupl, 
         LxTx.dupl-LxTx.dupl.error, 
         doses.dupl, 
         LxTx.dupl+LxTx.dupl.error, 
         length=0, #0.05,
         angle=90, 
         code=3)
  
  # Natural
  points(x=0,
         y=LxTx.nat,
         pch=18,
         col=3)
  # error
  arrows(x0 = 0, 
         y0 = LxTx.nat-LxTx.nat.error, 
         x1 = 0, 
         y1 = LxTx.nat+LxTx.nat.error, 
         length = 0, #0.05, 
         angle = 90, 
         code=3)
  
  # Q.GC
  points(x=Q.GC,
         y=0,
         pch=18,
         col=3) 
  # error
  arrows(Q.GC-Q.GC.error, 
         0, 
         Q.GC+Q.GC.error, 
         0, 
         length=0, #0.05,
         angle=90, 
         code=3)
  
  # Q
  points(x=Q.DP,
         y=0,
         pch=1,
         col=2) 
  
  #Growth curve
    abline(GC.Q.line)    
  
  # crossing
  points(x=Q.GC,
         y=LxTx.nat,
         pch=18,
         col=3)
  
  segments(x0=0,
           y0=LxTx.nat,
           x1=Q.GC,
           y1=LxTx.nat,
           col=3)
  
  segments(x0=Q.GC,
           y0=LxTx.nat,
           x1=Q.GC,
           y1=0,
           col=3)
  
  # Legend
  legend(x = "topleft", 
         legend = c("REG points", "REG point repeated", "REG point 0", "D\u2091 (GC)", "D\u2091"),
         pch = c(18,2,5,18,5),
         col = c(1,1,1,3,2),
         bty = "n")
  
  par(new = FALSE)
  
  # Plotting  Testdose response ----------------------------------------
  
  temp.x <- seq(from = 1,to = length(TxTn),by = 1)
  
  par(mar=c(5,4,4,2))
  plot(x=temp.x, 
       y=TxTn, 
       type="o", 
       col=1, 
       xlim=c(0, max(temp.x)),
       ylim=c(0, max(TxTn)), 
       main = "Testdose response",
       xlab = "SAR cycle", 
       ylab = "Tx/Tn"
  )
  par(new = TRUE)
  abline(h=1,col=2,lty=3)
  par(new = FALSE)
  
  # Plotting  rejection criteria ----------------------------------------
  rejection.title <- "Rejection criteria"
  
  #rejection test
  rejection.text <- vector()
  for(i in 1: length(recycling.ratio)){
    if(i==1){
      temp.text <- "Recycling ratio:"      
    }else{
      temp.text <- ""
    }
    
    temp.ratio <- format(recycling.ratio[i],digits = 2)
    
    temp.name <-  paste("(", names(recycling.ratio[i]), ")", sep="")
    
    rejection.text <- c(rejection.text, temp.text, temp.ratio, temp.name)
  }
  
  for(i in 1: length(recuperation.rate)){
    if(i==1){
      temp.text <- "Recuperation rate:"
    }else{
      temp.text <- ""
    }
    
    temp.rate <- paste(format(recuperation.rate[i]*100,digits = 2),"%",sep="")
    
    temp.name <-  paste("(", names(recuperation.rate[i]), ")", sep="")
    
    rejection.text <- c(rejection.text, temp.text, temp.rate, temp.name)
  }
  
  rejection.text <- c(rejection.text, 
                      "Lx error (max):",
                      paste(format(Lx.error.max*100,digits = 2),"%",sep = ""),
                      "")
  
  rejection.text <- c(rejection.text,                    
                      "Tx error (max):",
                      paste(format(Tx.error.max*100,digits = 2),"%",sep = ""),
                      "")
  
  #rejection color
  data.color <- vector()
  
  for(i in 1: length(test.recycling)){
    if(test.recycling[i]){
      temp.color <- 1
    }else{
      temp.color <- 2
    }
    data.color <- c(data.color,temp.color,temp.color,temp.color)
  }
  
  for(i in 1: length(test.recuperation)){
    if(test.recuperation[i]){
      temp.color <- 1
    }else{
      temp.color <- 2
    }
    data.color <- c(data.color,temp.color,temp.color,temp.color)
  }
  
  if(test.Lx.error){
    temp.color <- 1
  }else{
    temp.color <- 2
  }
  data.color <- c(data.color,temp.color,temp.color,temp.color)
  
  if(test.Tx.error){
    temp.color <- 1
  }else{
    temp.color <- 2
  }
  data.color <- c(data.color,temp.color,temp.color,temp.color)
  
  
  rejection.color <- matrix(data = data.color,
                            nrow = length(data.color)/3,
                            ncol = 3,
                            byrow=TRUE)
  
  rejection <- matrix(data = rejection.text, 
                      nrow = length(rejection.text)/3,
                      ncol = 3,
                      byrow = TRUE
  )
  
  textplot(object=rejection, 
           col.data=rejection.color,
           cex=1.2,
           halign="center",
           valign="top",
           show.colnames= FALSE,
           show.rownames= FALSE
  )
  title(rejection.title)
  
  # Curve fitting -----------------------------------------------------
  
  if(fit.method == "LIN"){
    fitting.title <- paste("Curve fitting (GC):", "\n",                           
                           "Linear", if(fit.weighted){"(weighted)"}, "\n",
                           "y = a + bx")
    
    fitting.text <- c("a =",
                      paste(format(GC.Q.slope$a, digits = 3, scientific = TRUE), 
                            "\u00b1", 
                            format(GC.Q.slope$a.error, digits = 3, scientific = TRUE)),
                      "b =",
                      paste(format(GC.Q.slope$b, digits = 3, scientific = TRUE), 
                            "\u00b1", 
                            format(GC.Q.slope$b.error, digits = 3, scientific = TRUE)))
    
    fitting <- matrix(data = fitting.text,
                      nrow = 2,
                      ncol = 2,
                      byrow=TRUE)
    
    fitting.color <- matrix(data = c(1,1,
                                     1,1),
                            nrow = 2,
                            ncol = 2,
                            byrow=TRUE)
    
  }else if(fit.method == "EXP"){
    fitting.title <- paste("Curve fitting (GC):", "\n",                           
                           "Exponential", if(fit.weighted){"(weighted)"}, "\n",
                           "y = a.(1-exp(-(x+c)/b))")
    
    fitting.text <- c("a =",
                      paste(format(GC.Q.slope$a, digits = 3, scientific = TRUE), 
                            "\u00b1", 
                            format(GC.Q.slope$a.error, digits = 3, scientific = TRUE)),
                      "b =",
                      paste(format(GC.Q.slope$b, digits = 3, scientific = TRUE), 
                            "\u00b1", 
                            format(GC.Q.slope$b.error, digits = 3, scientific = TRUE)),
                      "c =",
                      paste(format(GC.Q.slope$c, digits = 3, scientific = TRUE), 
                            "\u00b1", 
                            format(GC.Q.slope$c.error, digits = 3, scientific = TRUE)))
    
    
    fitting <- matrix(data = fitting.text,
                      nrow = 3,
                      ncol = 2,
                      byrow=TRUE)
    
    fitting.color <- matrix(data = c(1,1,
                                     1,1,
                                     1,1),
                            nrow = 3,
                            ncol = 2,
                            byrow=TRUE)
  }else{
    fitting.title <- paste("Curve fitting (GC):", "\n",                           
                           "Unknown")
    
    textplot("")
    
    title(main = fitting.title)
  }
  
  textplot(object= fitting, 
           col.data = fitting.color,
           cex = 1.2,
           halign = "center",
           valign="top",
           show.colnames= FALSE,
           show.rownames= FALSE)
  
  title(main = fitting.title)
  
  # Results ----------------------------------------
  results.text <- c("D\u2091 (DP):",
                    paste(round(Q.DP, digits = 2),
                          "\u00b1", 
                          round(Q.DP.error, digits = 2)),
                    " ",
                    paste("(\u00b1", 
                          round(Q.DP.error/Q.DP*100, digits = 2), 
                          "%)",
                          sep = ""),
                    "D\u2091 (GC):",
                    paste(round(Q.GC, digits = 2),
                          "\u00b1", 
                          round(Q.GC.error, digits = 2)),
                    " ",
                    paste("(\u00b1", 
                          round(Q.GC.error/Q.GC*100, digits = 2), 
                          "%)",
                          sep = ""))
  
  results.color <- matrix(data = c(2,2,
                                   2,2,
                                   3,3,
                                   3,3),
                          nrow = 4,
                          ncol = 2,
                          byrow=TRUE)
  
  results <- matrix(data = results.text, 
                    nrow = 4,
                    ncol = 2,
                    byrow = TRUE)
  
  #par(mar=c(1,1,4,1))
  
  textplot(object=results, 
           col.data=results.color,
           cex=1.2,
           halign="center",
           valign="top",
           show.colnames= FALSE,
           show.rownames= FALSE)
  
  title("Results")
  
  #Page title ----------------------------------------
  page.title <- paste("SAR: ",
                      sample.name,
                      ", disk ", 
                      sample.position, 
                      " - page 2: Paleodose estimation | fit:",
                      fit.method,
                      if(fit.weighted){" (weighted)"},
                      sep = "")
  mtext(page.title, outer=TRUE,font = 2)
  
  # --------------------------------------------------------------------
  # Clean layout
  layout(1)
  par(old.par)  
}