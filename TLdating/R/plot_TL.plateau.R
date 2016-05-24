#' plot plateau test result
#' 
#' This function plots the results for \link{analyse_TL.plateau}.
#' 
#' @param sample.name
#'  \link{character} (\bold{required}): Sample name.
#' @param temperatures
#'  \link{numeric} (\bold{required}): temperature vector 
#' @param names
#'  \link{character} (\bold{required}): Name vector for the additive curves.
#' @param doses
#'  \link{numeric} (\bold{required}): Dose vector for the additive curves.
#' @param Lx
#'  \link{numeric} (\bold{required}): Lx matrix for the additive curves.
#' @param Lx.a
#'  \link{numeric} (\bold{required}): Lx matrix for the average additive curves.
#' @param Lx.plateau
#'  \link{numeric} (\bold{required}): Ln/Lx matrix for the additive curves.
#' @param LxTx
#'  \link{numeric} (\bold{required}): Lx/Tx matrix for the additive curves.
#' @param LxTx.a
#'  \link{numeric} (\bold{required}): Lx/Tx matrix for the average additive curves.
#' @param LxTx.plateau
#'  \link{numeric} (\bold{required}): (Ln/Tn)/(Lx/Tx) matrix for the additive curves.
#' @param plotting.parameters
#'  \link{list} (with default): list containing the plotting parameters. See details.
#'  
#' @details 
#' 
#' \bold{Plotting parameters} \cr
#' The plotting parameters are:  \cr
#' \describe{
#'  \item{\code{plot.Tmin}}{
#'    \link{numeric}: Lowest temperature plotted.} 
#'  \item{\code{plot.Tmax}}{
#'    \link{numeric}: Highest temperature plotted.}
#'  \item{\code{no.plot}}{
#'    \link{logical}: If \code{TRUE}, the results will not be plotted.}
#' }
#' See also \link{analyse_TL.MAAD}. \cr
#'
#' @seealso 
#'  \link{analyse_TL.plateau},
#'  \link{calc_TL.MAAD.fit.Q},
#'  \link{calc_TL.MAAD.fit.I}.
#'  
#' @author David Strebler
#' 
## @export plot_TL.plateau

plot_TL.plateau <- function(
  
  sample.name,
  temperatures,
  names,
  doses,
  Lx,
  Lx.a,
  Lx.plateau,
  LxTx,
  LxTx.a,
  LxTx.plateau,
  plotting.parameters=list(plateau.Tmin=0,
                           plateau.Tmax=NA,
                           plot.Tmin=0,
                           plot.Tmax=NA)
  
){ 
  # ------------------------------------------------------------------------------  
  # Integrity Check
  # ------------------------------------------------------------------------------  
  if (missing(sample.name)){
    stop("[plot_TL.plateau] Error: Input 'sample.name' is missing.") 
  }else if (!is.character(sample.name)){
    stop("[plot_TL.plateau] Error: Input 'sample.name' is not of type 'character'.")
  } 
  
  # ...
  
  # ------------------------------------------------------------------------------  
  
  Tmax <- max(temperatures)
  nPoints <- length(temperatures)
  Tstep <- Tmax/nPoints
  
  plateau.Tmin <- plotting.parameters$plateau.Tmin
  plateau.Tmax <- plotting.parameters$plateau.Tmax
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
    if(!is.finite(plot.Tmax) || is.null(plot.Tmax)){
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
  
  if(!is.numeric(plateau.Tmin)){
    if(is.na(plateau.Tmin)  || is.null(plateau.Tmin)){
      plateau.Tmin <- plot.Tmin
    }else{
      stop("[plot_TL.MAAD] Error: plateau.Tmin is not numeric.")    
    }
  }
  
  if(!is.numeric(plateau.Tmax)){
    if(is.na(plateau.Tmax) || is.null(plateau.Tmax)){
      plateau.Tmax <- plot.Tmax
    }else{
      stop("[plot_TL.MAAD] Error: plateau.Tmax is not numeric.")      
    }
  }
  
  if(plateau.Tmin > plateau.Tmax){
    stop("[plot_TL.MAAD] Error: plateau.Tmin > plateau.Tmax")
  }
  
  if(plateau.Tmin < plot.Tmin){
    plateau.Tmin <- plot.Tmin
  }
  
  if(plateau.Tmax > plot.Tmax){
    plateau.Tmax <- plot.Tmax
  }
  # -------------------------------
  plateau.min <- ceiling(plateau.Tmin/Tstep)
  plateau.max <- ceiling(plateau.Tmax/Tstep)
  
  uDoses <- unique(doses)
  uNames <- unique(names)
  #----------------------------------------------------------------------------------------------------------------
  #Plot results
  #----------------------------------------------------------------------------------------------------------------
  old.par <- par( no.readonly = TRUE )
  par( oma = c(0.5, 0, 3, 0 ) )
  
  #---------------------------------------------------------------------------
  #Page 1
  #---------------------------------------------------------------------------
  if(length(Lx) > 0){
    #Layout
    layout(matrix(c(1,2,1,2,3,3), 3, 2, byrow = TRUE))
    
    #color
    ref_colors <- rainbow(n=length(uNames)-1)
    colors <- seq(length(uNames))
    names(colors) <- uNames
    colors[names(colors)=="N"] <- 1
    colors[names(colors)!="N"] <- ref_colors
    
    
    #Lx (additive)
    plot.Lx.max <- max(Lx[plateau.min:plateau.max,],na.rm = TRUE)*1.1
    
    for(i in 1 : length(uDoses)){
      temp.curve <- Lx.a[,i]
      temp.name <- colnames(Lx.a)[i]
      temp.color <- colors[temp.name]
      
      if(i == 1) {
        plot(x=temperatures, 
             y=temp.curve, 
             type="l", 
             col=temp.color, 
             xlim=c(plot.Tmin,plot.Tmax),
             ylim=c(0,plot.Lx.max), 
             main = "Lx",
             xlab = "Temperature (\u00b0C)", 
             ylab = "Luminescence signal"
        )
        par(new = TRUE)
        
      }else{
        lines(x=temperatures, 
              y=temp.curve, 
              col=temp.color,
              xlim=c(plot.Tmin,plot.Tmax),
              ylim=c(0,plot.Lx.max) 
        )
      }
    }
    
    for(i in 1: length(doses)){
      temp.curve <- Lx[,i]
      temp.name <- colnames(Lx)[i]
      temp.color <- colors[temp.name]
      
      lines(x=temperatures, type = "p", 
            y=temp.curve, 
            col=temp.color,
            xlim=c(plot.Tmin,plot.Tmax),
            ylim=c(0,plot.Lx.max),
            pch=18) 
    }
    
    par(new = FALSE)
    
    #Lx.plateau 
    plot.Lx.plateau.max <- max(Lx.plateau[plateau.min:plateau.max,],na.rm = TRUE)*1.2
    
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
              ylim=c(0,plot.Lx.plateau.max) 
        )
      }
    }
    par(new = FALSE)
    
    #Legend
    legend.size <- length(uDoses)
    legend.text <- c(uNames, uDoses)
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
    page.title <- paste("Plateau test:",
                        sample.name,
                        " - page 1",
                        sep = " ")
    mtext(page.title, outer=TRUE,font = 2)
  }
  
  #---------------------------------------------------------------------------
  #Page 2
  #---------------------------------------------------------------------------
  if(length(LxTx) > 0){
    #Layout
    layout(matrix(c(1,2,1,2,3,3), 3, 2, byrow = TRUE))
    
    #color
    ref_colors <- rainbow(n=length(uNames)-1)
    colors <- seq(length(uNames))
    names(colors) <- uNames
    colors[names(colors)=="N"] <- 1
    colors[names(colors)!="N"] <- ref_colors
    
    #LxTx (additive)
    plot.LxTx.max <- max(LxTx[plateau.min:plateau.max,],na.rm = TRUE)*1.1
    
    for(i in 1 : length(uDoses)){
      temp.curve <- LxTx.a[,i]
      temp.name <- colnames(LxTx.a)[i]
      temp.color <- colors[temp.name]
      
      if(i == 1) {
        plot(x=temperatures, 
             y=temp.curve, 
             type="l", 
             col=temp.color, 
             xlim=c(plot.Tmin,plot.Tmax),
             ylim=c(0,plot.LxTx.max), 
             main = "LxTx",
             xlab = "Temperature (\u00b0C)", 
             ylab = "Luminescence signal"
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
    for(i in 1: length(doses)){
      temp.curve <- LxTx[,i]
      temp.name <- colnames(LxTx)[i]
      temp.color <- colors[temp.name]
      
      lines(x=temperatures, type = "p", 
            y=temp.curve, 
            col=temp.color,
            xlim=c(plot.Tmin,plot.Tmax),
            ylim=c(0,plot.LxTx.max),
            pch=18) 
    }
    
    par(new = FALSE)
    
    #LxTx.plateau 
    plot.LxTx.plateau.max <- max(LxTx.plateau[plateau.min:plateau.max,],na.rm = TRUE)*1.2
    
    for(i in 1 : ncol(LxTx.plateau)){  
      temp.curve <- LxTx.plateau[,i]
      temp.name <- colnames(LxTx.plateau)[i]
      temp.color <- colors[temp.name]
      
      if(i == 1) {
        plot(x=temperatures, 
             y=temp.curve, 
             type="l", 
             col=temp.color, 
             xlim=c(plot.Tmin,plot.Tmax),
             ylim=c(0,plot.LxTx.plateau.max), 
             main = "Plateau test (LxTx)",
             xlab = "Temperature (\u00b0C)", 
             ylab = "Luminescence signal"
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
    par(new = FALSE)
    
    #Legend
    legend.size <- length(uDoses)
    legend.text <- c(uNames, uDoses)
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
    page.title <- paste("Plateau test:",
                        sample.name,
                        " - page 2",
                        sep = " ")
    mtext(page.title, outer=TRUE,font = 2)
  }
  
  #clean layout...
  layout(1)
  par(old.par)
}