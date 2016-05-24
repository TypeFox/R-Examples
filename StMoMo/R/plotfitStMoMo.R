

#' Plot fitted parameters from a stochastic mortality model
#' 
#' Plot fitted parameters of a stochastic mortality model of class 
#' \code{"fitStMoMo"}.
#' 
#' @usage 
#' \method{plot}{fitStMoMo}(x, nCol = 2, parametricbx = TRUE, type = "l", ...)
#' 
#' @param x an object of class \code{"fitStMoMo"} with the fitted 
#' parameters of a stochastic mortality model.
#' @param nCol number of columns to use in the plot.
#' @param parametricbx if \code{FALSE} parametric age-modulating terms, 
#' which don't need to be estimated, are not plotted.
#' @param ... additional arguments to control graphical appearance.
#' See \code{\link[graphics]{plot}}.
#' @param type what type of plot should be drawn. See 
#' \code{\link[graphics]{plot}}.
#' 
#' @examples
#' 
#' #Fit and plot a Lee-Carter model
#' LCfit <- fit(lc(), Dxt = EWMaleData$Dxt, Ext = EWMaleData$Ext, 
#'              ages = EWMaleData$ages, years = EWMaleData$years,
#'              ages.fit = 55:89)
#' plot(LCfit)
#' plot(LCfit, type = "p", pch = 19)
#' 
#' #Fit and plot a CBD model
#' CBDfit <- fit(cbd(), Dxt = EWMaleData$Dxt, Ext = EWMaleData$Ext, 
#'               ages = EWMaleData$ages, years = EWMaleData$years,
#'               ages.fit = 55:89)
#' plot(CBDfit)
#' plot(CBDfit, parametricbx = FALSE)
#' plot(CBDfit, nCol = 1, parametricbx = FALSE, lwd = 2)
#' 
#' @export 
#' @method plot fitStMoMo
plot.fitStMoMo <- function(x, nCol = 2, parametricbx = TRUE, type = "l", ...) {
  
  years <- x$years
  ages <- x$ages
  cohorts <- x$cohorts
  ax <- x$ax
  bx <- x$bx
  kt <- x$kt
  b0x <- x$b0x
  gc <- x$gc
  N <- x$model$N
    
  #Calculate number of plots and rows
  nPlots <- 2 * N + (!is.null(ax)) + 2 * (!is.null(gc))
  is.nonparametric <- function(bx) {
    is.character(bx) && bx == "NP"
  }  
  if (parametricbx == FALSE) {  #substract the parametric plots        
    nParametricbx <- ifelse( is.null(x$model$periodAgeFun), 0, 
                             sum(sapply(x$model$periodAgeFun, 
                                        FUN = function(x) !is.nonparametric(x)))) +
                      (!is.null(x$model$cohortAgeFun) && !is.nonparametric(x$model$cohortAgeFun))
                     
    nPlots <- nPlots - nParametricbx
  }   
  if (!is.null(ax) && nCol == 2 && parametricbx == TRUE) {
    nPlots <- nPlots + 1  # add and empty plot to the rigth of ax
  }  
  
  nRow <- ceiling(nPlots / nCol)  
  oldpar <- par(no.readonly = TRUE)
  par(mfrow = c(nRow, nCol))
  
  #ax
  if (!is.null(ax)) {
    plot(x = ages, y = ax, ylab = "", xlab = "age", 
         main = expression(paste(alpha[x], " vs. x", "")), type = type, ...)
  }
  
  if (!is.null(ax) && nCol == 2 && parametricbx == TRUE) {
    frame()  # add and empty plot to the rigth of ax
  }
      
  # bx, kt
  if (N > 0) {
    for (i in 1:N) {      
      #bx
      if (parametricbx == TRUE || is.nonparametric(x$model$periodAgeFun[[i]])) {
        plot(x = ages, y = bx[, i], ylab = "", xlab = "age", 
             main = substitute(paste(beta[x]^{(i)}, " vs. x", ""), 
                               list(i = i)), type = type, ...)
      }
      #kt
      plot(x = years,y = kt[i, ], ylab = "", xlab = "year", 
           main = substitute(paste(kappa[t]^{(i)}, " vs. t", ""), 
                             list(i = i)), type = type, ...) 
    }
  }
  
  
  if (!is.null(gc) == TRUE) {
    #bx0
    if (parametricbx == TRUE || is.nonparametric(x$model$cohortAgeFun)) {
      plot(x = ages, y = b0x, ylab = "", xlab = "age", 
           main = substitute(paste(beta[x]^{(i)}, " vs. x", ""), 
                             list(i = 0)), type = type, ...)    
    }    
    #gc  
    plot(x = cohorts, y = gc, ylab = "", xlab = "cohort", 
         main = expression(paste(gamma[t-x], " vs. t-x","")),
         type = type, ...)
  }
  par(oldpar)
}
