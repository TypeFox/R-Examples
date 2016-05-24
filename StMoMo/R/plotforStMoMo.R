

#' Plot a forecast from a Stochastic Mortality Model
#' 
#' Plot a forecasted Stochastic Mortality Model of class \code{"forStMoMo"}.
#' 
#' @usage 
#' \method{plot}{forStMoMo}(x, nCol = 2, parametricbx = TRUE, only.kt = FALSE,
#'  only.gc = FALSE, ...)
#' 
#' @param x an object of class \code{"forStMoMo"} with the forecast 
#'  of a stochastic mortality model.
#' @inheritParams plot.fitStMoMo
#' @param only.kt If \code{TRUE} only the period indexes of the model are 
#' plotted.
#' @param only.gc If \code{TRUE} only the cohort index of the model is 
#' plotted. This argument is ignored if \code{only.kt} is \code{TRUE}.
#'
#' @seealso \code{\link{plot.fitStMoMo}}
#'
#' @examples
#' wxt <- genWeightMat(55:89,  EWMaleData$years, clip = 3)
#' APCfit <- fit(apc(), Dxt = EWMaleData$Dxt, Ext = EWMaleData$Ext, 
#'               ages = EWMaleData$ages, years = EWMaleData$years, 
#'               ages.fit = 55:89, wxt = wxt)
#' APCfor <- forecast(APCfit)
#' plot(APCfor)
#' plot(APCfor, parametricbx = FALSE, nCol = 3)
#' plot(APCfor, only.kt = TRUE)
#' plot(APCfor, only.gc = TRUE, lwd = 2)
#' @export 
#' @method plot forStMoMo
plot.forStMoMo <- function(x, nCol = 2, parametricbx = TRUE, 
                           only.kt = FALSE, only.gc = FALSE, ...) {
  
  x.h <- x$model
  years.h <- x.h$years
  years.f <- x$kt.f$years
  ages <- x.h$ages
  cohorts.h <- x.h$cohorts
  cohorts.f <- x$gc.f$cohorts
  ax <- x.h$ax
  bx <- x.h$bx
  kt.h <- x.h$kt
  kt.f <- x$kt.f
  b0x <- x.h$b0x
  gc.h <- x.h$gc
  gc.f <- x$gc.f
  N <- x.h$model$N
    
  #Calculate number of plots and rows
  if (only.kt) 
    only.gc <- FALSE
  if (only.kt || only.gc) {
    nPlots <- N * only.kt + (!is.null(gc.h)) * only.gc
  } else {
    nPlots <- 2 * N + (!is.null(ax)) + 2 * (!is.null(gc.h))
    is.nonparametric <- function(bx) {
      is.character(bx) && bx == "NP"
    }  
    if (parametricbx == FALSE) {  #substract the parametric plots        
      nParametricbx <- ifelse(is.null(x.h$model$periodAgeFun), 0, 
                               sum(sapply(x.h$model$periodAgeFun, 
                                          FUN = function(x) !is.nonparametric(x)))) +
                        (!is.null(x.h$model$cohortAgeFun) && !is.nonparametric(x.h$model$cohortAgeFun))
                       
      nPlots <- nPlots - nParametricbx
    }   
    if (!is.null(ax) && nCol == 2 && parametricbx == TRUE) {
      nPlots <- nPlots + 1  # add and empty plot to the rigth of ax
    }  
  }

  nRow <- ceiling(nPlots / nCol)  
  oldpar <- par(no.readonly = TRUE)
  if (nPlots > 1L)
    par(mfrow = c(nRow, nCol))
  
  #ax
  if (!only.kt && !only.gc) {
    if (!is.null(ax)) {
      plot(x = ages, y = ax, ylab = "", xlab = "age", 
           main = expression(paste(alpha[x], " vs. x", "")), type = "l", ...)
    }
    if (!is.null(ax) && nCol == 2 && parametricbx == TRUE) {
      frame()  # add and empty plot to the rigth of ax
    }
  }
  # bx, kt
  if (N > 0) {
    for (i in 1:N) {      
      #bx
      if (!only.kt && !only.gc) {
        if (parametricbx == TRUE || is.nonparametric(x$model$periodAgeFun[[i]])) {
          plot(x = ages, y = bx[, i], ylab = "", xlab = "age", 
               main = substitute(paste(beta[x]^{(i)}, " vs. x", ""), 
                                 list(i = i)), type = "l", ...)
        
        }
      }
      #kt
      if (!only.gc) {
        kt.ylim <- range(kt.h[i, ], kt.f$mean[i, ], kt.f$upper[i, ], 
                         kt.f$lower[i, ], na.rm=TRUE)
        kt.xlim <- c(years.h[1],tail(years.f,1))
        plot(x = years.h,y = kt.h[i, ], ylab="", xlab = "year", 
             main = substitute(paste(kappa[t]^{(i)}, " vs. t", ""), 
                               list(i = i)), type = "l",
             xlim = kt.xlim, ylim = kt.ylim, ...)  
        lines(years.f, kt.f$mean[i, ], lty = 4, ...)
        lines(years.f, kt.f$upper[i, ], lty = 3, ...)
        lines(years.f, kt.f$lower[i, ], lty = 3, ...)
      }
    }
  }
  
  
  if (!is.null(gc.h) == TRUE) {
    #bx0
    if (!only.kt && !only.gc) {
      if (parametricbx == TRUE || is.nonparametric(x.h$model$cohortAgeFun)) {
        plot(x = ages, y = b0x, ylab = "", xlab = "age", 
             main = substitute(paste(beta[x]^{(i)}, " vs. x", ""), 
                               list(i = 0)), type = "l", ...)    
      }
    }
    #gc
    if (!only.kt) {
      gc.ylim <- range(gc.h, gc.f$mean, gc.f$upper, gc.f$lower, na.rm = TRUE)
      gc.xlim <- c(cohorts.h[1], tail(cohorts.f, 1))
      plot(x = cohorts.h, y = gc.h, ylab = "", xlab = "cohort", 
           main = expression(paste(gamma[t-x], " vs. t-x","")), type = "l", 
           xlim = gc.xlim, ylim = gc.ylim, ...)  
      lines(cohorts.f, gc.f$mean, lty = 4, ...)
      lines(cohorts.f, gc.f$upper, lty = 3, ...)
      lines(cohorts.f, gc.f$lower, lty = 3, ...)
    }
  }
  par(oldpar)
}


