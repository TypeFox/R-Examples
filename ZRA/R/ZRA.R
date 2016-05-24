#' @name ZRA
#'
#' @title Dynamic Plots for Time Series Forecasting
#'
#' @param data (ts) Time series data.
#' @param FP (numeric) Forecast period.
#' @param SL (numeric) significance levels.
#' @param ... further arguments passed to the forecast-function.
#'
#' @return An Object of class "ZRA".
#'
#' @examples zra <- ZRA(fdeaths)
#' plot(zra, zero = TRUE)
#'
#' @seealso forecast, dygraphs
#'
#' @import forecast
#' @import dygraphs
#' @importFrom stats end frequency ts is.ts
#'
#' @export

ZRA <- function(data,FP = 10, SL = c(0.80,0.95), ...) {

  startvalue <- end(data)[1]+1
  f  <- frequency(data)
  d <- data


  if (is.ts(d)== TRUE) {


    if (is.matrix(d)==TRUE) {
      result <- NULL
      stop("Only 1 Time Series can analyzed at once")
    }


    else {

      prognose <- forecast(d, h = FP, level=SL, ...)

      result <- list()
      result$series <- data
      result$SL <- SL
      result$FP <- FP

      if (length(SL)==1) {

        up1 <- ts(prognose$upper[,1],start=startvalue,frequency = f)
        low1 <- ts(prognose$lower[,1],start=startvalue,frequency = f)
        fit1 <- (up1 + low1)/2

        result$up1 <- up1
        result$low1 <- low1
        result$fit1 <- fit1
        result$piv1 <- cbind(fit1,up1,low1)

      }

      if (length(SL)==2) {

        up1 <- ts(prognose$upper[,1],start=startvalue,frequency = f)
        low1 <- ts(prognose$lower[,1],start=startvalue,frequency = f)
        fit1 <- (up1 + low1)/2

        up2 <- ts(prognose$upper[,2],start=startvalue,frequency = f)
        low2 <- ts(prognose$lower[,2],start=startvalue,frequency = f)
        fit2 <- (up2 + low2)/2

        result$up1 <- up1
        result$low1 <- low1
        result$fit1 <- fit1
        result$piv1 <- cbind(fit1,up1,low1)

        result$up2 <- up2
        result$low2 <- low2
        result$fit2 <- fit2
        result$piv2 <- cbind(fit2,up2,low2)

      }

      if (length(SL)!=1 & length(SL)!=2 ) {
        stop("Only 2 Significance levels can be plotted at once.")

      }
    }
  }

  else {
    result <- NULL
    stop("Data have to be a Time Series Obejct, with the Class ts.")
  }

  class(result) <- "ZRA"

  return(result)

}



#' @name print.ZRA
#'
#' @title Printig a ZRA-Class Obejct
#'
#' @description This print-Method will print out the original Time Series and the predicted Data.
#'
#' @param x (ZRA) An ZRA-Object.
#' @param ... further arguments passed to the print method.
#'
#' @method print ZRA
#'
#' @usage \method{print}{ZRA}(x, ...)
#'
#' @export

print.ZRA <- function(x, ...){

  sl <- as.character(x$SL)

  if (length(x$SL)==1) {
    colnames(x$piv1) <- c(paste("Fit",sl,sep = "-"),paste("Up",sl,sep = "-"),paste("Low",sl,sep = "-"))
    cat("Time Series: \n")
    print(x$series)
    cat("\n Forecast: \n")
    print(round(x$piv1,3))

  }

  if (length(x$SL)==2) {

    tabelle <- cbind(x$piv1,x$piv2)

    colnames(tabelle) <- c(paste("Fit",sl[1],sep = "-"),paste("Up",sl[1],sep = "-"),paste("Low",sl[1],sep = "-"),paste("Fit",sl[2],sep = "-"),paste("Up",sl[2],sep = "-"),paste("Low",sl[2],sep = "-"))
    cat("Time Series: \n")
    print(x$series)
    cat("\n Forecast: \n")
    print(round(tabelle,3))

  }

  if (length(x$SL)!=1 & length(x$SL)!=2 ) {
    stop("Only 2 Significance levels can be plotted at once.")

  }

}

#' @name plot.ZRA
#'
#' @title Ploting a ZRA-Class Obejct
#'
#' @description This plot-Method will plot the original Time Series and the predicted Data, using dygraphs.
#'
#' @param x (ZRA) An ZRA-Object.
#' @param zero (boolean) If zero=TRUE, they will be the 0 included in the Graph.
#' @param ... further arguments passed to the plot method.
#'
#' @method plot ZRA
#'
#' @usage \method{plot}{ZRA}(x,zero, ...)
#'
#' @export

plot.ZRA <- function(x, zero =TRUE, ...) {

  if (length(x$SL)==1 || length(x$SL)==2) {

    result <- ZRAplot(x, ...)
    return(result)

  }

  if (length(x$SL)!=1 & length(x$SL)!=2 ) {
    stop("Only 2 Significance levels can be plotted at once." )

  }

}

ZRAplot <- function(x,zero=TRUE, ...) {

  sl <- as.character(x$SL)
  h <- as.character(x$FP)

  if (length(x$SL)==1) {

    plotreihe1 <- cbind(x$series, x$piv1, x$up1, x$low1)

    plot1 <- dygraph(plotreihe1)  %>%
      dySeries("x$series",label = "Time Series",color="blue",fillGraph = TRUE) %>%
      dyLimit(limit=0) %>%

      dySeries("x$up1",label = "Upper Interval limit",color="chocolate",strokePattern = "dotted") %>%
      dySeries("x$low1",label = "Lower Interval limit",color="chocolate",strokePattern = "dotted") %>%

      dySeries(c("x$piv1.up1","x$piv1.fit1","x$piv1.low1"),label="Point Estimator", color="red")  %>%

      dyOptions(includeZero = TRUE,fillAlpha =0.1 )  %>%
      dyRangeSelector()

    print(paste("Significance level:",sl))
    print(paste("Prediction interval:",h,"Periods"))
    return(plot1)

  }

  if (length(x$SL)==2) {

    plotreihe2 <- cbind(x$series, x$piv1, x$up1, x$low1, x$piv2, x$up2, x$low2)

    plot2 <- dygraph(plotreihe2)  %>%
      dySeries("x$series",label = "Time Series",color="blue",fillGraph = TRUE) %>%
      dyLimit(limit=0) %>%

      dySeries("x$up1",label = paste("Up",sl[1]),color="chocolate",strokePattern = "dotted") %>%
      dySeries("x$low1",label = paste("Low",sl[1]),color="chocolate",strokePattern = "dotted") %>%

      dySeries("x$up2",label = paste("Up",sl[2]),color="darkseagreen",strokePattern = "dotted") %>%
      dySeries("x$low2",label = paste("Low",sl[2]),color="darkseagreen",strokePattern = "dotted") %>%

      dySeries(c("x$piv1.up1","x$piv1.fit1","x$piv1.low1"),label=paste("Point Estimator",sl[1]), color="red")  %>%
      dySeries(c("x$piv2.up2","x$piv2.fit2","x$piv2.low2"),label=paste("Point Estimator",sl[2]), color="green")  %>%

      dyOptions(includeZero = zero,fillAlpha =0.1 )  %>%
      dyRangeSelector()

    print(paste("Significance levels:",sl[1],"(red),",sl[2],"(green)"))
    print(paste("Prediction interval:",h,"Periods"))
    return(plot2)

  }

}
