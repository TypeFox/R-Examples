#' Extract deviance residuals of a Stochastic Mortality Model
#' 
#' Compute deviance residuals of a fitted Stochastic Mortality Model. 
#' These residuals can be plotted using \code{\link{plot.resStMoMo}}.
#' 
#' @param object an object of class \code{"fitStMoMo"} with the fitted 
#' parameters of a stochastic mortality model.
#' @param scale logical indicating whether the residuals should be scaled 
#' or not by dividing the deviance by the  overdispersion of the model.  
#' Default is \code{TRUE}.
#' @param ... other arguments.
#' 
#' @return An object of class \code{"resStMoMo"} with the residuals. This 
#' object has components:
#'   \item{residuals}{ a matrix with the residuals.}
#'   \item{ages}{ ages corresponding to the rows in \code{residuals}.}
#'   \item{years}{ years corresponding to the columns in \code{residuals}.}
#' 
#' @seealso \code{\link{plot.resStMoMo}}
#' 
#' @examples
#' CBDfit <- fit(cbd(), Dxt = EWMaleData$Dxt, Ext = EWMaleData$Ext, 
#'               ages = EWMaleData$ages, years = EWMaleData$years,
#'               ages.fit = 55:89)
#' CBDres <- residuals(CBDfit)
#' plot(CBDres)
#' @export 
residuals.fitStMoMo <- function(object, scale = TRUE, ...) {
  D <- object$Dxt + 0.00001 #Add a small amount to compensate for the 
                            #possibility of 0 deaths
  W <- object$wxt
  ind <- (W > 0)
  res <- array(NA, dim(W))
  if (object$model$link == "log") {
    Dhat <- fitted(object, type = "deaths")
    res[ind] <- 2 * W[ind] * (D[ind] * log(D[ind] / Dhat[ind]) - (D[ind] - Dhat[ind]))
    signRes <- sign(D - Dhat)
  } else if (object$model$link == "logit") {
    E <- object$Ext
    qxhat <- fitted(object, type = "rates")
    qx <- D / E
    res[ind] <- 2 * W[ind] * E[ind] *(qx[ind] * log(qx[ind] / qxhat[ind]) + 
          (1 - qx[ind]) * log((1 - qx[ind]) / (1 - qxhat[ind]))) 
    signRes <- sign(qx - qxhat)
  }
  if (scale) 
    phi <- sum(res[ind]) / (object$nobs - object$npar)
  else 
    phi <- 1
  res <- signRes * sqrt(abs(res) / phi) 
  colnames(res) <- object$years
  rownames(res) <- object$ages
  structure(list(residuals = res, ages = object$ages, years = object$years),
            class = "resStMoMo")  
}


#' Plot the residuals of a Stochastic Mortality Model
#' 
#' Plots the deviance residuals of a Stochastic Mortality Model which are 
#' of class \code{"resStMoMo"}. Three types of plots
#' are available: scatter plot of residuals by age, period and cohort,
#' colour map (heatmap) of the residuals, and a black and white signplot 
#' of the residuals.
#' 
#' @usage \method{plot}{resStMoMo}(x, type = c("scatter", "colourmap", 
#'                                             "signplot"), 
#'                                 reslim = NULL, plotAge = TRUE, 
#'                                 plotYear = TRUE, plotCohort  = TRUE, 
#'                                 pch = 20, col = NULL, ...)
#' 
#' @param x an object of class \code{resStMoMo} with the residuals of a 
#' Stochastic Mortality Model.
#' @param type the type of the plot. The alternatives are 
#' \code{"scatter"}(default), \code{"colourmap"}, and \code{"signplot"}.
#' @param reslim optional numeric vector of length 2, giving the range of the 
#' residuals.
#' @param plotAge logical value indicating if the age scatter plot should be 
#' produced. This is only used when \code{type = "scatter"}.
#' @param plotYear logical value indicating if the calendar year scatter plot 
#' should be produced. This is only used when \code{type = "scatter"}.
#' @param plotCohort logical value indicating if the cohort scatter plot 
#' should be produced. This is only used when \code{type = "scatter"}.
#' @param pch optional symbol to use for the points in a scatterplot. 
#' This is only used when \code{type = "scatter"}. See 
#' \code{\link[graphics]{plot}}.
#' @param col optional colours to use in plotting. If 
#' \code{type = "scatter"} this is a single colour to use in the points
#' in the scatter plots, while if \code{type = "colourmap"} this should
#' be a list of colours (see help in \code{\link[fields]{image.plot}} 
#' for details). This argument is ignored if \code{type = "signplot"}.
#' @param ... other plotting parameters to be passed to the plotting 
#' functions. This can be used to control the appearance of the plots.
#'
#' @details
#' When \code{type = "scatter"} scatter plots of the residuals against age, 
#' calendar year and cohort (year of birth) are produced. 
#'
#' When \code{type = "colourmap"} a two dimensional colour map of the 
#' residuals is plotted. This is produced using function 
#' \code{\link[fields]{image.plot}}. See \code{\link[fields]{image.plot}} 
#' for further parameters that can be passed to this type of plots.
#'
#' When \code{type = "signplot"} a two dimensional black and white map of the
#'  residuals is plotted with dark grey representing negative residuals and 
#'  light grey representing positive residuals. This is produced using 
#'  function \code{\link[graphics]{image.default}}. 
#'   
#'  @seealso \code{\link{residuals.fitStMoMo}}
#'   
#' @examples
#' CBDfit <- fit(cbd(), Dxt = EWMaleData$Dxt, Ext = EWMaleData$Ext, 
#'               ages = EWMaleData$ages, years = EWMaleData$years,
#'               ages.fit = 55:89)
#' CBDres <- residuals(CBDfit)
#' plot(CBDres)
#' plot(CBDres, type = "signplot")
#' plot(CBDres, type = "colourmap")
#'   
#' @export 
#' @method plot resStMoMo
plot.resStMoMo <- function(x, type = c("scatter", "colourmap", "signplot"), 
                           reslim = NULL, plotAge = TRUE, plotYear = TRUE, 
                           plotCohort  = TRUE, pch = 20, col = NULL, ...) {
  type <- match.arg(type)
  oldpar <- par(no.readonly = TRUE)
  
  if (is.null(reslim)) {
    maxRes <- max(abs(x$residuals), na.rm = TRUE)
    reslim <- c(-maxRes, maxRes)
  }
  if (is.null(col) & type == "colourmap") {
    col <- colorRampPalette(RColorBrewer::brewer.pal(10, "RdBu"))(64)
  }
  if (is.null(col) & type == "scatter") {
    col <- "black"
  }
  
  switch(type, 
         scatter = scatterplotAPC(x$residuals, x$ages, x$years, 
                                  plotAge = plotAge, plotYear = plotYear, 
                                  plotCohort = plotCohort, pch = pch, 
                                  ylab = "residuals", ylim = reslim, 
                                  col = col, ...),
         colourmap = fields::image.plot(x$year, x$age, t(x$residuals), 
                                        zlim = reslim, ylab = "age", 
                                        xlab = "calendar year", col = col, 
                                        ...),
         signplot = image.default(x$year, x$age, t(x$residuals), 
                                  zlim = reslim, ylab = "age", 
                                  xlab = "calendar year", 
                                  breaks = c(-10e10, 0, 10e10), 
                                  col = grey.colors(2), ...)
         )
  par(oldpar)
}


#'Do a scatter plot of a matrix according to age-period-cohorts
#'
#' @param mat matrix with the data to plot.
#' @param ages ages corresponding to the rows in \code{mat}.
#' @param years years corresponding to the columns in \code{mat}.  
#' @param plotAge logical value indicating if the age scatter plot should be 
#' produced.
#' @param plotYear logical value indicating if the calendar year scatter plot 
#' should be produced.
#' @param plotCohort logical value indicating if the cohort scatter plot 
#' should be produced.
#' @param zeroLine logical valae indicating if a horizontal line at zero
#' should be plotted.
#' @param ... other arguments to pass to the plot function.
#' @keywords internal
scatterplotAPC <- function(mat, ages, years, plotAge = TRUE, plotYear = TRUE, 
                           plotCohort  = TRUE, zeroLine  = TRUE, ...) {
  nAges <- length(ages)
  nYears <- length(years)  
  cohorts <- (years[1] - ages[nAges]):(years[nYears] - ages[1])
  nCohorts <- length(cohorts)
  
  mat <- as.matrix(mat)
  if ( nrow(mat) != nAges ||  ncol(mat) != nYears) {
    stop( "Mismatch between the dimensions of the plottin data and the 
          number of years or ages")
  }
  rownames(mat) <- ages
  colnames(mat) <- years
  data <- (reshape2::melt(mat, value.name = "y", varnames = c("x", "t")))
  x <- NULL #hack to remove note in CRAN check
  data <- transform(data, c = t - x) 
  
  N <- plotAge + plotYear + plotCohort
  if (N > 0) par(mfrow=c(1, N))
  
  if (plotAge) {
    plot(data$x, data$y, type = "p", xlab = "age", ...)
    if (zeroLine) 
      abline(h = 0)
  }
  if (plotYear) {
    plot(data$t, data$y, type = "p", xlab = "calendar year", ...)
    if (zeroLine) 
      abline(h = 0)
  }
  
  #cohort 
  if (plotCohort) {
    plot(data$c, data$y, type = "p", xlab = "year of birth", ...)
    if (zeroLine) 
      abline(h = 0)    
  }   
}
