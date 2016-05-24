
qq_get_p <- function(x, method=1, ties.method = "average")
{
  if (method < 1 | method > 7) 
    stop("'method' must be a number between 1 and 7:", 
         "\n\t1 = Blom \n\t2 = Rankit / Hazen \n\t3 = Tukey",
         "\n\t4 = van de Waerden / Weibull \n\t5 = Benard and Bos-Levenbach",
         "\n\t6 = Gringorten\n\t7 = Yu and Huang", call. = FALSE) 
  n <- length(x) 
  #i <- order(order(x))   # to recreate original order from sorted data
  #xs <- sort(x)
  r <- rank(x, ties.method = ties.method)
  p <- switch(method,
              "1" = (r - 3/8) / (n + 1/4),      # Blom
              "2" = (r - 1/2) / n,              # Rankit / Hazen
              "3" = (r - 1/3) / (n + 1/3),      # Tukey   
              "4" = r / (n + 1),                # van de Waerden / Weibull
              "5" = (r - 3/10) / (n + 4/10),    # Benard and Bos-Levenbach
              "6" = (r - 0.44) / (n + 0.12),    # Gringorten
              "7" = (r - 0.326) / (n + 0.348))  # Yu and Huang 
  #p[i]
  p
}


#' SPSS like QQ-plot
#'
#' The QQ-plot in SPSS and R looks very different. The points
#' points and the QQ-line are positioned differently. 
#' \code{qqnorm_spss} implements a version of the QQ-plot that resembles
#' the SPSS version. The function returns an object containing the 
#' processed data. The output can be plotted using the function \code{plot} 
#' and \code{ggplot}. The parameters that can be passed to the 
#' plotting functions are documented in \code{\link{plot.qqnorm.spss}} and
#' \code{\link{ggplot.qqnorm.spss}}.
#' 
#' @param x A numeric vector.
#' @param standardize Whether the quantiles of the standardized values
#'  should be displayed. The default is to display the quantiles using the
#'  original data.
#' @param method The method used to assign probabilties for the
#'  ranks that are then converted into quantiles.   
#'  The following methods are implemented (see Castillo-Gutiérrez, Lozano-Aguilera, 
#'  & Estudillo-Martínez, 2012): 
#'  \code{1 =} Blom (default), \code{2 =} Rankit / Hazen, \code{3 =} Tukey,
#'  \code{4 =} Van der Waerden / Weibull, \code{5 =} Benard and Bos-Levenbach,
#'  \code{6 =} Gringorten and \code{7 =} Yu and Huang.
#' @param  ties.method Method to assign ranks to ties. One of 
#'  \code{"average", "first", "random", "max", "min"}. See \code{ties.method} 
#'  argument from \code{\link{rank}} for more details.
#' 
#' @return An list object of class \code{qqnorm.spss} with the 
#'  following elements:
#'  \item{x}{The orginal data}
#'  \item{y}{Corresponding quantiles in original scaling}
#'  \item{x.std}{Standardized values}
#'  \item{y.std}{Corresponding quantiles for standardized values}
#'  \item{method.name}{Name of the method to assign probabilities to ranks}
#'  \item{ties.method}{Method to treat ties}
#'  \item{xname}{Name of the variable used to produce the plot}
#' @references 
#'  Castillo-Gutiérrez, S., Lozano-Aguilera, E., & Estudillo-Martínez, M. D. (2012). 
#'    Selection of a Plotting Position for a Normal Q-Q Plot. R Script. 
#'    \emph{Journal of Communication and Computer, 9}(3), 243–250.
#' @export
#' @section TODO:
#' Check output against SPSS results. 
#' @example inst/examples/example-qq-plot.R
#' 
qqnorm_spss <- function(x, standardize=FALSE, method=1, 
                        ties.method="average") 
{ 
  xname <- deparse(substitute(x))
  x <- na.omit(x)
  methods <- c('Blom'=1, 'Rankit / Hazen'=2, 'Tukey'=3, 'Van der Waerden / Weibull'=4,
               'Benard and Bos-Levenbach'=5, 'Gringorten'=6, 'Yu and Huang'=7)
  method.name <- names(methods[method])
  p <- qq_get_p(x, method=method, ties.method=ties.method)
  y.std <- qnorm(p) 
  x.std <- as.vector(scale(x))
  y <- y.std * sd(x) + mean(x)
  l <- list(x=x, y=y, x.std=x.std, y.std=y.std, 
            method.name=method.name,
            standardize=standardize,
            ties.method=ties.method,
            xname=xname)
  class(l) <- "qqnorm.spss"
  l
} 


#' Plot the output from \code{qqplot.spss}
#' 
#' @param x An object as returned by \code{\link{qqnorm_spss}}
#' @param plottype The type of plot created. 
#'  \code{1 =} Standard QQ-plot, \code{2 =} Detrended QQ-plot.
#' @param line Whether to plot a QQ-line (defaul is \code{TRUE})
#' @param l.col Color of the QQ-line.
#' @param ... Passed to \code{plot} method.
#' @export
#' @keywords internal
#' 
plot.qqnorm.spss <- function(x, plottype=1, line=TRUE,
                             l.col="black", ...) 
{
  qq <- x
  x <- qq$x
  y <- qq$y
  main <- paste("Normal Q-Q plot of", qq$xname) 
  xlab <- "Observed value"
  ylab <- "Expected normal value"
  if (qq$standardize) {
    x <- qq$x.std
    y <- qq$y.std
    xlab <- "Standardized observed value"
  }  
  if (plottype == 2) {        # convert to detrended data
    main <-  paste("Detrended normal Q-Q plot of", qq$xname) 
    ylab <- "Deviation from normal"
    y <- x - y
  }
  plot(x, y, main=main,   # cex=.8, pch=16, col="black"
       xlab=xlab, ylab=ylab, ...)
  if (line) {
    if (plottype == 2)          # detrended plot
      abline(h=0, col=l.col)    # zero line is shown
    else                        # standard plot
      abline(0, 1, col=l.col)   # slope of 1
  }  
}


#' Plot the output from \code{qqplot.spss} using \code{ggplot2}
#' 
#' @param x An object as returned by \code{\link{qqnorm_spss}}
#' @param plottype The type of plot created. 
#'  \code{1 =} Standard QQ-plot, \code{2 =} Detrended QQ-plot.
#' @param line Whether to plot a QQ-line (defaul is \code{TRUE})
#' @param l.col Color of the QQ-line.
#' @param ... Not evaluated.
#' @return A ggplot object.
#' @export
#' @keywords internal
#' 
ggplot.qqnorm.spss <- function(x, plottype=1, line=TRUE,
                      l.col="black", ...) 
{
  qq <- x
  x <- qq$x
  y <- qq$y
  main <- paste("Normal Q-Q plot of", qq$xname) 
  xlab <- "Observed value"
  ylab <- "Expected normal value"
  if (qq$standardize) {
    x <- qq$x.std
    y <- qq$y.std
    xlab <- "Standardized observed value"
  }  
  if (plottype == 2) {        # convert to detrended data
    main <-  paste("Detrended normal Q-Q plot of", qq$xname) 
    ylab <- "Deviation from normal"
    y <- x - y
  }
  d <- data.frame(x, y)
  g <- ggplot2::ggplot(data=d, aes(x,y)) + 
       ggplot2::geom_point() + 
       ggplot2::xlab(xlab) + ggplot2::ylab(ylab) + ggplot2::ggtitle(main)
  
  if (line) {
    if (plottype == 2)          # detrended plot
      gline <- ggplot2::geom_abline(intercept = 0, slope=0, colour=l.col)
    else                        # standard plot
      gline <- ggplot2::geom_abline(intercept=0, slope=1, colour=l.col)   # slope of 1
    g <- g + gline
  }  
  g
}

