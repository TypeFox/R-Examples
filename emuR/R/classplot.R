##' Produce a classification plot from discriminant or SVM modelling
##' 
##' The function classifies all point specified within the ranges of xlim and
##' ylim based on the training model specified in model. It then produces a
##' two-dimensional plot colour-coded for classifications.
##' 
##' 
##' @param model A two-dimensional training model output from qda(), lda() of
##' MASS package , or svm() of e1071 package
##' @param xlim A vector of two numeric elements specifying the range on the
##' x-axis (parameter 1) over which classifications should be made
##' @param ylim A vector of two elements specifying the range on the y-axis
##' (parameter 2) over which classifications should be made
##' @param N A vector of one numeric element which specifies the density of
##' classification (greater N gives higher density). The default is 100.
##' @param pch A single element numeric vector specifying the plotting symbol
##' to be used in the classification plot. Defaults to 15.
##' @param col Either Null in which case the colours for the separate classes
##' are col = c(1, 2, ...n) where n is the number of classes; or else a vector
##' specifying the desired colours that is the same length as there are
##' classes.
##' @param legend A single element logical vector specifying whether a legend
##' should be drawn. Defaults to T
##' @param position A single element vector specifying the position in the
##' figure where the legend should be drawn. Defaults to "topright"
##' @param bg A single element vector specifying the background colour on which
##' the legend should be drawn.
##' @param ... Further arguments to plot.
##' @author Jonathan Harrington
##' @seealso \code{\link[MASS]{qda}}, \code{\link[MASS]{lda}}, svm of e1071
##' package. There is a function plot.svm which produces a prettier plot for
##' SVMs.
##' @examples
##' 
##' library(MASS)
##' # Data from female speaker 68
##' temp = vowlax.spkr=="68"
##' # Quadratic discriminant analysis
##' fm.qda = qda(vowlax.fdat.5[temp,1:2], vowlax.l[temp])
##' # Linear discriminant analysis
##' fm.lda = lda(vowlax.fdat.5[temp,1:2], vowlax.l[temp])
##' 
##' xlim=c(0,1000)
##' ylim=c(0,3000)
##' 
##' par(mfrow=c(1,2))
##' classplot(fm.qda, xlim=xlim, ylim=ylim, main="QDA")
##' classplot(fm.lda, xlim=xlim, ylim=ylim, main="LDA")
##' 
##' 
##' # install.packages("e1071")
##' # library(e1071)
##' # Support vector machine
##' \dontrun{fm.svm = svm(vowlax.fdat.5[temp,1:2], factor(vowlax.l[temp]))}
##' \dontrun{xlim = range(vowlax.fdat.5[temp,1])}
##' \dontrun{ylim = range(vowlax.fdat.5[temp,2])}
##' \dontrun{classplot(fm.svm, xlim=xlim, ylim=ylim, xlab="F1", ylab="F2", main="SVM")}          
##' 
##' @import MASS
##' @export classplot
`classplot` <- function(model, xlim, ylim, N = 100,  
                        pch=15, col=NULL, legend=TRUE, 
                        position="topright", bg="gray90",  ...)
{
  if(any(class(model) %in% "svm"))
  {
    priorlabels <- model$levels
    if(ncol(model$SV)!=2)
      stop("data must be two-dimensional")
  }
  else if(any(class(model) %in% c("lda", "qda")))
  {
    priorlabels <- rownames(model$means)
    if(ncol(model$means)!=2)
      stop("data must be two-dimensional")
  }
  pnts <- cbind(sort(rep(1:N/N, N)), rep(1:N/N, N))
  pnts[, 1] <- pnts[, 1] * (xlim[2] - xlim[1]) + xlim[1]
  pnts[, 2] <- pnts[, 2] * (ylim[2] - ylim[1]) + ylim[1]
  if(any(class(model) %in% c("lda", "qda")))
    blabs <- as.character(predict(model, pnts)$class)
  else if(any(class(model) %in% "svm"))
    blabs <- as.character(predict(model, pnts))
  k <- 1
  if(is.null(col))
    colours <- mu.colour(priorlabels, TRUE, FALSE)$colour
  else
    colours <- col
  graphics::plot(pnts, xlim=xlim, ylim=ylim, ...)
  for (j in priorlabels) {
    temp <- muclass(blabs, j)
    graphics::points(pnts[temp, ],pch=pch, col = colours[k])
    k <- k + 1
  }
  if(legend)
    graphics::legend(position, legend=priorlabels, col=colours, fill=colours, bg=bg)
}
