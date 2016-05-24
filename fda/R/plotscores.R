plotscores <- function(pcafd, scores = c(1, 2), xlab = NULL, ylab = NULL,
                       loc = 1, matplt2 = FALSE, ...)
{
#
#   Plot a scatter plot of the pca scores from a pca.fd object
#   If loc >0, you can then click on the plot in loc places and you'll get
#    plots of the functions with these values of the principal component
#    coefficients.
#
#  The present implementation doesn't work for multivariate functional data
#
#    pcafd      a pca.fd object
#    scores     a two dimensional vector giving the indices of the two
#                  scores to be plotted; if scores is a single number then
#                  that score will be plotted against component 1; the default
#                  is to print the first two components.
#    xlab, ylab   labels for the principal components scores scatterplot
#    loc   if an integer, the number of sample functions to be plotted.
#                  This number of clicks on the first plot are needed.
#      if a list with components x and y, the coordinates of the
#         functions to be plotted (the output from a previous call
#                  of plotscores, for instance).  No prompting will be done.
#               if 0 or NULL nothing is plotted beyond the scatterplot.
#
#    matplt2    the matplt value for the plot of the sample functions;
#               if matplt=TRUE, the curves are plotted on the same plot;
#               if matplt=FALSE, they are plotted separately.
#
#   ...         arguments to be passed to the pc scores scatterplot
#
#  RETURNS:   a list containing the PC scores of the plotted functions

#  Last modified 26 October 2005

   if (!(inherits(pcafd, "pca.fd"))) stop('Argument PCAFD is not a pca.fd object.')

   if(length(scores) == 1) scores <- c(1, scores)
   if(length(scores) != 2)
      scores <- c(1, 2)
   if(max(scores) > dim(pcafd$harmonics$coefs)[2]) {
      stop(paste("The pca.fd object does not contain ", max(scores),
         "components"))
   }
   if(is.null(xlab))
      xlab <- paste("PCA score ", scores[1])
   if(is.null(ylab))
      ylab <- paste("PCA score ", scores[2])
   plot(pcafd$scores[, scores], xlab = xlab, ylab = ylab, ...)
   if(is.list(loc))
      zz <- loc
   else {
      if(is.na(loc) || is.null(loc) || loc == 0)
         return(NULL)
      zz <- locator(loc)
   }
   zzmat <- rbind(zz$x, zz$y)
   coefs <- pcafd$meanfd$coefs %*% rep(1, dim(zzmat)[2]) + pcafd$harmonics$
      coefs[, scores] %*% zzmat
   fdnames <- pcafd$meanfd$fdnames
   fdnames[[2]] <- paste("Score", scores[1], "=", signif(zz$x, 2),
      "; Score", scores[2], "=", signif(zz$y, 2))
   names(fdnames)[2] <- "Sample function"
   names(fdnames)[3] <- "Function value"
#   fd <- create.fd(coefs, pcafd$meanfd$basis, fdnames)
   fd <- fd(coefs, pcafd$meanfd$basis, fdnames)
   plot(fd, matplt = matplt2)
   return(zz)
}
