#-----------------------------------------------------------------------------------------------------------------#
# Function plotLSVM
#-----------------------------------------------------------------------------------------------------------------#

# plot the data, the hyperplanes and the estimate of the limit state
#' @export
plotLSVM <- function(X, Y, A.model.lsvm, hyperplanes = FALSE, limit.state.estimate = TRUE, convexity){

  rA <- ifelse( is.null(dim(A.model.lsvm)), 1, dim(A.model.lsvm)[1])
  x <- seq(from = min(X[,1]), to = max(X[,1]), length = 1e3)
  y <- matrix(0, ncol = rA, nrow = 1e3)

  XS <- X[which(Y > 0), ]
  XF <- X[which(Y <= 0), ]

  if(convexity == -1){
    plot(XF[,1], XF[,2], col= "blue", xlim = c(min(X[,1]), max(X[,1])), ylim = c(min(X[,2]), max(X[,2])), pch=19, xlab= "", ylab="")
  }else{
    plot(XS[,1], XS[,2], col= "red",xlim = c(min(X[,1]), max(X[,1])), ylim = c(min(X[,2]), max(X[,2])), pch=19, xlab= "", ylab="")
  }
  for(j in 1:rA){
    y[, j] <- -(A.model.lsvm[j, 3] + A.model.lsvm[j,1]*x)/A.model.lsvm[j, 2]
    #plot hyperplanes
    if(hyperplanes){
      lines(x, y[,j], col= gray(0.65))
    }
  }

  if(convexity == - 1){
    points(XS[,1], XS[,2], col= "red", pch=19)
  }else{
    points(XF[,1], XF[,2], col= "blue", pch=19)
  }

  #plot an estimate of the limit state
  if(limit.state.estimate){
    x.classifier.min <- apply(y, MARGIN = 1, min)
    lines(x, x.classifier.min, col= "purple", lwd = 3)
  }
}






