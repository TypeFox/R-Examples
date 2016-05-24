rbfKernGradX <-
function (kern, x1, x2) {

  gX = array(0, c(dim(as.array(x2))[1], dim(as.array(x2))[2], dim(as.array(x1))[1]))
  for (i in 1:dim(x1)[1]) {
    gX[, , i] = rbfKernGradXpoint(kern, x1[i, ], x2)
  }

  return (gX)
}
