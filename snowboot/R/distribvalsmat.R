distribvalsmat <- function(empd, vals, x) {
      if (any(vals == x)) {
            res <- empd[, which(vals == x)]
      } else {
            res <- rep(0, dim(empd)[1])
      }
      res
}

# empd 1 2 [1,] 0.5 0.5 [2,] 0.5 0.5 [3,] 0.5 0.5 [4,] 0.0 1.0 [5,] 0.5 0.5
# sapply(X=0:4,FUN=distribvalsmat,empd=empd,vals=vals) [,1] [,2] [,3] [,4] [,5] [1,] 0 0.5 0.5 0 0 [2,] 0 0.5 0.5 0 0
# [3,] 0 0.5 0.5 0 0 [4,] 0 0.0 1.0 0 0 [5,] 0 0.5 0.5 0 0
