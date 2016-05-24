table.row <- function(vecdat, vect) {
      # vecdat es el vector que se quiere arreglar segun los valores de vect
      table(c(vecdat, vect)) - 1
}
# example: m<-matrix(1:9,3,3) t(apply(m,1,table.row,1:9)) 1 2 3 4 5 6 7 8 9 [1,] 1 0 0 1 0 0 1 0 0 [2,] 0 1 0 0 1 0 0 1 0
# [3,] 0 0 1 0 0 1 0 0 1

# ---------------------------------------------------------------------------------------#

min.greater <- function(v, x, ge = TRUE) {
      # I construct this function to obtain percentiles from a empirical distribution.  This form is to allow obtaining it for
      # a matrix of distributions. v is a vector x is a value
      if (ge) {
            res <- min(which(v >= x))
      } else {
            res <- min(which(v > x))
      }
      res
}
