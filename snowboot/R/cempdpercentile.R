cempdpercentile <- function(cempd, perc, vals) {
      # cempd is a matrix perc is a scalar
      res <- vals[apply(X = cempd, 1, FUN = min.greater, x = perc, ge = TRUE)]
      res
}
# Ejemplo cempdpercentile(cempd,0.25,vals)

# sapply(X=c(.25,.5,.75),FUN=cempdpercentile,cempd=cempd,vals=vals)
# sapply(X=c(.25,.5,.75),FUN=cempdpercentile,cempd=cempd,vals=vals) [,1] [,2] [,3] [1,] 1 1 2 [2,] 1 1 2 [3,] 1 1 2 [4,]
# 2 2 2 [5,] 1 1 2
