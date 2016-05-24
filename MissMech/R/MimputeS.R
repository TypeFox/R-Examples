MimputeS <- function(data, patused, y1, s1, e)
{
# This function imputes the missing data by Srivastava method,
# given the observed and ybar and s from complete data set
# for a single pattern it goeas through each patterns and uses the
# linear regresion to predict missing given obsrved data and add the
# residual to impute missing data 
 x <- data
 ni <- nrow(x)
 pp <- ncol(x)
 indm <- which(is.na(patused))
 indo <- which(!is.na(patused))
 pm <- length(indm)
 po <- length(indo)
 a <- matrix(s1[indm, indo], pm, po) %*% solve(s1[indo, indo])
 dif <-  x[, indo] - matrix(y1[indo], ni, po, byrow = TRUE)
 z <- matrix(y1[indm], ni, pm, byrow = TRUE) + dif %*% t(a)
 etta <- matrix(e[, indm], ni, pm) - matrix(e[, indo], ni, po) %*% t(a)
 zij <- z + etta
 x[, indm] <- zij 
 x
}