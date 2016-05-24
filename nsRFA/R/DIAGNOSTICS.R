R2 <- function(x, y, na.rm=FALSE) {

 # INPUT
 # x = observed values
 # y = estimated values

 if (na.rm==TRUE) {
  y <- y[!is.na(x)]
  x <- x[!is.na(x)]
  x <- x[!is.na(y)]
  y <- y[!is.na(y)]
 }
 n <- length(x)
 if (!length(y)==n) stop("R2: x and y must have the same length")
 SST <- sum((x-mean(x))^2)
 SSRes <- sum((x-y)^2)
 R2 <- 1-SSRes/SST

 return(R2)
}


# ------------------------------------------------------------------ #

RMSE <- function(x, y, na.rm=FALSE) {

 # INPUT
 # x = observed values
 # y = estimated values

 if (na.rm==TRUE) {
  y <- y[!is.na(x)]
  x <- x[!is.na(x)]
  x <- x[!is.na(y)]
  y <- y[!is.na(y)]
 }
 n <- length(x)
 if (!length(y)==n) stop("RMSE: x and y must have the same length")
 res <- x-y
 RMSE <- sqrt(sum((res)^2)/n)

 return(RMSE)
}


# ------------------------------------------------------------------ #

MAE <- function(x, y, na.rm=FALSE) {

 # INPUT
 # INPUT
 # x = observed values
 # y = estimated values

 if (na.rm==TRUE) {
  y <- y[!is.na(x)]
  x <- x[!is.na(x)]
  x <- x[!is.na(y)]
  y <- y[!is.na(y)]
 }
 n <- length(x)
 if (!length(y)==n) stop("MAE: x and y must have the same length")
 res <- x-y
 MAE <- sum(abs(res))/n

 return(MAE)
}


# ------------------------------------------------------------------ #

RMSEP <- function(x, y, na.rm=FALSE) {

 # INPUT
 # x = observed values
 # y = estimated values

 if (na.rm==TRUE) {
  y <- y[!is.na(x)]
  x <- x[!is.na(x)]
  x <- x[!is.na(y)]
  y <- y[!is.na(y)]
 }
 n <- length(x)
 if (!length(y)==n) stop("RMSE: x and y must have the same length")
 res <- (x-y)/x
 RMSEP <- sqrt(sum((res)^2)/n)

 return(RMSEP)
}


# ------------------------------------------------------------------ #

MAEP <- function(x, y, na.rm=FALSE) {

 # INPUT
 # INPUT
 # x = observed values
 # y = estimated values

 if (na.rm==TRUE) {
  y <- y[!is.na(x)]
  x <- x[!is.na(x)]
  x <- x[!is.na(y)]
  y <- y[!is.na(y)]
 }
 n <- length(x)
 if (!length(y)==n) stop("MAE: x and y must have the same length")
 res <- (x-y)/x
 MAEP <- sum(abs(res))/n

 return(MAEP)
}

