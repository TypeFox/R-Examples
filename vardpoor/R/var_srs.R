
var_srs <- function(Y, w) {
 
  ### Checking
  # Y
  Y <- data.table(Y, check.names=TRUE)
  n <- nrow(Y)   
  if (!all(sapply(Y, is.numeric))) stop("'Y' must be numerical")
  if (any(is.na(Y))) print("'Y' has unknown values")
  if (is.null(colnames(Y))) stop("'Y' must be colnames")
    
  # w 
  w <- data.frame(w)
  if (nrow(w) != n) stop("'w' must be equal with 'Y' row count")
  if (ncol(w) != 1) stop("'w' must be vector or 1 column data.frame, matrix, data.table")
  w <- w[, 1]
  if (!is.numeric(w)) stop("'w' must be numerical")
  if (any(is.na(w))) stop("'w' has unknown values")
    
  ### Calculation

  # N
  Nn <- sum(w)
  
  konst <- Nn^2 * (1 - n/Nn) / n

  wyN <- as.numeric(Y[, lapply(.SD, function(x) sum(x*w)/Nn)])

  varsrs <- data.table(t(as.numeric(konst*Y[, mapply(function(x, x.vid) sum(w * (x - x.vid)^2),
           .SD, wyN)] / (Nn - 1))))
  setnames(varsrs, names(varsrs), names(Y))

  varsrs2 <- data.table(t(as.numeric(konst*(Y[, lapply(.SD, function(x) sum(w * x^2))] / (Nn - 1) -
    Nn / (Nn - 1) * (wyN)^2))))
  setnames(varsrs2, names(varsrs2), names(Y))

 return(varsrs2)
}

