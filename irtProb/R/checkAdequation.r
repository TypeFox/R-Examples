`checkAdequation`<-
function(x) {
 res    <- c(isNumbers        = FALSE,
             correctClass     = FALSE,
             squareMatrix     = FALSE,
             diagPositive     = FALSE,
             positiveDefinite = FALSE,
             nonSingular      = FALSE)
 xVector <- as.numeric(x)
 res[1]  <- !(any(is.infinite(xVector)) || any(is.nan(xVector)) || any(is.na(xVector)))
 if (is.data.frame(x) || is.matrix(x)) {
  res[2] <- TRUE
  if (res[1] == TRUE) {
   if (dim(x)[1] == dim(x)[2]) {
    res[3] <- TRUE
    res[4] <- !all(diag(x) < 0)
    res[5] <- all((eigen(x)$values) > 0) # Positive definiteness check: TRUE
    res[6] <- abs(det(x)) > 0            # Non sigularity check: (det(x) != 0
    }
   }
  }
 return(res)
 }




