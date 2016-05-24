func.extendDataframeForRegValues <-
function( extFrame, coeff ){
  tmp <- NULL
  tmp$noOfRows <- NCOL(extFrame)
  tmp$intercept <- 0
  
  for(i in 1:NROW(coeff)){
    param <- letters[i - tmp$intercept]
    extFrame[FALSE, i+tmp$noOfRows] <- numeric(0)
    
    if(names(coeff)[i] == "(Intercept)"){
      param <- "const"
      increment(tmp$intercept)
    }
  
    colnames(extFrame)[i+tmp$noOfRows] <- param
  }
  
  return( extFrame )
}

