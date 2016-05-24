align <- function(wt, coe=FALSE, inverse=FALSE){

  # error checking
  if(!inverse){
    if(wt@aligned) stop("Unnecessary proceedure: 'wt' object is already aligned.")
  } else {
    if(!wt@aligned) stop("Unnecessary proceedure: 'wt' object is already unaligned.")
  }
  if(is.na(match(class(wt), c("dwt", "modwt"))))
    stop("Invalid argument: 'wt' must be of class 'dwt' or 'modwt'")
  
  filter <- wt@filter
  if(filter@transform == "modwt") modwt <- TRUE
  J <- length(wt@W)

  # a function to align the wavelet and scaling coefficients
  align.coef <- function(wt.coef, wavelet, coe, modwt, inverse){
    shift.coef <- lapply(1:J, function(j,coef,filter,wavelet,coe,modwt,inverse){
      shift <- wt.filter.shift(filter, j, wavelet, coe, modwt)
      N <- dim(coef[[j]])[1]
      if(shift >= N) shift <- shift - floor(shift/N)*N
      if(shift == 0){
        shiftcoef <- coef[[j]]
      } else {
        if(dim(coef[[j]])[2] == 1){
          if(!inverse)
            shiftcoef <- as.matrix(c(coef[[j]][(shift+1):N,], coef[[j]][1:shift,]))
          else
            shiftcoef <- as.matrix(c(coef[[j]][(N-shift+1):N,], coef[[j]][1:(N-shift),]))
        } else {
          if(!inverse)
            shiftcoef <- rbind(coef[[j]][(shift+1):N,], coef[[j]][1:shift,])
          else{
            shiftcoef <- rbind(coef[[j]][(N-shift+1):N,], coef[[j]][1:(N-shift),])
          }
        }
      }
      return(shiftcoef)
    }, coef=wt.coef, filter=filter, wavelet=wavelet, coe=coe, modwt=modwt, inverse=inverse)
    return(shift.coef)
  }

  # align the coefficients
  wt.shifted <- wt
  wt.shifted@W <- align.coef(wt@W, wavelet=TRUE, coe=coe, modwt=modwt, inverse=inverse)
  wt.shifted@V <- align.coef(wt@V, wavelet=FALSE, coe=coe, modwt=modwt, inverse=inverse)
  if(!inverse) wt.shifted@aligned <- TRUE else wt.shifted@aligned <- FALSE

  return(wt.shifted)
}
