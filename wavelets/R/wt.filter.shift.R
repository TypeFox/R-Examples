wt.filter.shift <- function(filter, J, wavelet=TRUE, coe=FALSE, modwt=FALSE){

  # error checking
  if(is.na(match(class(filter), c("character", "wt.filter", "numeric", "integer"))))
    stop("Invalid argument: 'filter' must be of class 'character', 'wt.filter', 'numeric' or 'integer'")
  if((length(setdiff(J, round(J)) != 0) | (sum(J <= 0) != 0)))
    stop("Invalid argument: all elements of 'J' must be a positive integers.")

  # create filter if necessary
  if(!is.na(match(class(filter), c("numeric","integer","character"))))
    filter <- wt.filter(filter, modwt=modwt)
  L <- filter@L

  # calculate level 1 shift (equation 112e)
  if(!coe){
    if((filter@wt.class == "Daubechies") & !is.na(match(L, c(2,4)))){
      if(L == 2) nu <- 0
      if(L == 4) nu <- 1
    } else if(filter@wt.class == "Least Asymmetric"){
      if(!is.na(match(L,c(8,12,16,20)))) delta <- 1
      if(L == 10 | L == 18) delta <- 0
      if(L == 14) delta <- 2
      nu <- abs(-(L/2) + delta)
    } else if(filter@wt.class == "Coiflet") {
      nu <- abs(-2*L/3 + 1)
    } else if(filter@wt.class == "Best Localized") {
      if(L == 14) nu <- 5
      if(L == 18) nu <- 11
      if(L == 20) nu <- 9
    }
    else nu <- sum((1:(L-1))*filter@g[-1]^2)/sum(filter@g^2)
  } else {
    nu <- sum((1:(L-1))*filter@g[-1]^2)/sum(filter@g^2)
  }

  # calculate shift for each j > 1 (equation 114c)
  if(filter@wt.name != "haar"){
    shift <- sapply(J, function(j, wavelet){
      if(j > 1){
        if(wavelet) shift <- (2^(j-1))*(L-1) - nu
        else shift <- (2^j -1)*nu
      } else {
        if(wavelet) shift <- L - nu - 1 else shift <- nu
      }
    }, wavelet=wavelet)
  } else {
    if(!modwt){
      shift <- 0^J
    } else {
      shift <- 2^(J-1)
    }
  }

  # calculate shift for dwt
  if(filter@transform == "dwt") shift <- ceiling(((shift+1)/(2^J))-1)
     
  return(shift)
}
