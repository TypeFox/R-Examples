kernelFun <- function(type=c('truncated','bartlett','daniell','QS','parzen'),z){ 
 type <- match.arg(type)
 if (missing(type)) stop("Kernel type is missing. Type must be one of 'truncated','bartlett','daniell','QS','parzen'")
 if (type != "truncated" && type != "bartlett" && type != "daniell" && type != "QS" && type != "parzen") stop("Wrong type. Type must be one of 'truncated','bartlett','daniell','QS','parzen'")
 if (length(z)>1) stop("z must be of length 1")
 if (!is.numeric(z)) 
      stop("'z' must be numeric") 
 if (type=="truncated"){
  k <- ifelse(abs(z)<=1,1,0)
 }
 else if (type=="bartlett"){
  k <- ifelse(abs(z)<=1,1-abs(z),0)
 }
 else if (type=="daniell"){
  if (z==0) stop("z cannot be zero")
  k <- sin(pi*z)/(pi*z)
 }
 else if (type=="QS"){
  k <- (9/(5*pi^2*z^2))*(sin(sqrt(5/3)*pi*z)/(sqrt(5/3)*pi*z) - cos(sqrt(5/3)*pi*z))
 }
 else if (type=="parzen"){
  if (abs(z) <= (3/pi)){
   k <- 1-6*((pi*z)/6)^2 + 6*abs((pi*z)/6)^3
  }
  else if ((abs(z) >= (3/pi)) && (abs(z) <= (6/pi))){
   k <- 2*(1-abs((pi*z)/6))^3
  }
  else {
   k <- 0
  }
 }
 return(k)
}
