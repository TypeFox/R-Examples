ssd <- function(x, na.rm=TRUE){			#BeginFnc
  if(na.rm==TRUE){
    x <- na.omit(x)
  }
  output <- sqrt(1/length(x)*sum((x-mean(x))^2))
  return(output)
}					#EndFnc