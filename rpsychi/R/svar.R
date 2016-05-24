svar <- function(x, na.rm=TRUE){			#BeginFnc
  if(na.rm==FALSE){
    x <- na.omit(x)
  }
  output <- 1/length(x)*sum((x-mean(x))^2)
  return(output)
}					#EndFnc

