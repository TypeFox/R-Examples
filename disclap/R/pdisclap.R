pdisclap <- 
function(x,p,lower.tail=TRUE){
  if (any(p < 0) | any(p >= 1))
    stop("0 <= p < 1 is required")
  if (length(p) != length(x)) {
    if (length(p) == 1) {
      p <- rep(p, length(x))
    }
    else {
      stop("length(p) != 1 and length(p) != length(x)")
    }
  }
  oneP <- function(x,p,lower.tail){
    GTabsx <- function(x,p) p^(abs(x)+1)/(1+p)
    x <- floor(x)
    if(x>=0 & !lower.tail) return(GTabsx(x,p))
    if(x<0 & !lower.tail) return(1-GTabsx(-x-1,p))
    if(x>=0 & lower.tail) return(1-GTabsx(x,p))
    if(x<0 & lower.tail) return(GTabsx(-x-1,p))
  }
  mapply(oneP,x=x,p=p,lower.tail=lower.tail)
}

