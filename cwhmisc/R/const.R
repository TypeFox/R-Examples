c3Q <- .Machine$double.xmax^0.75
c38 <- sqrt(c3Q) 

GreatestIntAsRealF <- function( ) {
  ## Find the greatest integer K which is distiguishable from (K+1), both represented as real
##    p <- 0.0
##    q <- 1.0
##    while (p != q) { q <- q+q;  p <- q+1.0 };
##    return( q/2.0)
##  
  1.0/ .Machine$double.eps
}  ## END  GreatestIntAsRealF;

ASCII <- c(NA, sapply(1:255, function(i) parse(text=paste("\"\\", structure(i,class="octmode"), "\"", sep=""))[[1]]));

HexDig <- c('0','1','2','3','4','5','6','7','8','9',LETTERS[1:6])

TAU <- (1+sqrt(5))/2
