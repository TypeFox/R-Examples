vonNeumann <-
function(x,n)  {
  rx <- NULL
  d <- max(2,length(unlist(strsplit(as.character(x),""))));
  getNext <- function(x,d)  {
    temp <- x^2
    tbs <- as.numeric(unlist(strsplit(as.character(temp),""))) # to be split
    tbs_n <- length(tbs);
    diff_n <- 2*d - tbs_n;
    dn <- ceiling(d/2)
    ifelse(diff_n == 0, tbs <- tbs, tbs <- c(rep(0,diff_n),tbs))
    tbs_n <- length(tbs)
    NEXT <- tbs[-c(1:dn,((tbs_n-dn+1):tbs_n))]
    return(as.numeric(paste(NEXT,collapse="")))
  }
  rx[1] <- x
  for(i in 2:(n+1)) rx[i] <- getNext(rx[i-1],d)  
  return(rx)
}
