url2data <- function(url.){
  fili <- names(url.)
  dati <- sapply(strsplit(fili, '\\.'), function(x)x[1:(length(x)-1)])
  n <- length(fili)
  found <- rep(NA, n)
  if(n>0)for(i in 1:n){
    chk <- try(download.file(url.[i], fili[i]))
    found[i] <- (c("TRUE", "FALSE")[1 +
               (class(chk)=="try-error")])
  }
  cbind(data=dati, file=fili, url=url., found=found) 
}
