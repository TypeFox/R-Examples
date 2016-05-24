gLog.ep <-
function(e, zero.adj=0.1, p=2, na.replace=NA)
{
#check:
if(p <= 0) stop("Power must be greater than zero")

if(is.zoo(e)){
  zoo.chk <- TRUE
  e.index <- index(e)
  e <- as.vector(e)
}else{zoo.chk <- FALSE}

na.where <- which(is.na(e))
if(length(na.where) > 0){
  e[na.where] <- 0
}
zero.where <- which(e==0)
eabs <- abs(e)
if(length(zero.where) > 0){
  eabsq <- quantile(eabs[-zero.where], zero.adj)
  eabsadj <- eabs + as.numeric(eabs==0)*eabsq
  eabsp <- eabsadj^p
}else{eabsp <- eabs^p}
logep <- log(eabsp)
if(length(na.where) > 0){
  if(is.na(na.replace)){
    logep[na.where] <- NA
  }else{
    logep[na.where] <- na.replace
  }
}
if(zoo.chk){
  logep <- zoo(logep, order.by=e.index)
}

return(logep)
}
