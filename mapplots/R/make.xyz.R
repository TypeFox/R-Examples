make.xyz <-
function(x,y,z,group,FUN=sum,...){
  Z <- tapply(z,list(paste(x,y,sep=', '),group),FUN,...)
  Z <- ifelse(is.na(Z),0,Z)
  XY <- rownames(Z)
  tempfun <- function(XY,i) {
    as.numeric(unlist(lapply(strsplit(XY,', '),function(x) x[i])))
  }
  X <- tempfun(XY,1)
  Y <- tempfun(XY,2)
  return(list(x=X,y=Y,z=Z))
}

