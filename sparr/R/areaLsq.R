areaLsq <- function(xh,WIN,iter=1000){
    if(is.na(xh[3])||is.na(xh[1])||is.na(xh[2])) return(NA)
    return(AREA_BASIC(ndim=2,lower=c(xh[1]-4*xh[3],xh[2]-4*xh[3]),upper=c(xh[1]+4*xh[3],xh[2]+4*xh[3]),functn=Lsq,minpts=iter,uh=xh,WIN=WIN)$val)
}