geoshape <- function(longvar,latvar,linefile=NULL,pointfile=NULL,coormatrix=NULL) {
  library(spatstat)

  if (!identical(linefile,NULL)) {
    filetype = "line"
    shpfile <- linefile
  }
  if (!identical(pointfile,NULL)) {
    filetype = "point"
    shpfile <- pointfile
  }
  if (!identical(coormatrix,NULL)) {
    filetype = "matrix"
    shpfile <- coormatrix
  }

  amat <- cbind(longvar,latvar)
  amatnew <- unique(amat)
  shpmat <- coordinates(shpfile)
  if (filetype=="point"|filetype=="matrix") {bmat <- shpmat}
  if (filetype=="line") {
    n = length(shpmat)
    btemp <- unlist(shpmat[1])
    n1 = length(btemp)/2
    n2 = 2*n1
    bmat <- cbind(btemp[1:n1],btemp[(n1+1):n2])
    if (n>1) {
      for (i in seq(2,n)) {
        btemp <- unlist(shpmat[i])
        n1 = length(btemp)/2
        n2 = 2*n1
        bmat <- rbind(bmat,cbind(btemp[1:n1],btemp[(n1+1):n2]))
      }
    }
  }
  bmat <- unique(bmat)

  lo1 <- amatnew[,1]
  la1 <- amatnew[,2]
  lo2 <- bmat[,1]
  la2 <- bmat[,2]
  lobar = mean(lo1)
  labar = mean(la1)
  fit <- geodistance(lo1,la1,lobar,labar,dcoor=TRUE)
  dnorth1 <- fit$dnorth
  deast1 <- fit$deast
  fit <- geodistance(lo2,la2,lobar,labar,dcoor=TRUE)
  dnorth2 <- fit$dnorth
  deast2 <- fit$deast
  xlim = c(min(deast1,deast2),max(deast1,deast2))
  ylim = c(min(dnorth1,dnorth2),max(dnorth1,dnorth2))

  amap <- ppp(deast1,dnorth1,xrange=xlim,yrange=ylim)
  bmap <- ppp(deast2,dnorth2,xrange=xlim,yrange=ylim) 

  amatnew <- cbind(amatnew,nncross(amap,bmap)[,1])
  n = nrow(amat)
  amat <- cbind(seq(1:n),amat)
  colnames(amat) <- c("obs","longitude","latitude")
  colnames(amatnew) <- c("longitude","latitude","dist")
  amat <- merge(amat,amatnew,by=c("longitude","latitude"))
  o <- order(amat[,"obs"])
  amat <- amat[o,]
  return(amat[,"dist"])
}

