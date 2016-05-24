
####################################
#### Plot data
####################################
plotData <- function(Y=NULL, loc=NULL, Yp=NULL, locp=NULL, Yt=NULL, loct=NULL, bdry=NULL, cols=1:3, pchs=1:3, size=c(0.3, 2.7), ...){
  if(!is.null(bdry)){
    if(!is.list(bdry)) bdry <- list(bdry)
    n <- length(bdry)
    tt <- sapply(bdry, function(t) apply(t, 2, range))
    tt1 <- apply(tt, 1, range)
    xrange <- range(tt1[1:4])
    yrange <- range(tt1[5:8])
    plot(xrange, yrange, type="n", ...)
    for(i in 1:n) lines(bdry[[i]])
  } else{
    plot( c(loc[,1],locp[,1],loct[,1]), c(loc[,2],locp[,2],loct[,2]), 
          type="n", ...)
  }
  if(is.null(Y)) Y <- 0
    points(loc[,1], loc[,2], col=cols[1], pch=pchs[1],
          cex = size[1]+size[2]*(Y-min(Y))/(max(Y)-min(Y)) )
  if(!is.null(Yp)){
    points(locp[,1], locp[,2], col=cols[2], pch=pchs[2],
          cex = size[1]+size[2]*(Yp-min(Yp))/(max(Yp)-min(Yp)) )
  }
  if(!is.null(Yt)){
    points(loct[,1], loct[,2], col=cols[3], pch=pchs[3],
          cex = size[1]+size[2]*(Yt-min(Yt))/(max(Yt)-min(Yt)) )
  }
}
####################################
#### Plot Texas
####################################
plotTexas <- function(TexasCounty.boundary, 
                      ind.col = sample(2:5, 254, replace=T)){ 
bdry <- TexasCounty.boundary
n <- length(bdry)
tt <- sapply(bdry, function(t) apply(t, 2, range))
tt1 <- apply(tt, 1, range)
xrange <- range(tt1[1:4])
yrange <- range(tt1[5:8])
plot(xrange, yrange, type="n", xlab = "Longitude", ylab = "Latitude")
for(i in 1:n) polygon(bdry[[i]], col=ind.col[i])
}
####################################
#### END
####################################
