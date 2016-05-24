`tapply.stat` <-
function (y, x, stat = "mean")
{
k<-0
numerico<- NULL
if(is.null(ncol(x))){
if(is.numeric(x)){
k<-1
numerico[1]<-1
}
}
else {
ncolx<-ncol(x)
for (i in 1:ncolx) {
if(is.numeric(x[,i])){
k<-k+1
numerico[k]<-i
}}}
    cx <- deparse(substitute(x))
    cy <- deparse(substitute(y))
    x <- data.frame(c1 = 1, x)
    y <- data.frame(v1 = 1, y)
    nx <- ncol(x)
    ny <- ncol(y)
    namex <- names(x)
    namey <- names(y)
    if (nx == 2)
        namex <- c("c1", cx)
    if (ny == 2)
        namey <- c("v1", cy)
    namexy <- c(namex, namey)
    for (i in 1:nx) {
        x[, i] <- as.character(x[, i])
    }
    z <- NULL
    for (i in 1:nx) {
        z <- paste(z, x[, i], sep = "&")
    }
    w <- NULL
    for (i in 1:ny) {
        m <- tapply(y[, i], z, stat)
        m <- as.matrix(m)
        w <- cbind(w, m)
    }
    nw <- nrow(w)
    c <- rownames(w)
    v <- rep("", nw * nx)
    dim(v) <- c(nw, nx)
    for (i in 1:nw) {
        for (j in 1:nx) {
            v[i, j] <- strsplit(c[i], "&")[[1]][j + 1]
        }
    }
    rownames(w) <- NULL
    junto <- data.frame(v[, -1], w)
    junto <- junto[, -nx]
    names(junto) <- namexy[c(-1, -(nx + 1))]
    if(k==1 & nx==2) {
    junto[,numerico[1]]<-as.character(junto[,numerico[1]])
    junto[,numerico[1]]<-as.numeric(junto[,numerico[1]])    
    junto<-junto[order(junto[,1]),]
    }
    if (k>0 & nx > 2) {
    for (i in 1:k){
    junto[,numerico[i]]<-as.character(junto[,numerico[i]])
    junto[,numerico[i]]<-as.numeric(junto[,numerico[i]])
    }
    junto<-junto[do.call("order", c(junto[,1:(nx-1)])),]
    }
    rownames(junto)<-1:(nrow(junto))
    return(junto)
}

