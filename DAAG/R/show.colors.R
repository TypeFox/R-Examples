"show.colors" <-
function(type=c("singles", "shades", "gray"), order.cols=TRUE){
    if(order.cols){
        MASS.out <- try(requireNamespace("MASS", quietly=TRUE), silent = TRUE)
    MASSCheck <- !is.logical(MASS.out) | (MASS.out == FALSE)
    if(MASSCheck)print("Error: package MASS is not installed properly")
    if(MASSCheck) return()
    }
type <- type[1]
oldpar <- par(mar=c(.75, .75,1.5, .75))
on.exit(par(oldpar))
order.cols <- order.cols & requireNamespace("stats")
unique.colors <- function(){
    colnam <- colors()
    vector.code <- apply(col2rgb(colnam),2,function(x)x[1]+x[2]*1000+x[3]*10000)
    unique.code <- unique(vector.code)
    sub <- match(unique.code, vector.code)
    colnam[sub]
}
plotshades <- function(x=1, start=1, nlines=14, numlabels=FALSE, colmat, colnam){
    endlines <- min(start+nlines-1, length(colnam))
    colrange <- start:endlines
    nlines <- length(colrange)
    points(rep(x:(x+4), rep(nlines,5)), nlines+1.25-rep(1:nlines, 5),
        col=as.vector(colmat[colrange,]), pch=15, cex=2.95)
    text(rep(x-0.25,nlines), nlines+1.25-(1:nlines), colnam[colrange],
        adj=0, col=paste(colnam[1:10],"4",sep=""), cex=0.8, xpd=TRUE)
    text((x+1):(x+4), rep(nlines+0.95,4), 1:4, cex=0.75)
}
plotcols <- function(x=1, start=1, wid=5, nlines, numlabels=FALSE, colvec=loners){
    nlines <- min(nlines, length(colvec)-start+1)
    colrange <- start:(start+nlines-1)
    xleft <- rep(x, nlines)
    xright <- xleft+wid
    ybottom <- nlines+1-(1:nlines)
    ytop <- ybottom+1
    rect(xleft, ybottom, xright, ytop, col=colvec[colrange], xpd=TRUE)
    colvals <- lapply(colvec[colrange], function(x){z<-col2rgb(x)/256; 0.4*(1-(1-z)^2)+0.6*(1-z)^2})
    colvals <- sapply(colvals, function(x)rgb(x[1],x[2],x[3]))
    text(rep(x+0.25, nlines), nlines-(1:nlines)+1.5, colvec[colrange],
        col=colvals,  adj=0, cex=0.8, xpd=TRUE)
}
classify.colors <- function(colr, colset=loners){
    gsub <- grep("green",colr)
    rsub <- grep("red",colr)
    bsub <- grep("blue",colr)
    colxyz <- t(col2rgb(colr[c(rsub,gsub,bsub)]))
    colxyz <- data.frame(colxyz, rep(c("red","green","blue"), c(length(rsub),length(gsub),length(bsub))))
    names(colxyz)<- c("red","green","blue","gp")
    col.lda <- MASS::lda(gp ~ red+green+blue, data=colxyz)
    colrgb <- data.frame(t(col2rgb(colset)))
    names(colrgb) <- c("red", "green", "blue")
    newcol <- predict(col.lda, newdata=colrgb)
    newcol
}
allcols <- unique.colors()
gray <- as.logical(match(substring(allcols,1,4), "gray", nomatch=0))
grayshades <- allcols[gray]
nongray <- allcols[!gray]
nlast <- nchar(nongray)
five <- substring(nongray,nlast,nlast) %in% c("1","2","3","4")
fivers <- unique(substring(nongray[five],1,nlast[five]-1))
fiveshades <- outer(fivers,c("","1","2","3","4"),
    function(x,y)paste(x,y,sep=""))
subs <- match(nongray, fiveshades, nomatch=0)
loners <- nongray[subs==0]
print(c(length(loners),length(fiveshades)))
ncolm <- switch(type, singles=3, shades=4, gray=4)
nlines <- switch(type, singles=ceiling(length(loners)/3),
    shades=ceiling(length(fivers)/4), gray=ceiling(length(grayshades)/4))


plot(c(1,21.5), c(1,nlines+1), type="n", axes=FALSE, xlab="", ylab="")
heading <- switch(type, singles="Colors that do not have shades",
 shades="Colors that have 4 or 5 shades", gray="Shades of gray")
mtext(side=3, line=-0.25, heading, at=1, adj=0)

# arrange <- function(colvec){
#    xyz <- t(sweep(col2rgb(colvec),1,c(.2126, .7152, .0722),"*"))
#    red <- xyz[,1]
#    green <- xyz[,2]
#    blue <- xyz[,3]
#    scores <- (red+blue+400)*(green>165)+ (red+green+200)*(red>25)*(green<165)
#        +(green+blue)*(red<25)*(green>165)
#    ord <- order(scores)
#    ord}
arrange <- function(colvec){
newcols <- classify.colors(colr=c(loners,fiveshades), colset=colvec)
n1 <- 1:length(colvec)
blue <- n1[newcols$class=="blue"]
green <- n1[newcols$class=="green"]
red <- n1[newcols$class=="red"]
colblue <- colvec[blue]
colred <- colvec[red]
colgreen <- colvec[green]
ordblue <- order(apply(sweep(col2rgb(colblue),1,c(.2126, .7152, .0722),"*"),2,sum))
ordred <- order(apply(sweep(col2rgb(colred),1,c(.2126, .7152, .0722),"*"),2,sum))
ordgreen <- order(apply(sweep(col2rgb(colgreen),1,c(.2126, .7152, .0722),"*"),2,sum))
c(red[ordred], green[ordgreen], blue[ordblue])
}

if(order.cols){
z <- arrange(colvec=loners)
loners <- loners[z]
z <- arrange(colvec=fiveshades[,3])
fivers <- fivers[z]
fiveshades <- fiveshades[z, ]
}

if(type=="singles"){
plotcols(nlines=nlines, wid=6.5)
plotcols(x=8, nlines=nlines, wid=6.5, start=nlines+1)
plotcols(x=15, nlines=nlines, wid=6.5, start=2*nlines+1, numlabels=TRUE)
}
if(type=="gray"){
plotcols(colvec=grayshades, wid=5, nlines=nlines)
plotcols(x=6.25, colvec=grayshades, wid=5, nlines=nlines, start=nlines+1)
plotcols(x=11.5, colvec=grayshades, wid=5, nlines=nlines, start=2*nlines+1, numlabels=TRUE)
plotcols(x=16.75, colvec=grayshades, wid=5, nlines=nlines, start=3*nlines+1, numlabels=TRUE)
}
if(type=="shades"){
plotshades(nlines=nlines, colmat=fiveshades, colnam=fivers)
plotshades(x=6.5, start=nlines+1, nlines=nlines,colmat=fiveshades,colnam=fivers)
plotshades(x=12, start=2*nlines+1, nlines=nlines, numlabels=TRUE, colmat=fiveshades, colnam=fivers)
plotshades(x=17.5, start=3*nlines+1, nlines=nlines, numlabels=TRUE,colmat=fiveshades, colnam=fivers)
}
invisible(list(singles=loners, shades=fiveshades, grayshades=grayshades))
}
