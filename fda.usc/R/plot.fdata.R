
################################################################################
plot.fdata<-function(x,type,main,xlab,ylab,mfrow=c(1,1),time=1,...) {
if (any(class(x)=="fdata2d"))  {
#stop("Object is not fdata2d class")
if (missing(type)) type="image.contour"
#if (missing(main)) main=x[["names"]][["main"]]
#if (missing(xlab)) xlab=x[["names"]][["xlab"]]
#if (missing(ylab)) ylab=x[["names"]][["ylab"]]
dm<-dim(x$data)
j<-1
len.dm<-length(dm)
rng<-range(x$data)
#if (!ask) {
#   if (dm[1]>9)  ask=TRUE
#   else {
#       ask=FALSE
#       if (dm[1]==1) par(mfrow=c(1,1),ask=FALSE) 
#       if (dm[1]==2) par(mfrow=c(1,2),ask=FALSE)        
#       if (dm[1]==3) par(mfrow=c(1,3),ask=FALSE)        
#       if (dm[1]==4) par(mfrow=c(2,2),ask=FALSE)        
#       if (dm[1]>4) par(mfrow=c(2,3),ask=FALSE)                             
#       if (dm[1]>6) par(mfrow=c(3,3),ask=FALSE)                                    
#   }}
par(mfrow=mfrow)
#npar<-par()$mfrow[1]*par()$mfrow[2]
npar<-mfrow[1]*mfrow[2]
for (i in 1:dm[1]) {
#if (ask) {   par(mfrow=c(1,1),ask=ask)     }
z <- x[["data"]][i,,]
if (len.dm==3) {
switch (type,
"persp"={                          
par(bg = "white")
xx <- x[["argvals"]][[1]]
y <- x[["argvals"]][[2]]


#nrz <- nrow(z);
#ncz <- ncol(z)
#jet.colors <- colorRampPalette( c("yellow", "red") ) 
#nbcol <- length(xx)
#color <- jet.colors(nbcol)
#zfacet <- z[-1, -1] + z[-1, -ncz] + z[-nrz, -1] + z[-nrz, -ncz]
#facetcol <- cut(zfacet, nbcol)
 persp(x=xx,y=y,z=z,xlim=x[["rangeval"]][[1]],ylim=x[["rangeval"]][[2]],zlim=rng,
 main = paste(x$names$main," ",dimnames(x$data)[[1]][i],sep=""),xlab = x$names$xlab[1], ylab =  x$names$xlab[2],...) 
# main =  x$names$main[i],xlab = x$names$xlab[1], ylab =  x$names$xlab[2],...)
#, col=color[facetcol],...)
# par(op)
},
"filled.contour"={
filled.contour(x=x[["argvals"]][[1]],y=x[["argvals"]][[2]],z=z,
xlim=x[["rangeval"]][[1]],ylim=x[["rangeval"]][[2]],
plot.title=title(main = paste(x$names$main," ",dimnames(x$data)[[1]][i],sep=""),
xlab =x$names$xlab[1], ylab =  x$names$xlab[2]),...)},

"contour"={contour(x=x[["argvals"]][[1]],y=x[["argvals"]][[2]],z=z,zlim=rng,
xlim=x[["rangeval"]][[1]],ylim=x[["rangeval"]][[2]],
plot.title=title(main =  paste(x$names$main," ",dimnames(x$data)[[1]][i],sep=""),
    xlab = x$names$xlab[1], ylab =  x$names$xlab[2]),...)},#labels repetidos
    
"image"={image(x = x[["argvals"]][[1]],y = x[["argvals"]][[2]],z=z,
xlim = x[["rangeval"]][[1]],ylim = x[["rangeval"]][[2]],zlim=rng,
main = paste(x$names$main," ",dimnames(x$data)[[1]][i],sep=""),
    xlab = x$names$xlab[1], ylab =  x$names$xlab[2],...)},
"image.contour"={image(x = x[["argvals"]][[1]],y = x[["argvals"]][[2]],z=z,
xlim = x[["rangeval"]][[1]],ylim = x[["rangeval"]][[2]],zlim=rng,
main = paste(x$names$main," ",dimnames(x$data)[[1]][i],sep=""),
    xlab = x$names$xlab[1], ylab =  x$names$xlab[2],...)
    contour(x=x[["argvals"]][[1]],y=x[["argvals"]][[2]],z=z,add = TRUE, drawlabels = FALSE,...)
    }#lattice plot    
#"contourplot"={contourplot(data=z,
#xlim = x[["rangeval"]][[1]],ylim = x[["rangeval"]][[2]],...)}    
)
}


if (names(dev.cur())!="pdf" & j==npar) {
   Sys.sleep(time)
   j<-1
   }
else j<-j+1   
}  
#else {
#  if (len.dm==3) { contourplot(data=z,
#xlim = x[["rangeval"]][[1]],ylim = x[["rangeval"]][[2]],...)}
#    }
if (len.dm>3) stop("Not implemented plot for arrays of more than 3 dimension yet")
}
else {
if (!is.fdata(x))  stop("Object is not fdata class")
if (missing(type)) type="l"
if (missing(main)) main=x[["names"]][["main"]]
if (missing(xlab)) xlab=x[["names"]][["xlab"]]
if (missing(ylab)) ylab=x[["names"]][["ylab"]]
if (is.vector(x[["data"]])) matplot.default(x[["argvals"]],(x[["data"]]),type=type,main=main,ylab=ylab,xlab=xlab,...)
else matplot.default(x[["argvals"]],t(x[["data"]]),type=type,main=main,ylab=ylab,xlab=xlab,...)
}
}

lines.fdata=function(x,...){plot(x,add=TRUE,...)}

title.fdata<-function(x,main=NULL,xlab=NULL,ylab=NULL,rownames=NULL) {
if (!is.fdata(x))  stop("Object is not fdata class")
if (!is.null(rownames)) rownames(x[["data"]])<-rownames
if (!is.null(main)) x[["names"]][["main"]]<-main
if (!is.null(xlab)) x[["names"]][["xlab"]]<-xlab
if (!is.null(ylab)) x[["names"]][["ylab"]]<-ylab
x
}
################################################################################
