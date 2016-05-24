plot.IdMap <- function(x,xlim,ylim,levels.contour=NULL,color=FALSE,inverse=FALSE, ...){
    res <- x
    if (!inherits(res,"IdMap"))
        stop("Non convenient data")
    if (!length(xlim)==2)
        stop("'xlim' should contain the minimum and maximum coordinates on the x axis. Two coordinates should be indicated")
    if (!length(ylim)==2)
        stop("'ylim' should contain the minimum and maximum coordinates on the y axis. Two coordinates should be indicated")
    x1=min(xlim)
    x2=max(xlim)
    y1=min(ylim)
    y2=max(ylim)
    if (!is.numeric(x1) || !is.numeric(x2) || !is.numeric(y1) || !is.numeric(y2))
        stop("'x1','x2','y1' and/or 'y2' should be numerical values.")
    pos.x1 <- grep(paste("X_",x1,sep=""),colnames(res$idmap$data))
    if (length(pos.x1)>1){
        pos.temp <- NULL
        for (l in 1:length(pos.x1))
            if (paste("X_",x1,sep="")==colnames(res$idmap$data)[pos.x1[l]])
                pos.temp <- pos.x1[l]
        pos.x1 <- pos.temp
    }
    if (length(pos.x1)==0 || is.null(pos.x1))
        stop(paste("Not convenient 'x1' definition. Moreover, 'x1' should be a multiplicator of ",res$idmap$precision,".",sep=""))
    pos.x2 <- grep(paste("X_",x2,sep=""),colnames(res$idmap$data))
    if (length(pos.x2)>1){
        pos.temp <- NULL
        for (l in 1:length(pos.x2))
            if (paste("X_",x2,sep="")==colnames(res$idmap$data)[pos.x2[l]])
                pos.temp <- pos.x2[l]
        pos.x2 <- pos.temp
    }
    if (length(pos.x2)==0 || is.null(pos.x2))
        stop(paste("Not convenient 'x2' definition. Moreover, 'x2' should be a multiplicator of ",res$idmap$precision,".",sep=""))
    pos.y1 <- grep(paste("Y_",y1,sep=""),rownames(res$idmap$data))
    if (length(pos.y1)>1){
        pos.temp <- NULL
        for (l in 1:length(pos.y1))
            if (paste("Y_",y1,sep="")==rownames(res$idmap$data)[pos.y1[l]])
                pos.temp <- pos.y1[l]
        pos.y1 <- pos.temp
    }
    if (length(pos.y1)==0 || is.null(pos.y1))
        stop(paste("Not convenient 'y1' definition. Moreover, 'y1' should be a multiplicator of ",res$idmap$precision,".",sep=""))
    pos.y2 <- grep(paste("Y_",y2,sep=""),rownames(res$idmap$data))
    if (length(pos.y2)>1){
        pos.temp <- NULL
        for (l in 1:length(pos.y2))
            if (paste("Y_",y2,sep="")==rownames(res$idmap$data)[pos.y2[l]])
                pos.temp <- pos.y2[l]
        pos.y2 <- pos.temp
    }
    if (length(pos.y2)==0 || is.null(pos.y2))
        stop(paste("Not convenient 'y2' definition. Moreover, 'y2' should be a multiplicator of ",res$idmap$precision,".",sep=""))
    dev.new()
    res2 <- res$idmap$data[pos.y1:pos.y2,pos.x1:pos.x2]
    f1 <- seq(x1,x2,res$idmap$precision)
    f2 <- seq(y1,y2,res$idmap$precision)
    if (!is.null(levels.contour)){
        if (min(levels.contour)<0 || max(levels.contour)>100 || length(levels.contour)<2){
            warning("not convenient 'levels.contour' definition: the default value will be used")
            levels.contour=NULL
        } else {
            oo <- order(levels.contour)
            levels.contour <- levels.contour[oo]
        }
    }
    if (is.null(levels.contour))
        levels.contour <- c(10,15,20,25,30,35,40,45,50)
    if (!sd(as.vector(res$idmap$j.weight))==0){
        titre="Weighted Ideal Mapping"
    } else {
        titre="Ideal Mapping"
    }
    coord <- res$PCA$dim
    if (color){
        image(f1,f2,t(res2),col=terrain.colors(200),xlab=paste("Dim",coord[1]," (",round(res$PCA$eig[coord[1],2],2),"%)",sep=""),ylab=paste("Dim",coord[2]," (",round(res$PCA$eig[coord[2],2],2),"%)",sep=""),main=titre)
        contour(f1,f2,t(res2),nlevels=length(levels.contour),levels=levels.contour,add=T,labex=0)
        for (i in 1:nrow(res$PCA$ind$coord)) {
            points(res$PCA$ind$coord[i,res$PCA$dim[1]],res$PCA$ind$coord[i,res$PCA$dim[2]],pch=15,cex=0.7)
            text(res$PCA$ind$coord[i,res$PCA$dim[1]],res$PCA$ind$coord[i,res$PCA$dim[2]],rownames(res$PCA$ind$coord)[i],pos=4,offset=0.2,cex=0.7)
        }
        abline(v=0,lty=2)
        abline(h=0,lty=2)
    } else {
        if (inverse){
            image(f1,f2,t(res2),col=grey(100:(100-max(res2))/100),xlab=paste("Dim",coord[1]," (",round(res$PCA$eig[coord[1],2],2),"%)",sep=""),ylab=paste("Dim",coord[2]," (",round(res$PCA$eig[coord[2],2],2),"%)",sep=""),main=titre)
            contour(f1,f2,t(res2),nlevels=length(levels.contour),levels=levels.contour,add=T,labex=0,col="black")
            for (i in 1:nrow(res$PCA$ind$coord)) {
                points(res$PCA$ind$coord[i,res$PCA$dim[1]],res$PCA$ind$coord[i,res$PCA$dim[2]],pch=15,cex=1,col="black")
                text(res$PCA$ind$coord[i,res$PCA$dim[1]],res$PCA$ind$coord[i,res$PCA$dim[2]],rownames(res$PCA$ind$coord)[i],pos=4,offset=0.4,cex=1.1,col="black")
            }
            axis(2)
            axis(3,labels=F)
            abline(v=0,lty=2,col="black")
            abline(h=0,lty=2,col="black")
        } else {
            image(f1,f2,t(res2),col=grey(0:max(res2)/100),xlab=paste("Dim",coord[1]," (",round(res$PCA$eig[coord[1],2],2),"%)",sep=""),ylab=paste("Dim",coord[2]," (",round(res$PCA$eig[coord[2],2],2),"%)",sep=""),main=titre)
            contour(f1,f2,t(res2),nlevels=length(levels.contour),levels=levels.contour,add=T,labex=0,col="white")
            for (i in 1:nrow(res$PCA$ind$coord)) {
                points(res$PCA$ind$coord[i,res$PCA$dim[1]],res$PCA$ind$coord[i,res$PCA$dim[2]],pch=15,cex=0.7,col="white")
                text(res$PCA$ind$coord[i,res$PCA$dim[1]],res$PCA$ind$coord[i,res$PCA$dim[2]],rownames(res$PCA$ind$coord)[i],pos=4,offset=0.2,cex=0.7,col="white")
            }
            abline(v=0,lty=2,col="white")
            abline(h=0,lty=2,col="white")
        }
    }
}