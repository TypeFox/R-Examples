cpa<- function(senso, hedo, coord=c(1,2),center=TRUE,scale=TRUE,nb.clusters=0,scale.unit=FALSE,col=terrain.colors(45)[1:41]) {

colplot<-function(mat, k=0,coord, z, level=41, col = terrain.colors(level+level%/%10)[1:level], xlab="", ylab="") { #heat.colors(level)

    abs <- coord[1]
    ord <- coord[2]
    x <- mat[,abs]
    y <- mat[,ord]
    z <- mat[,z]
#    x1 <- min(z)
#    x2 <- max(z)
    x1 <- -1
    x2 <- 1
    plot(mat[,abs],mat[,ord],xlab=xlab, ylab=ylab,asp=1,type="n")
    legend("topleft",legend=c(1,0.75,0.5,0.25,0,-0.25,-0.5,-0.75,-1),fill=c(col[level],col[(level%/%2)+(level%/%4)+(level%/%8)+1],col[(level%/%2)+(level%/%4)+1],col[(level%/%2)+(level%/%8)+1],col[(level%/%2)+1],col[(level%/%4)+(level%/%8)+1],col[(level%/%4)+1],col[(level%/%8)+1],col[1]),cex=0.7)
    abline(v=0,lty=2)
    abline(h=0,lty=2)
####rect(0, levels[-length(levels)], 1, levels[-1], col = col)
    n <- nrow(mat)
    h <- (x2-x1)/level

    for (ind in 1:(n-k)) points(x[ind],y[ind],col=col[max(1,(z[ind]-x1)%/%h)],pch=20)
    for (ind in (n-k+1):n) points(x[ind],y[ind],col=col[max(1,(z[ind]-x1)%/%h)],pch=15,cex=1)
    for (ind in (n-k+1):n) text(x[ind],y[ind],col=col[max(1,(z[ind]-x1)%/%h)],rownames(mat)[ind],cex=1,pos = 1, offset = 0.05)
}

### Main program

    if (max(coord) > (nrow(hedo)-1)) {
      print (paste("Problem with coord. Max (coord) must be less than",nrow(hedo)-1," Axes 1-2 will be taken",sep=""))
      coord=c(1,2)
    }
    senso <- scale(senso,center=center,scale=scale)[,]
    hedo <- scale(hedo,center=center,scale=scale)[,]
    if (scale) senso=senso*sqrt(nrow(senso)/(nrow(senso)-1))
    if (scale) hedo=hedo*sqrt(nrow(hedo)/(nrow(hedo)-1))
    op <- par(no.readonly = TRUE)
    on.exit(par(op))
    senso <- data.frame(senso)
    hedo <- data.frame(hedo)
    nbjuge <- ncol(hedo)
    nbdesc <- ncol(senso)

  classif <- cluster::agnes(dist(t(hedo)),method="ward")
  plot(as.dendrogram(classif),main="Cluster Dendrogram",xlab="Panelists") 
    if (nb.clusters==0){
       classif2 <- as.hclust(classif)
       nb.clusters = which.max(rev(diff(classif2$height))) + 1
#      classif=hopach(t(MatH),d="euclid",K=10,mss="mean")
#      nb.clusters=classif$clustering$k
    }
  clusters=kmeans(t(hedo),centers=nb.clusters)$cluster
  mat <- matrix(0,nb.clusters,nrow(hedo))
  dimnames(mat) <- list(1:nb.clusters,rownames(hedo))
  for (i in 1:nb.clusters){
    mat[i,] <- apply(t(hedo[,clusters==i]),2,mean)
    rownames(mat)[i] <- paste("cluster",i)
  }  
  desc.clusters=cor(senso,t(mat),use="pairwise.complete.obs")

    A <- rbind.data.frame(t(hedo),mat,t(senso))
    colnames(A) <- row.names(hedo)
    result <- A
    auxil = cbind.data.frame(A,as.factor(c(clusters,rep(1,nrow(mat)+ncol(senso)))))
    colnames(auxil)[ncol(A)+1]="cluster"
    hedo.pca <- PCA(auxil,quali.sup=ncol(A)+1,ind.sup=(nbjuge+1):nrow(A),scale.unit=scale.unit,graph=FALSE,ncp = min(nbjuge-1,ncol(A)))
    plot(hedo.pca,choix="ind",axes=coord,cex=0.7,habillage=ncol(A)+1)
    plot(hedo.pca,choix="var",axes=coord)
    TA <- t(A)
    coef <- matrix(NA,nbjuge+nb.clusters,nbdesc)
    for (d in 1:nbdesc) {
      coef[1:nbjuge,d] <- cor(TA[,1:nbjuge],TA[,nbjuge+nb.clusters+d],use="pairwise.complete.obs")
      coef[(nbjuge+1):(nbjuge+nb.clusters),d] <- cor(TA[,(nbjuge+1):(nbjuge+nb.clusters)],TA[,nbjuge+nb.clusters+d],use="pairwise.complete.obs")
    }
    coef <- data.frame(coef)
    colnames(coef) <- colnames(senso)
    B <- cbind.data.frame(rbind.data.frame(hedo.pca$ind$coord,hedo.pca$ind.sup$coord[1:nb.clusters,]),coef)
    for (d in 1:nbdesc) {
dev.new()
      par(mar = c(4.2,4.1,3.5,2))
      colplot(as.matrix(B), k=nb.clusters,coord, (nrow(hedo)+d),col=col, xlab=paste("Dim",coord[1]," (",signif(hedo.pca$eig[coord[1],2],4),"%)",sep=""), ylab=paste("Dim",coord[2]," (",signif(hedo.pca$eig[coord[2],2],4),"%)",sep=""))
      points(hedo.pca$ind.sup$coord[nb.clusters+d,coord[1]],hedo.pca$ind.sup$coord[nb.clusters+d,coord[2]],col="red",pch=15,cex=0.8)
      text(hedo.pca$ind.sup$coord[nb.clusters+d,coord[1]],hedo.pca$ind.sup$coord[nb.clusters+d,coord[2]],col="red",labels=colnames(B)[nrow(hedo)+d],pos = 1, offset = 0.05)
      title(main = paste("Consumers' preferences analysed by",colnames(B)[nrow(hedo)+d]),cex.main = 1.1, font.main = 2)
    }
    don <- cbind.data.frame(as.factor(clusters),t(hedo))
    colnames(don) <- c("clusters",paste("Prod",rownames(hedo),sep="."))
    resdecat <- decat(don,formul="~clusters",firstvar=2,proba=1,graph=FALSE)
    res <- list()
    res$clusters <- clusters
    res$result <- result
    res$prod.clusters <- resdecat$resT
    res$desc.clusters <- desc.clusters
    return(res)
}

