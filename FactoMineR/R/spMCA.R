spMCA <- function(X,excl=NULL,ncp=5,quali.sup=NULL,graph=TRUE) {
    moy.p <- function(V, poids) {
        res <- sum(V * poids)/sum(poids)
    }
    X <- as.data.frame(X)
	X <- droplevels(X)
	if (is.null(excl)) {
	  warning("Function MCA has been used instead of spMCA. You should use the function MCA")
	  return(res=MCA(X,quali.sup=quali.sup,ncp=ncp))
	}
    n <- nrow(as.data.frame(X))
    Q <- ncol(as.data.frame(X))
	if (!is.null(quali.sup)){
      nmod <- c(0,cumsum(sapply(X,nlevels)))
	  deb <- nmod[quali.sup]+1
	  fin <- nmod[quali.sup+1]
      for (i in 1:length(quali.sup)) excl=c(excl,deb[i]:fin[i])
	}
	Qact=Q-length(quali.sup)
    Z <- tab.disjonctif(X)
    K <- ncol(Z)
    Z0 <- Z-(1/n)*matrix(1,ncol=n,nrow=n)%*%Z
	H0 <- (1/sqrt(Qact))*sweep(Z0,2,sqrt(colSums(Z)),FUN="/")
	poids <- rep(1,ncol(Z0))
	poids[excl]=1e-15
    svd2 <- svd.triplet(H0,col.w=poids)
	poids[excl]=0
    dims <- paste('dim',1:ncp,sep='.')
    noms <- NULL
    for (j in 1:ncol(X)) noms = c(noms, levels(X[, j]))
    for (j in 1:ncol(X)) {
      if (sum(noms %in% levels(X[, j])) != nlevels(X[, j])) levels(X[, j]) = paste(colnames(X)[j], levels(X[, j]), sep = "_")
    }
    YIt <- sweep(svd2$U,2,svd2$vs[1:ncol(svd2$U)],FUN="*")*sqrt(n)
    YKpt <- n*sqrt(Qact)*sweep(sweep(svd2$V,2,svd2$vs[1:ncol(svd2$V)],FUN="*"),1,sqrt(colSums(Z)),FUN="/")
    vp <- svd2$vs^2*n
    eig <- as.data.frame(matrix(NA, length(vp), 3))
    rownames(eig) <- paste("comp", 1:length(vp))
    colnames(eig) <- c("eigenvalue", "percentage of variance", "cumulative percentage of variance")
    eig[, "eigenvalue"] <- vp
    eig[, "percentage of variance"] <- (vp/sum(vp)) * 100
    eig[, "cumulative percentage of variance"] <- cumsum(eig[,"percentage of variance"])
    coord <- YIt[,1:ncp]
    contrib <- 100/n*sweep(as.matrix(coord^2), 2, vp[1:ncol(coord)], FUN = "/")
    dimnames(coord) <- dimnames(contrib) <- list(1:n,dims)
    ind <- list(coord=coord,contrib=round(contrib,6))
    coord <- YKpt[,1:ncp]
    fK <- colSums(Z)/n
    contrib <- 100*(fK/Qact)*sweep(sweep(coord^2,2,vp[1:ncol(coord)],FUN="/"),1,poids,FUN="*")
    s <- vector()
    for(i in 1:Q) s <- c(s,rep(i,times=nlevels(X[,i])))
    v.contrib <- aggregate(contrib,list(s),sum)[,-1]
    dimnames(v.contrib) <- list(colnames(X)[unique(s)],dims)
    ctr.cloud <- data.frame(100*(1-fK)/(ncol(Z)-length(excl)-Qact))
    colnames(ctr.cloud) <- 'ctr.cloud'
    vctr.cloud <- aggregate(ctr.cloud,list(s),FUN=sum)[-1]
    dimnames(vctr.cloud) <- list(colnames(X)[unique(s)],'vctr.cloud')
    cos2 <- coord*coord/((1/fK)-1)
    dimnames(coord) <- dimnames(contrib) <- dimnames(cos2) <- list(noms,dims)
    eta2 <- matrix(nrow=Q,ncol=ncp)
    for(j in 1:Q) {
      vrc <- aggregate(ind$coord,list(X[,j]),var)[,-1]
      wi <- apply(vrc,2,moy.p,poids=as.numeric(table(X[,j])))
      be <- eig[[1]][1:ncp]-wi
      eta2[j,] <- be/eig[[1]][1:ncp]
      }
    dimnames(eta2) <- list(colnames(X),dims)
    v.test <- sqrt(cos2)*sqrt(n-1)*sign(coord)
    if (!is.null(quali.sup)){
	  var <- list(coord=coord[-excl,,drop=FALSE],contrib=round(contrib[-excl,,drop=FALSE],6),cos2=round(cos2[-excl,,drop=FALSE],6),v.test=round(v.test[-excl,,drop=FALSE],6),eta2=round(eta2[-quali.sup,,drop=FALSE],6),v.contrib=round(v.contrib[-quali.sup,,drop=FALSE],6))
      quali.sup <- list(coord=coord[excl,,drop=FALSE],cos2=round(cos2[excl,,drop=FALSE],6),v.test=round(v.test[excl,,drop=FALSE],6),eta2=round(eta2[quali.sup,,drop=FALSE],6))
    } else {
	  var <- list(coord=coord[-excl,,drop=FALSE],contrib=round(contrib[-excl,,drop=FALSE],6),cos2=round(cos2[-excl,,drop=FALSE],6),v.test=round(v.test[-excl,,drop=FALSE],6),eta2=round(eta2,6),v.contrib=v.contrib)
      quali.sup <- list(coord=coord[excl,,drop=FALSE],cos2=round(cos2[excl,,drop=FALSE],6),v.test=round(v.test[excl,,drop=FALSE],6))
	}
	marge.col <- colSums(Z[,-excl])/(n*Qact)
    names(marge.col) <- noms[-excl]
    marge.row <- rep(1/(n*Qact),times=n)
    names(marge.row) <- 1:n
    quali <- 1:Q
    call <- list(X=X,mar.col=marge.col,marge.row=marge.row,ncp=ncp,quali=quali,excl=excl)
    RES <- list(eig=eig,call=call,ind=ind,var=var,quali.sup=quali.sup,svd=list(vs=svd2$vs*sqrt(n),U=svd2$U/sqrt(n),V=svd2$V))
    attr(RES,'class') <- c('spMCA','list')
	if (graph){
	  plot(RES)
	  plot(RES,invisible="ind",new.plot=TRUE)
	  plot(RES,choix="var",new.plot=TRUE)
	}
    return(RES)
}

