csMCA <- function(data,subcloud=rep(TRUE,times=nrow(data)),excl=NULL,ncp=5,row.w=rep(1,times=nrow(data))) {
    row.w <- row.w/sum(row.w)*nrow(data)
    row.wc <- row.w[subcloud]
    if(is.null(excl)) excl <- 99999
    N <- nrow(data)
    n <- nrow(data[subcloud,])
    n.w <- sum(row.wc)
    Q <- ncol(data)
    Z <- as.matrix(dichotom(data,out='numeric'))
    fK <- colSums((row.w*Z)[subcloud,])[-excl]/n.w
    FK <- colSums(row.w*Z)[-excl]/N
    K <- ncol(Z)
    Kp <- ncol(Z[,-excl])
    eIp <- matrix(rep(1,length=n),ncol=1)
    eKp <- matrix(rep(1,length=Kp),ncol=1)
    NKc <- diag(colSums((row.w*Z)[subcloud,]))[-excl,-excl]
    Zc <- Z[subcloud,-excl]
    Z0c <- Zc-(1/n.w)*eIp%*%t(eKp)%*%NKc
    Hc <- sqrt(row.wc)*sqrt(N/(n.w*Q))*Z0c%*%diag(1/sqrt(colSums(row.w*Z)[-excl]))
    svd <- svd(Hc)
    YIpc <- (1/sqrt(row.wc))*sqrt(n.w)*svd$u%*%diag(svd$d)
    YKc <- sqrt(N*Q)*diag(1/sqrt(colSums(row.w*Z)[-excl]))%*%svd$v%*%diag(svd$d)
    dims <- paste('dim',1:ncp,sep='.')
    noms <- vector(length=ncol(Z))
    id=0
    for(i in 1:Q) {
      for(j in 1:length(levels(data[,i]))) {
        id=id+1
        noms[id] <- paste(colnames(data)[i],levels(data[,i])[j],sep='.')
      }}
    eig <- list(svd$d*svd$d)
    eig[[2]] <- round(eig[[1]]/sum(eig[[1]])*100,2)
    eig[[3]] <- cumsum(eig[[2]])
    seuil <- 1/Q
    e <- eig[[1]][eig[[1]]>=seuil]
    pseudo <- (Q/(Q-1)*(e-seuil))^2
    eig[[4]] <- pseudo/sum(pseudo)*100
    eig[[5]] <- cumsum(eig[[4]])
    names(eig) <- c('eigen','rate','cum.rate','mrate','cum.mrate')
    weight <- n*fK
    coord <- YIpc[,1:ncp]
    contrib <- 100*row.wc/n.w*coord*coord/matrix(rep(eig[[1]][1:ncp],times=n),ncol=ncp,nrow=n,byrow=T)
    dimnames(coord) <- list(rownames(data)[subcloud],dims) # new
    dimnames(contrib) <- list(rownames(data)[subcloud],dims) # new
    ind <- list(coord=coord,contrib=round(contrib,6))
    coord <- YKc[,1:ncp]
    #Vspe <- sum(fK*(1-fK)/FK)/Q
    contrib <- 100*(FK/Q)*coord*coord/matrix(rep(eig[[1]][1:ncp],times=Kp),ncol=ncp,nrow=Kp,byrow=T)
    s <- vector()
    for(i in 1:Q) s <- c(s,rep(i,times=length(levels(data[,i]))))
    s <- s[-excl]
    v.contrib <- aggregate(contrib,list(s),sum)[,-1]
    dimnames(v.contrib) <- list(colnames(data),dims)
    ctr.cloud <- data.frame(100*(fK*(1-fK)/(Q*FK))/sum(eig[[1]]))
    colnames(ctr.cloud) <- 'ctr.cloud'
    vctr.cloud <- aggregate(ctr.cloud,list(s),FUN=sum)[-1]
    dimnames(vctr.cloud) <- list(colnames(data),'vctr.cloud')
    cos2 <- coord*coord*FK*FK/(fK*(1-fK))
    dimnames(coord) <- list(noms[-excl],dims)
    dimnames(contrib) <- list(noms[-excl],dims)
    dimnames(cos2) <- list(noms[-excl],dims)
    eta2 <- matrix(nrow=Q,ncol=ncp)
    for(j in 1:Q) eta2[j,] <- apply(ind$coord,2,function(x) summary(lm(x~data[subcloud,j],weights=row.wc))$r.squared)
    dimnames(eta2) <- list(colnames(data),dims)
    v.test <- sqrt(cos2)*sqrt(n.w-1)*(((abs(coord)+coord)/coord)-1)
    var <- list(weight=round(weight,1),coord=coord,contrib=round(contrib,6),ctr.cloud=round(ctr.cloud,6),cos2=round(cos2,6),v.test=round(v.test,6),eta2=round(eta2,6),v.contrib=v.contrib,vctr.cloud=vctr.cloud)
    marge.col <- colSums((row.w*Z)[subcloud,])[-excl]/(n.w*Q) # new
    names(marge.col) <- noms[-excl]
    marge.row <- rep(1/(n.w*Q),times=n)
    names(marge.row) <- 1:n
    quali <- 1:Q
    call <- list(X=data,marge.col=marge.col,marge.row=marge.row,ncp=ncp,quali=quali,subcloud=subcloud,excl=excl,row.w=row.w)
    RES <- list(eig=eig,call=call,ind=ind,var=var,svd=list(vs=svd$d,U=svd$u,V=svd$v))
    attr(RES,'class') <- c('csMCA','list')
    RES
}
