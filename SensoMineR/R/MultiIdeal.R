MultiIdeal <- function(dataset,col.p,col.j,id.recogn,level.search.desc=0.2,correct=FALSE,nbchoix=NULL,nbsimul=500,coord=c(1,2)){

################################################################################
    hotelling2 <- function(d1, d2, n1 = nrow(d1), n2 = nrow(d2)) {
        k <- ncol(d1)
        xbar1 <- apply(d1, 2, mean)
        xbar2 <- apply(d2, 2, mean)
        dbar <- xbar2 - xbar1
        if (n1 + n2 < 3)
            return(NA)
        v <- ((n1 - 1) * var(d1) + (n2 - 1) * var(d2))/(n1 +n2 - 2)
        if (sum(v^2) < 1/10^10){
            return(NA)
        } else {
            t2 <- n1 * n2 * dbar %*% solve(v) %*% dbar/(n1 +n2)
        }
        f <- (n1 + n2 - k - 1) * t2/((n1 + n2 - 2) * k)
        return(pf(f, k, n1 + n2 - k - 1, lower.tail = FALSE))
    }
    "ellipse2" <- function(loc, cov, alpha) {
        A <- cov
        detA <- A[1, 1] * A[2, 2] - A[1, 2]^2
        dist <- sqrt(qchisq(1 - alpha/2, 2))
        ylimit <- sqrt(A[2, 2]) * dist
        y <- seq(-ylimit, ylimit, 0.01 * ylimit)
        sqrt.discr <- sqrt(detA/A[2,2]^2*abs(A[2,2]*dist^2-y^2))
        sqrt.discr[c(1, length(sqrt.discr))] <- 0
        b <- loc[1] + A[1, 2]/A[2, 2] * y
        x1 <- b - sqrt.discr
        x2 <- b + sqrt.discr
        y <- loc[2] + y
        return(rbind(cbind(x1, y), cbind(rev(x2), rev(y))))
    }
################################################################################

    dataset[,col.j] <- as.factor(dataset[,col.j])
    juge <- levels(dataset[,col.j])
    nbjuge <- length(juge)
    dataset[,col.p] <- as.factor(dataset[,col.p])
    product <- levels(dataset[,col.p])
    nbprod <- length(product)
    info <- dataset[,c(col.j,col.p)]
    id.pos <- grep(id.recogn,colnames(dataset))
    intensity <- dataset[,id.pos-1]
    attribut <- colnames(intensity)
    nbatt <- length(attribut)
    ideal <- dataset[,id.pos]
    id.data <- cbind(info,ideal)
    int.data <- cbind(info,intensity)
    if (level.search.desc < 1){
        att.rm <- NULL
        res.aov <- vector("list",1)
        for (i in 1:nbatt){
            res.aov[1] <- as.matrix(summary(aov(int.data[,i+2]~int.data[,2]+int.data[,1])))
            if (res.aov[[1]][1,5]>level.search.desc){
                att.rm <- c(att.rm,i+2)
                print(paste("The attribute ",attribut[i]," is removed from the analysis.",sep=""))
            }
        }
        if (!is.null(att.rm)){
            int.data <- int.data[,-att.rm]
            id.data <- id.data[,-att.rm]
            attribut <- attribut[-(att.rm-2)]
            nbatt <- length(attribut)
            intensity <- intensity[,-(att.rm-2)]
            ideal <- ideal[,-(att.rm-2)]
        }
    }
    if (!correct){
        int.avg <- averagetable(int.data,formul=paste("~",colnames(info)[2],"+",colnames(info)[1],sep=""),firstvar=3)
        colnames(ideal) <- colnames(int.avg)
        data.pca <- rbind(int.avg,ideal)
    } else {
        int.p.avg <- averagetable(int.data,formul=paste("~",colnames(info)[2],"+",colnames(info)[1],sep=""),firstvar=3)
        int.j.avg <- averagetable(int.data,formul=paste("~",colnames(info)[1],"+",colnames(info)[2],sep=""),firstvar=3)
        ideal.juge <- vector("list",nbjuge)
        names(ideal.juge) <- juge
        for (j in 1:nbjuge){
            ideal.j <- id.data[id.data[,1]==juge[j],]
            rownames(ideal.j) <- ideal.j[,2]
            ideal.j <- ideal.j[,-c(1,2)]
            temp <- as.matrix(int.j.avg[j,])
            ideal.juge[[j]] <- sweep(ideal.j,2,as.vector(as.matrix(int.j.avg[j,])),FUN="-")
        }
        data.j.cplt <- matrix(0,0,nbatt)
        colnames(data.j.cplt) <- attribut
        l=0
        for (j in 1:nbjuge){
            data.j.cplt <- rbind(data.j.cplt,ideal.juge[[j]])
            rownames(data.j.cplt)[c((l+1):nrow(data.j.cplt))] <- paste(rownames(ideal.juge[[j]]),"_",juge[j],sep="")
            l=nrow(data.j.cplt)
        }
        colnames(data.j.cplt) <- colnames(int.p.avg)
        data.pca <- rbind.data.frame(scale(int.p.avg),data.j.cplt)
    }
    res.pca <- PCA(data.pca,ind.sup=c((nbprod+1):nrow(data.pca)),graph=F)
    ncp=0
    eig=res.pca$eig[1,1]
    while(eig>1 && ncp<nrow(res.pca$eig)){
        ncp=ncp+1
        eig <- res.pca$eig[ncp+1,1]
    }
    res.pca <- PCA(data.pca,ind.sup=c((nbprod+1):nrow(data.pca)),graph=F, ncp=max(ncp,2))
    for (i in 1:2)
        if (!coord[i] %in% c(1:ncp)){
            warning("The dimensions 1 and 2 will be used in this analysis!")
            coord=c(1,2)
        }
    ind.sup <- cbind(info,res.pca$ind.sup$coord)

    ideal.p.dim <- vector("list",nbprod)
    names(ideal.p.dim) <- product
    for (p in 1:nbprod){
        ideal.p.dim[[p]] <- ind.sup[ind.sup[,2]==product[p],-c(1,2)]
        rownames(ideal.p.dim[[p]]) <- paste(product[p],"_",juge,sep="")
    }
    if (is.null(nbchoix))
        nbchoix=nbjuge
    simul <- matrix(0,nbchoix,0)
    rownames(simul) <- paste("prod",1:nbchoix,sep="")
    for (sim in 1:nbsimul)
        simul <- cbind(simul,as.matrix(sample(nbprod,nbchoix,replace=T)))
    colnames(simul) <- paste("Simul.",1:nbsimul,sep="")
    coord.ellipse <- array(0,dim=c(nbprod,nbsimul,length(coord)))
    for (p in 1:nbprod)
        for (sim in 1:nbsimul)
            for (l in 1:length(coord))
                coord.ellipse[p,sim,l] <- mean(ideal.p.dim[[p]][simul[,sim],coord[l]])
    res.simul <- list()
    res.simul$sample <- t(simul)
    matP <- matrix(0,nbprod,length(coord)+1)
    rownames(matP) <- product
    colnames(matP) <- c(colnames(ideal.p.dim[[1]])[coord],"name.prod")
    for (p in 1:nbprod){
        matP[p,1:length(coord)] <- round(apply(ideal.p.dim[[p]][,coord],2,mean),3)
        matP[p,ncol(matP)] <- names(ideal.p.dim)[p]
    }
    matP <- as.data.frame(matP)
    for (l in 1:length(coord))
        matP[,l] <- as.numeric(as.character(matP[,l]))
    res.simul$moy$P <- matP
    matJP <- matrix(0,0,length(coord)+1)
    for (p in 1:nbprod)
        matJP <- rbind(matJP,cbind(ideal.p.dim[[p]][,coord],as.matrix(rep(product[p],nrow(ideal.p.dim[[p]])))))
    colnames(matJP) <- colnames(matP)
    matJP <- as.data.frame(matJP)
    for (l in 1:length(coord))
        matJP[,l] <- as.numeric(as.character(matJP[,l]))
    res.simul$moy$JP <- as.data.frame(matJP)
    res.simul2 <- matrix(0,0,length(coord)+1)
    for (p in 1:nbprod){
        res.simul2.temp <- matrix(0,nbsimul,0)
        for (l in 1:length(coord))
            res.simul2.temp <- cbind(res.simul2.temp,coord.ellipse[p,,l])
        res.simul2 <- rbind(res.simul2,cbind(res.simul2.temp,as.matrix(rep(product[p],nrow(res.simul2.temp)))))
    }
    colnames(res.simul2) <- colnames(matP)
    res.simul2 <- as.data.frame(res.simul2)
    for (l in 1:length(coord))
        res.simul2[,l] <- as.numeric(as.character(res.simul2[,l]))
    res.simul$moy$simul <- res.simul2
    mat <- res.simul
    alpha=0.05
    eig=res.pca$eig
    matP = cbind.data.frame(mat$moy$P[,1:2],mat$moy$P[,ncol(mat$moy$P)])
    matJP = cbind.data.frame(mat$moy$JP[,1:2],mat$moy$JP[,ncol(mat$moy$JP)])
    matsimul = cbind.data.frame(mat$moy$simul[,1:2],mat$moy$simul[,ncol(mat$moy$simul)])
    nbj <- nrow(matP)
    coord.ellipse.a.tracer <- matrix(0, 402, 2 * nbj)
    p <- 2
    for (i in 1:nbj) {
        VX <- var(matsimul[((i-1)*nbsimul+1):(i*nbsimul),1:2])
        coord.ellipse.a.tracer[,(1+2*(i-1)):(2*i)] <- ellipse2(as.numeric(t(matP[i,1:2])),VX,alpha)
    }
    minx <- min(min(coord.ellipse.a.tracer[,1+2*(0:(nbj-1))],na.rm=T),min(res.pca$ind$coord[,1]))
    maxx <- max(max(coord.ellipse.a.tracer[,1+2*(0:(nbj-1))],na.rm=T),max(res.pca$ind$coord[,1]))
    miny <- min(min(coord.ellipse.a.tracer[,2*(1:nbj)],na.rm=T),min(res.pca$ind$coord[,2]))
    maxy <- max(max(coord.ellipse.a.tracer[,2*(1:nbj)],na.rm=T),max(res.pca$ind$coord[,2]))
    dev.new()
    plot(res.pca,choix="ind",invisible="ind.sup",xlim=c(minx,maxx),ylim=c(miny,maxy),title="Single vs. Multiple Ideal",axes=coord)
    text(matP[,1],matP[,2],matP[,ncol(matP)],cex=0.7,pos=4,offset=0.2,col="blue3")
    for (p in 1:nbprod) {
        points(matP[p,1],matP[p,2],cex=1,col="blue3",pch=20)
        lines(coord.ellipse.a.tracer[,(1+2*(p-1)):(2*p)],col="blue3",lty=2)
    }
    res.hotelling <- matrix(1,nbprod,nbprod)
    rownames(res.hotelling) <- colnames(res.hotelling) <- product
    prod.sup <- averagetable(ind.sup,formul=paste("~",colnames(dataset)[col.p],"+",colnames(dataset)[col.j],sep=""),firstvar=3,method="mean")
    prod.sup <- cbind(prod.sup,as.matrix(rownames(prod.sup)),matrix("mean",nbprod,1))
    colnames(prod.sup)[(ncol(prod.sup)-1):ncol(prod.sup)] <- c("product","juge")
#    ind.sup <- cbind(ind.sup[,coord+2],ind.sup[,2],ind.sup[,1])
    ind.sup <- cbind(ind.sup[,-c(1,2)],ind.sup[,2],ind.sup[,1])
    colnames(ind.sup)[(ncol(ind.sup)-1):ncol(ind.sup)] <- c("product","juge")
    data.hotelling <- rbind(prod.sup,ind.sup)
    labprod <- data.hotelling[data.hotelling[,ncol(data.hotelling)]=="mean",ncol(data.hotelling)-1]
    aa = data.hotelling[-(1:nbprod), ]
    for (i in 1:(nbprod-1))
        for (j in (i+1):nbprod)
            res.hotelling[i,j]=res.hotelling[j,i]=hotelling2(aa[aa[,ncol(aa)-1]==labprod[i],1:(ncol(aa)-2)],aa[aa[,ncol(aa)-1]==labprod[j],1:(ncol(aa)-2)],nbchoix,nbchoix)
    return(res.hotelling)
}