IdMap <- function(dataset,col.p,col.j,col.lik,id.recogn,nbsimul=500,nbchoix=NULL,alpha=0.05,coord=c(1,2),precision=0.1,levels.contour=NULL,color=FALSE,cons.eq=FALSE){
################################################################################
forprefmap <- function(dataset,col.p,col.j,col.lik,id.recogn,rm.j=TRUE){
    dataset[,col.p] <- as.factor(dataset[,col.p])
    product <- levels(dataset[,col.p])
    nbprod <- length(product)
    dataset[,col.j] <- as.factor(dataset[,col.j])
    juge <- levels(dataset[,col.j])
    nbjuge <- length(juge)
    desc <- dataset[,c(col.j,col.p)]
    id.pos <- grep(id.recogn,colnames(dataset))
    intensity <- dataset[,id.pos-1]
    int.data <- cbind(desc,intensity)
    liking <- as.matrix(dataset[,col.lik])
    colnames(liking) <- "liking"
    lik.data <- cbind(desc,liking)
    int.avg <- averagetable(int.data,formul=paste("~",colnames(dataset)[col.p],"+",colnames(dataset)[col.j],sep=""),firstvar=3,method="coeff")
    lik.j <- matrix(0,nbprod,0)
    rownames(lik.j) <- product
    for (j in 1:nbjuge){
        lik.j.temp <- lik.data[lik.data[,1]==juge[j],]
        lik.j.temp.rn <- lik.j.temp[,2]
        lik.j.temp <- as.matrix(lik.j.temp[,3])
        rownames(lik.j.temp) <- lik.j.temp.rn
        lik.j <- merge(lik.j,lik.j.temp,all=T,by=0,sort=F)
        lik.j.rn <- lik.j[,1]
        lik.j <- as.matrix(lik.j[,-1])
        rownames(lik.j) <- lik.j.rn
        colnames(lik.j)[ncol(lik.j)] <- juge[j]
    }
    if (rm.j){
        rm.c <- NULL
        for (j in 1:ncol(lik.j))
            if (var(na.omit(lik.j[,j]))==0)
                rm.c <- c(rm.c,j)
        if (!is.null(rm.c)){
            lik.j <- lik.j[,-rm.c]
            juge.rm <- juge[rm.c]
        } else {
            juge.rm <- NULL
        }
    }
    res <- list()
    res$senso <- int.avg
    res$hedo <- as.data.frame(lik.j)
    if (rm.j)
        juge.rm
    return(res)
}
################################################################################
procrustes <- function(amat, target, orthogonal = FALSE, translate = FALSE, magnify = FALSE) {
    for (i in nrow(amat):1) {
        if (any(is.na(amat)[i, ]) | any(is.na(target)[i,])) {
            amat <- amat[-i, ]
            target <- target[-i, ]
        }
    }
    dA <- dim(amat)
    dX <- dim(target)
    if (length(dA) != 2 || length(dX) != 2)
        stop("arguments amat and target must be matrices")
    if (any(dA != dX))
        stop("dimensions of amat and target must match")
    if (length(attr(amat, "tmat")))
        stop("oblique loadings matrix not allowed for amat")
    if (orthogonal) {
        if (translate) {
            p <- dX[1]
            target.m <- (rep(1/p, p) %*% target)[, ]
            amat.m <- (rep(1/p, p) %*% amat)[, ]
            target.c <- scale(target, center = target.m,scale = FALSE)
            amat.c <- scale(amat, center = amat.m, scale = FALSE)
            j <- svd(crossprod(target.c, amat.c))
        }
        else {
            amat.c <- amat
            j <- svd(crossprod(target, amat))
        }
        rot <- j$v %*% t(j$u)
        if (magnify)
            beta <- sum(j$d)/sum(amat.c^2)
        else beta <- 1
        B <- beta * amat.c %*% rot
        if (translate)
            B <- B + rep(as.vector(target.m), rep.int(p,dX[2]))
        value <- list(rmat = B, tmat = rot, magnify = beta)
        if (translate)
            value$translate <- target.m - (rot %*% amat.m)[,]
    } else {
        b <- solve(amat, target)
        gamma <- sqrt(diag(solve(crossprod(b))))
        rot <- b * rep(gamma, rep.int(dim(b)[1], length(gamma)))
        B <- amat %*% rot
        fcor <- solve(crossprod(rot))
        value <- list(rmat = B, tmat = rot, correlation = fcor)
    }
    return(value)
}
################################################################################
    oo <- order(dataset[,col.p])
    dataset <- dataset[oo,]
    oo <- order(dataset[,col.j])
    dataset <- dataset[oo,]
    dataset[,col.p] <- as.factor(dataset[,col.p])
    product <- levels(dataset[,col.p])
    nbprod <- length(product)
    dataset[,col.j] <- as.factor(dataset[,col.j])
    juge <- levels(dataset[,col.j])
    nbjuge <- length(juge)
    info <- dataset[,c(col.j,col.p)]
    id.pos <- grep(id.recogn,colnames(dataset))
    intensity <- dataset[,id.pos-1]
    attribut <- colnames(intensity)
    nbatt <- length(attribut)
    ideal <- dataset[,id.pos]
    int.data <- cbind(info,intensity)
    id.data <- cbind(info,ideal)
    data.cut <- forprefmap(dataset,col.p=col.p,col.j=col.j,id.recogn=id.recogn,col.lik=col.lik)
    int.p.avg <- scale(data.cut$senso,scale=F)
    int.j.avg <- averagetable(int.data,formul=paste("~",colnames(dataset)[col.j],"+",colnames(dataset)[col.p]),firstvar=3,method="coeff")
    id.j.avg <- averagetable(id.data,formul=paste("~",colnames(dataset)[col.j],"+",colnames(dataset)[col.p]),firstvar=3,method="coeff")
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
    data.pca <- rbind.data.frame(int.p.avg,data.j.cplt)
    res.pca <- PCA(data.pca,ind.sup=c((nbprod+1):nrow(data.pca)),graph=F,ncp=Inf)

    id.j.avg.cor <- id.j.avg-int.j.avg
    colnames(id.j.avg.cor) <- colnames(int.p.avg)
    data.pcab <- rbind.data.frame(int.p.avg,id.j.avg.cor)
    res.pcab <- PCA(data.pcab,ind.sup=c((nbprod+1):nrow(data.pcab)),graph=F,ncp=Inf)
    plot.PCA(res.pcab,choix="ind",cex=0.8,label="ind.sup",new.plot=T,title="Projection of the individual averaged ideal profiles",axes=coord)
    eig <- res.pca$eig

    data.pca2 <- merge(data.cut$senso,data.cut$hedo,all=T,by=0,sort=F)
    rownames(data.pca2) <- data.pca2[,1]
    data.pca2 <- data.pca2[,-1]
    res.pca2 <- PCA(data.pca2,quanti.sup=(ncol(data.cut$senso)+1):ncol(data.pca2),graph=F)
    plot.PCA(res.pca2,choix="var",invisible="var",label="quanti.sup",cex=0.9,new.plot=T,title="Projection of the individual hedonic scores",axes=coord)
#    layout(matrix(1:2,1,2))                                                                            
#    plot.PCA(res.pca,choix="ind",label="none",new.plot=T)
#    plot.PCA(res.pca2,choix="var",invisible="var",label="none",new.plot=T)

    ideal.j.dim <- vector("list",nbjuge)
    names(ideal.j.dim) <- juge
    l=0
    for (j in 1:nbjuge){
        ideal.j.dim[[j]] <- res.pca$ind.sup$coord[c((l+1):(l+nrow(ideal.juge[[j]]))),]
        rownames(ideal.j.dim[[j]]) <- paste(product,"_",juge[j],sep="")
        l <- l+nrow(ideal.juge[[j]])
    }
    if (is.null(nbchoix))
        nbchoix=nbprod
    simul <- matrix(0,nbchoix,0)
    rownames(simul) <- paste("prod",1:nbchoix,sep="")
    for (sim in 1:nbsimul)
        simul <- cbind(simul,as.matrix(sample(nbprod,nbchoix,replace=T)))
    colnames(simul) <- paste("Simul.",1:nbsimul,sep="")

    ponder=res.pcab$call$col.w
    estim.ncp <- max(max(coord),estim_ncp(sweep(int.p.avg, 2, sqrt(ponder), FUN = "*"),scale = FALSE, ncp.min = 0, ncp.max = min(10, ncol(int.p.avg)))$ncp)
    estim.ncp <- max(estim.ncp,max(coord))
        
    target.pca <- res.pcab$ind.sup$coord[,1:estim.ncp]
    jdd <- target.pca
    for (sim in 1:nbsimul){
        juge.sample <- matrix(0,nbjuge,nbatt)
        rownames(juge.sample) <- juge
        colnames(juge.sample) <- colnames(int.p.avg)
        for (j in 1:nbjuge)
            juge.sample[j,] <- apply(ideal.juge[[j]][as.vector(simul[,sim]),],2,mean)
        data.pca.temp <- rbind.data.frame(int.p.avg,juge.sample)
        res.pca.temp <- PCA(data.pca.temp,ind.sup=c((nbprod+1):nrow(data.pca.temp)),graph=F,ncp=estim.ncp)$ind.sup$coord
        aux <- procrustes(res.pca.temp,target.pca, orthogonal = TRUE, translate = TRUE, magnify = FALSE)$rmat
        colnames(aux) <- colnames(jdd)
        jdd = rbind.data.frame(jdd, aux)    
    }
        
    truc = cbind.data.frame(jdd[-(1:nrow(target.pca)),],rep(rownames(target.pca),nbsimul))
    res.simul = list()
    res.simul$moy$simul = truc[order(truc[,ncol(truc)]),]
    res.simul$moy$J = cbind.data.frame(target.pca, rownames(target.pca))
    res.simul$moy$J = res.simul$moy$J[order(res.simul$moy$J[, ncol(res.simul$moy$J)]),]

################################################################################
plotellipse2 <- function(mat,alpha=0.05,coord=c(1,2),eig,cex=1,color=NULL){
    res <- plotellipseinter2(mat,alpha=alpha,coord=coord,nbgroup=1,moy=T,eig=eig,color=color,cex=cex)
    if (length(mat$partiel) != 0) {
        dev.new()
        nbgroup=length(levels(mat$partiel$simul[,ncol(mat$partiel$simul)]))/length(levels(mat$moy$simul[,ncol(mat$moy$simul)]))
        plotellipseinter2(mat,alpha=alpha,coord=coord,nbgroup=nbgroup,moy=F,eig=eig,color=color,cex=cex)
    }
    return(res)
}
plotellipseinter2 <- function(mat,alpha=0.05,coord=c(1,2),nbgroup=1,moy=TRUE,eig,cex=1,color=NULL){
    if (moy == T){
        matJ = cbind.data.frame(mat$moy$J[,coord],mat$moy$J[,ncol(mat$moy$J)])
        matJP = cbind.data.frame(mat$moy$JP[,coord],mat$moy$JP[,ncol(mat$moy$JP)])
        matsimul = cbind.data.frame(mat$moy$simul[,coord],mat$moy$simul[,ncol(mat$moy$simul)])
    }
    if (moy == F){
        matmoyJ = cbind.data.frame(mat$moy$J[,coord],mat$moy$J[,ncol(mat$moy$J)])
        matmoyJP = cbind.data.frame(mat$moy$JP[,coord],mat$moy$JP[,ncol(mat$moy$JP)])
        matmoysimul = cbind.data.frame(mat$moy$simul[,coord],mat$moy$simul[,ncol(mat$moy$simul)])
        matJ=cbind.data.frame(mat$partiel$J[,coord],mat$partiel$J[,ncol(mat$partiel$J)])
        matJP = cbind.data.frame(mat$partiel$JP[,coord],mat$partiel$JP[,ncol(mat$partiel$JP)])
        matsimul = cbind.data.frame(mat$partiel$simul[,coord],mat$partiel$simul[,ncol(mat$partiel$simul)])
    }
    nbp <- nrow(matJ)
    nbjuge <- nbp/nbgroup
    coord.ellipse.a.tracer <- matrix(0, 402, 2 * nbp)
    p <- 2
    nbprod <- nrow(matJP)/nrow(matJ)
    nbsimul <- nrow(matsimul)/nrow(matJ)
    for (i in 1:nbp) {
        VX <- var(matsimul[((i-1)*nbsimul+1):(i*nbsimul),1:2])
        coord.ellipse.a.tracer[,(1+2*(i-1)):(2*i)] <- ellipse2(as.numeric(t(matJ[i,1:2])),VX,alpha)
    }
    minx <- min(coord.ellipse.a.tracer[,1+2*(0:(nbp-1))],na.rm=T)
    maxx <- max(coord.ellipse.a.tracer[,1+2*(0:(nbp-1))],na.rm=T)
    miny <- min(coord.ellipse.a.tracer[,2*(1:nbp)],na.rm=T)
    maxy <- max(coord.ellipse.a.tracer[,2*(1:nbp)],na.rm=T)
    plot(0,0,xlab=paste("Dim ",coord[1]," (",round(eig[coord[1],2],2),"%)",sep=""),ylab=paste("Dim ",coord[2]," (",round(eig[coord[2],2],2),"%)",sep=""),xlim=c(minx*1.05,maxx*1.05),ylim=c(1.05*miny,1.05*maxy),col="white",asp=1)
    if (moy == T)
        title(main = "Individual ideal confidence ellipses")
    abline(v = 0, lty = 2)
    abline(h = 0, lty = 2)
    if (moy == F){
        points(matmoyJ[,1],matmoyJ[,2],cex=0.8*cex,col="blue3",pch=15)
        text(matmoyJ[,1],matmoyJ[,2],matmoyJ[,ncol(matmoyJ)],cex=0.8*cex,pos=4,offset=0.2,col="blue3")
    }
    if (moy == T)
        text(matJ[,1],matJ[,2],matJ[,ncol(matJ)],cex=0.8*cex,pos=4,offset=0.2,col="blue3")
    for (j in 1:nbgroup) {
        for (i in 1:nbjuge) {
            points(matJ[(j-1)*nbjuge+i,1],matJ[(j-1)*nbjuge+i,2],cex=0.8*cex,col="blue3",pch=20)
            if (moy == F)
                lines(c(matJ[(j-1)*nbjuge+i,1],matmoyJ[i,1]),c(matJ[(j-1)*nbjuge+i,2],matmoyJ[i,2]),col="lightblue",lty=j)
            lines(coord.ellipse.a.tracer[,(1+2*((i+(j-1)*nbjuge)-1)):(2*(i+(j-1)*nbjuge))],col="lightblue",lty=j)
        }
    }
    return(coord.ellipse.a.tracer)
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

    dev.new()
    res.ellipse <- round(plotellipse2(res.simul,alpha=0.05,coord=coord,eig=res.pca$eig,cex=1,color=NULL),1)
    minx <- min(min(res.ellipse[,1+2*(0:(nbjuge-1))],na.rm=T),floor(min(res.pca$ind$coord[,coord[1]])))
    maxx <- max(max(res.ellipse[,1+2*(0:(nbjuge-1))],na.rm=T),ceiling(max(res.pca$ind$coord[,coord[1]])))
    miny <- min(min(res.ellipse[,2*(1:nbjuge)],na.rm=T),floor(min(res.pca$ind$coord[,coord[2]])))
    maxy <- max(max(res.ellipse[,2*(1:nbjuge)],na.rm=T),ceiling(max(res.pca$ind$coord[,coord[2]])))
    juge.mat <- vector("list",nbjuge)
    names(juge.mat) <- juge
    lim.minx <- floor(minx)
    lim.maxx <- ceiling(maxx)
    lim.miny <- floor(miny)
    lim.maxy <- ceiling(maxy)
    nbrow <- length(seq(lim.minx,lim.maxx,precision))
    nbcol <- length(seq(lim.miny,lim.maxy,precision))
    juge.tot <- matrix(0,nbcol,nbrow)
    rownames(juge.tot) <- round(seq(lim.miny,lim.maxy,precision),1)
    colnames(juge.tot) <- round(seq(lim.minx,lim.maxx,precision),1)
    cons.wgt <- matrix(0,1,nbjuge)
    rownames(cons.wgt) <- "weight"
    colnames(cons.wgt) <- juge
    for (j in 1:nbjuge){
        juge.mat[[j]] <- matrix(0,nbcol,nbrow)
        rownames(juge.mat[[j]]) <- round(seq(lim.miny,lim.maxy,precision),1)
        colnames(juge.mat[[j]]) <- round(seq(lim.minx,lim.maxx,precision),1)
        ellipse.x <- res.ellipse[,(1+2*(j-1))]
        ellipse.y <- res.ellipse[,(2*j)]
        for (i in 1:length(ellipse.x)){
            posx <- grep(ellipse.x[i],colnames(juge.mat[[j]]))
            if (length(posx)>1){
                posx.temp=NULL
                for (l in 1:length(posx))
                    if (colnames(juge.mat[[j]])[posx[l]]==ellipse.x[i])
                        posx.temp=posx[l]
                if (!is.null(posx.temp)){
                    posx=posx.temp
                } else {
                    stop("Not convenient posx definition")
                }
            }
            posy <- grep(ellipse.y[i],rownames(juge.mat[[j]]))
            if (length(posy)>1){
                posy.temp=NULL
                for (l in 1:length(posy))
                    if (rownames(juge.mat[[j]])[posy[l]]==ellipse.y[i])
                        posy.temp=posy[l]
                if (!is.null(posy.temp)){
                    posy=posy.temp
                } else {
                    stop("Not convenient posy definition")
                }
            }
            juge.mat[[j]][posy,posx]=1
        }
        for (i in 1:nrow(juge.mat[[j]])){
            pos1 <- grep(1,juge.mat[[j]][i,])
            if (length(pos1)>=2)
                juge.mat[[j]][i,c(pos1[1]:pos1[length(pos1)])]=1
        }
        if (cons.eq){
            if (sum(juge.mat[[j]]) > max(10,1/precision))
                cons.wgt[1,j] <- 1/sum(juge.mat[[j]])
        } else {
            cons.wgt[1,j] <- 1
        }
        juge.mat[[j]] <- juge.mat[[j]]*cons.wgt[1,j]
        juge.tot=juge.tot+juge.mat[[j]]
    }
    juge.tot.rn <- paste("Y_",rownames(juge.tot),sep="")
    juge.tot.cn <- paste("X_",colnames(juge.tot),sep="")
    rownames(juge.tot) <- juge.tot.rn
    colnames(juge.tot) <- juge.tot.cn
    juge.tot <- 100*round(juge.tot/sum(cons.wgt[1,]),3)
    f1 <- seq(lim.minx,lim.maxx,precision)
    f2 <- seq(lim.miny,lim.maxy,precision)
    if (!is.null(levels.contour)){
        if (min(levels.contour)<0 || max(levels.contour)>100 || length(levels.contour)<2){
            warning("Not convenient 'levels.contour' definition: the default value will be used")
            levels.contour=NULL
        } else {
            oo <- order(levels.contour)
            levels.contour <- levels.contour[oo]
        }
    }
    if (is.null(levels.contour))
        levels.contour <- seq(10,5*floor(max(juge.tot)/5),5)
    dev.new()
    if (cons.eq){
        titre <- "Weighted Ideal Mapping"
    } else {
        titre <- "Ideal Mapping"
    }
    if (color){
        image(f1,f2,t(juge.tot),col=terrain.colors(200),xlab=paste("Dim ",coord[1],"(",round(res.pca$eig[coord[1],2],2),"%)",sep=""),ylab=paste("Dim ",coord[2],"(",round(res.pca$eig[coord[2],2],2),"%)",sep=""),main=titre)
        contour(f1,f2,t(juge.tot),nlevels=length(levels.contour),levels=levels.contour,add=T,labex=0)
        for (i in 1:nrow(res.pca$ind$coord)) {
            points(res.pca$ind$coord[i,coord[1]], res.pca$ind$coord[i,coord[2]],pch=15)
            text(res.pca$ind$coord[i,coord[1]],res.pca$ind$coord[i,coord[2]],rownames(res.pca$ind$coord)[i],pos=4,offset=0.2,cex=0.7)
        }
        abline(v=0,lty=2)
        abline(h=0,lty=2)
    } else {
        image(f1,f2,t(juge.tot),col=grey(1:max(juge.tot)/100),xlab=paste("Dim ",coord[1],"(",round(res.pca$eig[coord[1],2],2),"%)",sep=""),ylab=paste("Dim ",coord[2],"(",round(res.pca$eig[coord[2],2],2),"%)",sep=""),main=titre)
        contour(f1,f2,t(juge.tot),nlevels=length(levels.contour),levels=levels.contour,add=T,labex=0,col="white")
        for (i in 1:nrow(res.pca$ind$coord)) {
            points(res.pca$ind$coord[i,coord[1]], res.pca$ind$coord[i,coord[2]],pch=15,col="white")
            text(res.pca$ind$coord[i,coord[1]],res.pca$ind$coord[i,coord[2]],rownames(res.pca$ind$coord)[i],pos=4,offset=0.2,cex=0.7,col="white")
        }
        abline(v=0,lty=2,col="white")
        abline(h=0,lty=2,col="white")
    }
    maxval <- max(juge.tot)
    res.id <- matrix(0,0,2)
    colnames(res.id) <- c("X","Y")
    for (i in 1:nrow(juge.tot))
        for (j in 1:ncol(juge.tot))
            if (juge.tot[i,j]==maxval)
                res.id <- rbind(res.id,matrix(c(i,j),1,2))
#                res.id <- rbind(res.id,matrix(c(explode.list(colnames(juge.tot)[j],separator="_")[2],explode.list(rownames(juge.tot)[i],separator="_")[2]),1,2))
    rownames(res.id) <- paste("Ideal_",LETTERS[1:nrow(res.id)],sep="")
    id.profile <- matrix(0,0,nbatt)
    colnames(id.profile) <- attribut
    juge.max <- vector("list",nrow(res.id))
    names(juge.max) <- rownames(res.id)
    for (i in 1:nrow(res.id))
        for (j in 1:nbjuge)
            if (juge.mat[[j]][res.id[i,1],res.id[i,2]]==1)
                juge.max[[i]] <- c(juge.max[[i]],j)
    id.j.avg <- averagetable(id.data,formul=paste("~",colnames(dataset)[col.j],"+",colnames(dataset)[col.p]),firstvar=3,method="mean")
    for (i in 1:nrow(res.id))
        id.profile <- rbind(id.profile,t(as.matrix(apply(id.j.avg[juge.max[[i]],],2,mean))))
#        id.profile <- rbind(id.profile,ideaux.carto(int.p.avg,score=as.numeric(as.character(res.id[i,])),loading=res.pca$var$coord[,coord],eigenvalues=res.pca$eig[coord,1],correlation=T))
    id.profile <- t(as.matrix(apply(id.profile,2,mean)))
    rownames(id.profile) <- "Ideal"
    res <- vector("list")
    res$PCA <- res.pca
    res$PCA$data <- data.pca
    res$PCA$dim <- coord
#    res$simul <- simul
    res$idmap$data <- juge.tot
    res$idmap$j.weight <- cons.wgt
    res$idmap$precision <- precision
#    res$ideal$position <- res.id
    res$ideal$profiles <- id.profile
    res$ideal$pct.conso <- maxval
    class(res) <- c("IdMap","list")
    return(res)
}