plot.MIPCA <- function(x,choice="all",axes=c(1,2),new.plot=TRUE,main=NULL,level.conf=0.95, ...){

####
 procrustes <- function(amat, target, orthogonal = FALSE, translate = FALSE,
        magnify = FALSE) {
        for (i in nrow(amat):1) {
            if (any(is.na(amat)[i, ]) | any(is.na(target)[i,
                ])) {
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
                target.c <- scale(target, center = target.m,
                  scale = FALSE)
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
                B <- B + rep(as.vector(target.m), rep.int(p,
                  dX[2]))
   
       value <- list(rmat = B, tmat = rot, magnify = beta)
            if (translate)
                value$translate <- target.m - (rot %*% amat.m)[,
                  ]
    
  }

        else {
            b <- solve(amat, target)
            gamma <- sqrt(diag(solve(crossprod(b))))
            rot <- b * rep(gamma, rep.int(dim(b)[1], length(gamma)))
            B <- amat %*% rot
            fcor <- solve(crossprod(rot))
            value <- list(rmat = B, tmat = rot, correlation = fcor)
        }

        return(value)
    }
####
  res <- x
  if (!inherits(res, "MIPCA")) stop("non convenient data")
  ncp <- res$call$ncp
  reference <- PCA(res$res.imputePCA,scale.unit=res$call$scale,graph=FALSE,ncp=ncp)
  rec.pca <- res$res.imputePCA
#  rec <- reconst(reference,ncp)
#  rec.pca <- as.matrix(res$call$X)
#  rec.pca[res$call$missing] <- rec[res$call$missing]
  
  res.var <- res.supp <- rec.pca
  res.procrustes <- reference$ind$coord[,1:ncp]
  res.dim <- as.matrix(res$res.imputePCA)

##for (i in 1:dim(res$res.MI)[3]){
## rec.pca <- res$res.MI[,,i]
for (i in 1:length(res$res.MI)){
 rec.pca <- res$res.MI[[i]]
 acpfin <- PCA(rec.pca, scale.unit=res$call$scale,graph=FALSE,ncp=ncp)

 tourne <- procrustes(acpfin$ind$coord[,1:ncp], reference$ind$coord[,1:ncp],orthogonal = TRUE, translate = TRUE, magnify = TRUE)$rmat

 colnames(tourne) <- colnames(res.procrustes)
 colnames(rec.pca) <- colnames(res.supp)
 res.procrustes <- rbind.data.frame(res.procrustes,tourne)
 res.supp <- rbind.data.frame(res.supp,rec.pca)
 res.var <- cbind.data.frame(res.var,rec.pca)
 res.dim <- cbind.data.frame(res.dim,acpfin$ind$coord[,1:ncp])
}

####
if (!is.null(main)) title <- main
if ((choice=="all")|(choice=="ind.proc")){
  if (new.plot) dev.new()
  oo=PCA(res.procrustes,ind.sup=c((nrow(res$call$X)+1):nrow(res.procrustes)),scale.unit=FALSE,graph=FALSE)
  oo$eig=reference$eig
el=coord.ellipse(cbind.data.frame(as.factor(rep(rownames(res$call$X),res$call$nboot)),oo$ind.sup$coord[,axes]),level.conf=level.conf) 
  if (is.null(main)) title="Multiple imputation using Procrustes" 
  plot(oo,axes=axes,col.ind.sup=rep(1:nrow(res$call$X),res$call$nboot),label="ind",ellipse=el,col.quali="black", title=title,invisible="ind.sup",new.plot=FALSE)

#  if (!is.null(add.tab)){
#    vrai = PCA(add.tab,graph=FALSE,scale=res$call$scale)
#    tourne <- procrustes(vrai$ind$coord[,axes], reference$ind$coord[,axes],orthogonal = TRUE, translate = TRUE, magnify = TRUE)$rmat
#    points(tourne[,axes],cex=0.9,col=2)
#  }
}

if ((choice=="all")|(choice=="dim")){
  if ((new.plot)&!nzchar(Sys.getenv("RSTUDIO_USER_IDENTITY"))) dev.new()
  colnames(res.dim)=paste("V",1:ncol(res.dim))
  ooo=PCA(res.dim,quanti.sup=(ncol(res$call$X)+1):ncol(res.dim),scale.unit=res$call$scale,graph=FALSE)
  ooo$eig=reference$eig
  if (is.null(main)) title="Projection of the Principal Components"  
  plot(ooo,choi="var",axes=axes,title=title,label="none",new.plot=FALSE,invisible="var")
}

if ((choice=="all")|(choice=="ind.supp")){
  if ((new.plot)&!nzchar(Sys.getenv("RSTUDIO_USER_IDENTITY"))) dev.new()
  oo=PCA(res.supp,ind.sup=c((nrow(res$call$X)+1):nrow(res.supp)),scale.unit=res$call$scale,graph=FALSE)
  el=coord.ellipse(cbind.data.frame(as.factor(rep(rownames(res$call$X),res$call$nboot)),oo$ind.sup$coord[,1:2]),level.conf = level.conf) 
  if (is.null(main)) title="Supplementary projection"    
  plot(oo,axes=axes,col.ind.sup=rep(1:nrow(res$call$X),res$call$nboot),label="ind",ellipse=el,col.quali="black",
    title=title,invisible="ind.sup",new.plot=FALSE)
#  if (!is.null(add.tab)){
#    dele = PCA(rbind.data.frame(rec.pca,add.tab),ind.sup=c((nrow(res$call$X)+1):(2*nrow(res$call$X))),scale=scale,graph=FALSE)
#    points(dele$ind.sup$coord[,axes],col=2)
#  }
}

if ((choice=="all")|(choice=="var")){
  if ((new.plot)&!nzchar(Sys.getenv("RSTUDIO_USER_IDENTITY"))) dev.new()
  color = c("black", "red", "green3", "blue", "cyan", "magenta", 
            "darkgray", "darkgoldenrod", "darkgreen", "violet", 
            "turquoise", "orange", "lightpink", "lavender", "yellow", 
            "lightgreen", "lightgrey", "lightblue", "darkkhaki", 
            "darkmagenta", "darkolivegreen", "lightcyan", "darkorange", 
            "darkorchid", "darkred", "darksalmon", "darkseagreen", 
            "darkslateblue", "darkslategray", "darkslategrey", 
            "darkturquoise", "darkviolet", "lightgray", "lightsalmon", 
            "lightyellow", "maroon")
  colnames(res.var)=paste("V",1:ncol(res.var))
  colnames(res.var)[1:ncol(res$call$X)]=colnames(res$call$X)
  oo=PCA(res.var,quanti.sup=c((ncol(res$call$X)+1):ncol(res.var)),scale.unit=res$call$scale,graph=FALSE)
  if (is.null(main)) title="Variable representation"    
  plot(oo, axes=axes, choix = "var", title=title,invisible = "quanti.sup", col.hab = color[1:ncol(res$call$X)],new.plot=FALSE)
  for (k in 1:res$call$nboot) points(oo$quanti.sup$coord[((k-1)*ncol(res$call$X)+1):(k*ncol(res$call$X)),axes[1]], oo$quanti.sup$coord[((k-1)*ncol(res$call$X)+1):(k*ncol(res$call$X)),axes[2]], col = color[1:ncol(res$call$X)], pch = 15, cex = 0.3)
}
}
