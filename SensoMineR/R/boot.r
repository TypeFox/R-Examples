################ Function
boot <- function(X,method="sorting",axes=1:2,scale=TRUE,ncp=NULL,group=NULL,nbsim=200,level.conf = 0.95,nbchoix = NULL,color=NULL,cex=0.8, title=NULL,new.plot=TRUE){

 procrustes <- function(amat, target, orthogonal = FALSE, translate = FALSE, magnify = FALSE) {
        for (i in nrow(amat):1) {
            if (any(is.na(amat)[i, ]) | any(is.na(target)[i, ])) {
                amat <- amat[-i, ]
                target <- target[-i, ]
            }
        }
        dA <- dim(amat)
        dX <- dim(target)
        if (length(dA) != 2 || length(dX) != 2) stop("arguments amat and target must be matrices")
        if (any(dA != dX)) stop("dimensions of amat and target must match")
        if (length(attr(amat, "tmat"))) stop("oblique loadings matrix not allowed for amat")       
if (orthogonal) {
            if (translate) {
                p <- dX[1]
                target.m <- (rep(1/p, p) %*% target)[, ]
                amat.m <- (rep(1/p, p) %*% amat)[, ]
                target.c <- scale(target, center = target.m, scale = FALSE)
                amat.c <- scale(amat, center = amat.m, scale = FALSE)
                j <- svd(crossprod(target.c, amat.c))
            }
            else {
                amat.c <- amat
                j <- svd(crossprod(target, amat))
            }

            rot <- j$v %*% t(j$u)
            if (magnify)  beta <- sum(j$d)/sum(amat.c^2)
            else beta <- 1

            B <- beta * amat.c %*% rot
            if (translate)  B <- B + rep(as.vector(target.m), rep.int(p, dX[2]))
   
       value <- list(rmat = B, tmat = rot, magnify = beta)
            if (translate) value$translate <- target.m - (rot %*% amat.m)[, ]
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
##########End procrustes
  if (is.null(rownames(X))) rownames(X) <- 1:nrow(X)
  method <- tolower(method)
  ponder <- NULL
  
  if (method=="sorting") {
#    vrai <- MFA(X,group=rep(1,ncol(X)),type=rep("n",ncol(X)),graph=FALSE,ncp=Inf)
#    group <- vrai$call$group.mod
#    vrai$ind$coord = vrai$ind$coord/sqrt(length(group))
#    if (is.null(nbchoix)) nbchoix <- length(group)
#    X <- vrai$call$XTDC/sqrt(nbchoix)
    vrai <- MCA(X,ncp=Inf,graph=FALSE)
    group <- unlist(lapply(X,nlevels))
    if (is.null(nbchoix)) nbchoix <- length(group)
    tab.disj <- tab.disjonctif(X)
    X <- scale(tab.disj)*sqrt(nrow(X)/(nrow(X)-1))/sqrt(nbchoix)
    ponder <- 1-apply(tab.disj/nrow(X), 2, sum)
    if (is.null(title)) title <- "Confidence ellipses for sorting task"
  }

  if (method=="napping") {
    group <- rep(2,ncol(X)/2)
    type <- rep("c",length(group))
    if (is.null(title)) title <- "Confidence ellipses for the napping configuration"
    if (is.null(nbchoix)) nbchoix <- length(group)
    vrai <- MFA(X,group=group,type=type,graph=FALSE,ncp=Inf)
  }

  if (method=="sortnapping") {
    tab <- MFA(X[,1:3],group=c(2,1),type=c("c","n"),graph=FALSE,ncp=Inf)$ind$coord
    group <- ncol(tab)
    for (j in 2:(ncol(X)/3)){
      tab.aux <- MFA(X[,((j-1)*3+1):(j*3)],group=c(2,1),type=c("c","n"),graph=FALSE,ncp=Inf)$ind$coord
      tab <- cbind.data.frame(tab,tab.aux)
      group <- c(group,ncol(tab.aux))
    }
    X <- tab
    type <- rep("c",length(group))
    if (is.null(title)) title <- "Confidence ellipses for the sorted napping configuration"
    if (is.null(nbchoix)) nbchoix <- length(group)
    vrai <- MFA(X,group=group,type=type,graph=FALSE,ncp=Inf)
  }
  if (method=="freechoice") {
    if (scale) X <- scale(X)*sqrt(nrow(X)/(nrow(X)-1))
    type=rep("c",length(group))
    if (is.null(nbchoix)) nbchoix <- length(group)
    if (is.null(title)) title <- "Confidence ellipses for the free choice profiling"
    vrai <- MFA(X,group=group,type=type,graph=FALSE,ncp=Inf)
  }
  if (method=="hsort") {
    if (is.null(nbchoix)) nbchoix <- length(group)
    if (is.null(title)) title <- "Confidence ellipses for the hierarchical sorting task"
    type=rep("n",length(group))
    vrai <- MFA(X,group=group,type=type,graph=FALSE,ncp=Inf)
    X <- vrai$call$XTDC
    group <- vrai$call$group.mod
  }
  if (is.null(ponder)) ponder <- vrai$call$col.w
  estim.ncp <- estim_ncp(sweep(X,2,sqrt(ponder),FUN="*"),scale=FALSE,ncp.min=0,ncp.max=min(10,ncol(X)))
  if (is.null(ncp))  ncp <- max(estim.ncp$ncp,2,max(axes))
  else ncp <- min(ncp,ncol(vrai$ind$coord))
    
  listvar <- list()
  for (j in 1:length(group)) listvar[[j]] <- (cumsum(group)[j]-group[j]+1):cumsum(group)[j]
  jdd=vrai$ind$coord[,1:ncp]
        
  for (k in 1:nbsim){
    choix <- sample(1:length(group),nbchoix,replace=TRUE)
    auxi <- X[,unlist(listvar[choix])]
    ponder.auxi <- ponder[unlist(listvar[choix])]
    aux <- PCA(auxi,scale.unit=FALSE,graph=FALSE,col.w=ponder.auxi,ncp=ncp)$ind$coord
    aux <- procrustes(as.matrix(aux),as.matrix(vrai$ind$coord[,1:ncp]),orthogonal = TRUE, translate = TRUE,magnify=FALSE)$rmat
    colnames(aux) <- colnames(vrai$ind$coord)[1:ncp]
    jdd = rbind.data.frame(jdd,aux)
  }

#  res.pca=PCA(jdd,ind.sup=(nrow(X)+1):nrow(jdd),ncp=ncp,scale.unit=FALSE,graph=FALSE)

  if (new.plot) dev.new()
#  truc=cbind.data.frame(res.pca$ind.sup$coord,rep(rownames(X),nbsim))
  truc=cbind.data.frame(jdd[-(1:nrow(X)),],rep(rownames(X),nbsim))
  simul=list()
  simul$moy$simul= truc[order(truc[,ncol(truc)]),]
  simul$moy$P=cbind.data.frame(vrai$ind$coord,rownames(X))
  simul$moy$P=simul$moy$P[order(simul$moy$P[,ncol(simul$moy$P)]),]
  plotellipse (simul, alpha = 1-level.conf, eig = signif(vrai$eig,4),coord=axes,cex=cex,color=color,title=title)
  return(list(simul=simul,estim.ncp=estim.ncp))
}
