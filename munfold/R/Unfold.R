unfold <- function(x,...) UseMethod("unfold")

unfold.formula <- function(x,data=parent.frame(),...){

    latent.dims <- all.vars(x[c(1,3)])
    ndims <- length(latent.dims)
    formula <- x[1:2]
    formula <- update(formula,~.-1)
    cl <- match.call()
    m <- match.call(expand.dots=FALSE)
    mf <- m[c(1,match(c("formula","data","subset"),
            names(m),nomatch=0L))]
    mf[[1]] <- as.name("model.frame")
    mf$formula <- formula
    mf <- eval(mf,parent.frame())
    X <- model.matrix(formula,mf)

    res <- unfold.matrix(X,ndims=length(latent.dims),...)
    ndim <- ncol(res$A)
    colnames(res$A) <- latent.dims[1:ndim]
    colnames(res$B) <- latent.dims[1:ndim]

    attr(res,"plot_discrete") <- c(FALSE,TRUE)
    attr(res,"biplot_type") <- c("density","text")
    attr(res,"procrustes_use") <- "B"
    res
}

unfold.matrix <- function(x,ndims=NULL,squared=FALSE,tol=1e-7,
          method=c("Schoenemann","CG"),...){
    D <- x
    method <- match.arg(method)
    if(ncol(D)>nrow(D)) {
            D<-t(D)
            D.transposed <- TRUE
        }
    else D.transposed <- FALSE
    if(!squared) D <- D^2
    I <- nrow(D)
    J <- ncol(D)
    D.star <- sweep(D,1,rowMeans(D),"-")
    D.star <- sweep(D.star,2,colMeans(D.star),"-")/2
    D.svd <- svd(D.star)
    usable <- which(D.svd$d/max(D.svd$d) > tol)
    if(length(ndims)){
      if(length(usable) <ndims){
        warning("Number of singular values above threshold less than requested dimension")
        ndims <- length(usable)
        }
      }
    else ndims <- length(usable)

    opts <- options(warn=-1)
    uf.start <- schoenemann.unfolding(D,m=ndims)
    options(opts)
    if(method=="Schoenemann") return(uf.start)
    if(missing(ndims))
      ndims <- uf.start$dim
    A.start <- uf.start$A
    B.start <- uf.start$B

    SSQ.start <- full.uf.ssq(D,A=A.start,B=B.start)
    cat("\nStress for starting values",SSQ.start/sum(D^2))

    translate.objective <- function(par){
      t0 <- par[1:ndims]
      uf.ssq(D,t=rep(1,ndims),t0,A.start,B.start)
    }

    translate.gradient <- function(par){
      t0 <- par[1:ndims]
      trans.grad(D,t=rep(1,ndims),t0,A.start,B.start)
    }
    t0.start <- rep(0,ndims)
    translate.optim.res <- optim(t0.start,
        fn=translate.objective,
        gr=translate.gradient,
        method=method,
        ...
      )

    if(translate.optim.res$convergence>0){
      warning("Translation algorithm did not converge. Return code ",translate.optim.res$convergence)
      if(length(translate.optim.res$message))
        warning("Message: ",translate.optim.res$message)
      }
    t0 <- translate.optim.res$par
    A.trans <- sweep(A.start,2,t0,"+")
    B.trans <- sweep(B.start,2,t0,"-")

    SSQ.trans <- uf.ssq(D,t=rep(1,ndims),t0,A.start,B.start)
    cat("\nStress after translation",SSQ.trans/sum(D^2))

    rescale.cg.objective <- function(par){
      uf.ssq(D,par,t0,A.start,B.start)
    }

    rescale.cg.gradient <- function(par){
      t.t0 <- uf.grad(D,par,t0,A.start,B.start)
      t.t0[[1]]
    }


    t.start <- rep(1,ndims)

    rescale.optim.res <- optim(t.start,
        fn=rescale.cg.objective,
        gr=rescale.cg.gradient,
        method=method,
        ...
      )

    if(rescale.optim.res$convergence>0){
      warning("Algorithm did not converge. Return code ",rescale.optim.res$convergence)
      if(length(rescale.optim.res$message))
        warning("Message: ",rescale.optim.res$message)
      }
    t <- rescale.optim.res$par

    A.rescaled <- sweep(sweep(A.start,2,t,"*"),2,t0,"+")
    B.rescaled <- sweep(sweep(B.start,2,t,"/"),2,t0,"-")

    SSQ.rescaled <- uf.ssq(D,t=t,t0,A.start,B.start)
    cat("\nStress after rescaling",SSQ.rescaled/sum(D^2))
    cat("\nGradient:",rescale.cg.gradient(rescale.optim.res$par))

    cg.objective <- function(par){
      t <- par[1:ndims]
      t0 <- par[ndims + 1:ndims]
      uf.ssq(D,t,t0,A.start,B.start)
    }

    cg.gradient <- function(par){
      t <- par[1:ndims]
      t0 <- par[ndims + 1:ndims]
      t.t0 <- uf.grad(D,t,t0,A.start,B.start)
      unlist(t.t0)
    }

    optim.res <- optim(c(t,t0),
        fn=cg.objective,
        gr=cg.gradient,
        method=method,
        ...
      )

    if(optim.res$convergence>0){
      warning("Algorithm did not converge. Return code ",optim.res$convergence)
      if(length(optim.res$message))
        warning("Message: ",optim.res$message)
      }
    t <- optim.res$par[1:ndims]
    t0 <- optim.res$par[ndims + 1:ndims]

    SSQ.rigid <- uf.ssq(D,t=t,t0,A.start,B.start)
    cat("\nStress after rigid transformation",SSQ.rigid/sum(D^2))
    cat("\nGradient:",cg.gradient(optim.res$par))

    A <- sweep(sweep(A.start,2,t,"*"),2,t0,"+")
    B <- sweep(sweep(B.start,2,t,"/"),2,t0,"-")

    cat("\n")
      Delta <- sumalong(makedist(A,B)^2,1:2)
      stress <- sum((sqrt(D)-sqrt(Delta))^2)/sum(D)
      return(structure(list(
            A=A,
            B=B,
            fitted=makeD(A,B),
            dim=ndims,
            stress=stress),
        class="unfolding"))
}

makedist <- function(A,B){
  I <- nrow(A)
  J <- nrow(B)
  K <- ncol(A)
  stopifnot(K == ncol(B))

  d.ijk <- array(0,dim=c(I,J,K))

  ijk <- as.matrix(expand.grid(1:I,1:J,1:K))

  ij <-  ijk[,c(1,2)]
  ik <-  ijk[,c(1,3)]
  jk <-  ijk[,c(2,3)]

  d.ijk[ijk] <- A[ik] - B[jk]
  d.ijk
}

makeD <- function(A,B){
  I <- nrow(A)
  J <- nrow(B)
  K <- ncol(A)
  stopifnot(K == ncol(B))

  d.ijk <- array(0,dim=c(I,J,K))

  ijk <- as.matrix(expand.grid(1:I,1:J,1:K))

  ij <-  ijk[,c(1,2)]
  ik <-  ijk[,c(1,3)]
  jk <-  ijk[,c(2,3)]

  d.ijk[ijk] <- A[ik] - B[jk]
  sumalong(d.ijk^2,c(1,2))
}

full.uf.ssq <- function(D,A,B){
  Dhat <- sumalong(makedist(A,B)^2,1:2)
  R <- D - Dhat
  #browser()
  sum(R^2)
}

uf.ssq <- function(D,t,t0,A0,B0){
  A <- sweep(sweep(A0,2,t,"*"),2,t0,"+")
  B <- sweep(sweep(B0,2,t,"/"),2,t0,"-")
  Dhat <- sumalong(makedist(A,B)^2,1:2)
  R <- D - Dhat
  sum(R^2)
}


uf.grad <- function(D,t,t0,A0,B0){
  A <- sweep(sweep(A0,2,t,"*"),2,t0,"+")
  B <- sweep(sweep(B0,2,t,"/"),2,t0,"-")
  d.ijk <- makedist(A,B)
  R <- D - sumalong(d.ijk^2,1:2)
  rd.ijk <- d.ijk*c(R)
  tmp1 <- sumalong(rd.ijk,c(1,3))*c(A)
  tmp1 <- colSums(tmp1)
  tmp2 <- sumalong(rd.ijk,c(2,3))*c(B)
  tmp2 <- colSums(tmp2)/t^2
  grad.t <- -4*(tmp1+tmp2)
  grad.t0 <- -8*sumalong(rd.ijk,3)
  list(t=grad.t,t0=grad.t0)
}

trans.grad <- function(D,t,t0,A0,B0){
  A <- sweep(sweep(A0,2,t,"*"),2,t0,"+")
  B <- sweep(sweep(B0,2,t,"/"),2,t0,"-")
  d.ijk <- makedist(A,B)
  R <- D - sumalong(d.ijk^2,1:2)
  rd.ijk <- d.ijk*c(R)
  -8*sumalong(rd.ijk,3)
}

sumalong <- function(X,dims=dim(X)[-1]){
  all.dims <- seq(length(dim(X)))
  drop.dims <- setdiff(all.dims,dims)
  ans <- aperm(X,c(drop.dims,dims))
  dim(ans) <- c(prod(dim(X)[drop.dims]),dim(X)[dims])
  colSums(ans)
}

make.A.B <- function(parm,I,J,K,fixed=NULL){
  if(is.null(fixed)){
    A <- array(parm[1:(I*K)],dim=c(I,K))
    B <- array(parm[(I*K)+1:(J*K)],dim=c(J,K))
  }
  else {
    endofA <- (I*K)
    if(is.null(fixed$A))
      A <- array(parm[1:(I*K)],dim=c(I,K))
    else {
      A <- array(NA,dim=c(I,K))
      A[fixed$A$indices] <- fixed$A$start
      endofA <- endofA - length(fixed$A$start)
      A[is.na(A)] <- parm[1:endofA]
    }
    if(is.null(fixed$B))
      B <- array(parm[endofA+1:(J*K)],dim=c(J,K))
    else {
      B <- array(NA,dim=c(J,K))
      B[fixed$B$indices] <- fixed$B$start
      B[is.na(B)] <- parm[-(1:endofA)]
    }
  }
  list(A=A,B=B)
}

make.grad <- function(A.B,I,J,K,fixed=NULL){
  if(is.null(fixed)){
    return(unlist(A.B))
  }
  else {
    endofA <- (I*K)
    if(is.null(fixed$A))
      grad <- c(A.B$A)
    else {
      A.B$A[fixed$A$indices] <- NA
      grad <- A.B$A[!is.na(A.B$A)]
    }
    if(is.null(fixed$B))
      grad <- c(grad,A.B$A)
    else {
      A.B$B[fixed$B$indices] <- NA
      grad <- A.B$A[!is.na(A.B$A)]
    }
  }
  grad
}



biplot.unfolding<-function(x,dimen=c(1,2),type=attr(x,"biplot_type"),
    xlim,ylim,tpos=c(4,2),tposdim=1,asp=1,
    lty=c(1,2),
    lwd=c(1,1),
    pch=c(1,3),
    cex=c(1,1),
    col=c("black","black"),
    contour.col="black",
    contour.lty=1,
    xlab=paste("Dimension ",dimen[1]),
    ylab=paste("Dimension ",dimen[2]),
    ...){

    if(!length(type)) type <- "p"

    xdim<-dimen[1]
    ydim<-dimen[2]
    if(length(type)<2) type <- rep(type,2)
    if(length(lty)<2) lty <- rep(lty,2)
    if(length(lwd)<2) lwd <- rep(lwd,2)
    if(length(pch)<2) pch <- rep(pch,2)
    xy.range <- range(x$A[,dimen],x$B[,dimen])
    if(missing(xlim)) xlim <- xy.range
    if(missing(ylim)) ylim <- xy.range

    plot.default(rbind(x$A[,dimen],x$B[,dimen]),
            type="n",
            xlim=xlim,
            ylim=ylim,
            asp=asp,
            xlab=xlab,
            ylab=ylab,
            ...)
    abline(v=0,h=0,lty=3)
    tmatch1 <- charmatch(type[1],c("density","text","lines","both","points"))
    tmatch2 <- charmatch(type[2],c("density","text","lines","both","points"))
    if(is.na(tmatch1)) tmatch1 <- 0
    if(is.na(tmatch2)) tmatch2 <- 0

    if(tmatch1==1) {
        D <- kde2d(x=x$A[,dimen[1]],y=x$A[,dimen[2]],lims=c(xlim,ylim))
        #D <- bkde2D(x=O$A[,dimen],range.x=list(xlim,ylim),bandwidth=c(0.3,0.3))
        #names(D) <- c("x","y","z")
        contour(D,add=TRUE,lty=contour.lty,col=contour.col)
    }
    else if(tmatch1==2) {
            if(is.null(tpos)) text(x$A[,dimen],labels=rownames(x$A),col=col[1],cex=cex[1])
            else {
                points(x$A[,dimen],pch=20,col=col[1])
                tpos <- ifelse(x$A[,dimen[tposdim]]<0,tpos[1],tpos[2])
                text(x$A[,dimen],labels=rownames(x$A[,dimen]),pos=tpos,offset=0.2,col=col[1],cex=cex[1])
            }
        }
    else if(tmatch1==3) lines(x$A[,dimen],lty=lty[1],lwd=lwd[1],col=col[1])
    else if(tmatch1==4) lines(x$A[,dimen],type="b",lty=lty[1],pch=pch[1],lwd=lwd[1],col=col[1],cex=cex[1])
    else points(x$A[,dimen],pch=pch[1],col=col[1],cex=cex[1])

    if(tmatch2==1) {
        D <- kde2d(x=x$B[,dimen[1]],y=x$B[,dimen[2]],lims=c(xlim,ylim))
        contour(D,add=TRUE,lty=contour.lty,col=contour.col)
    }
    else if(tmatch2==2) {
            if(is.null(tpos)) text(x$B[,dimen],labels=rownames(x$A),col=col[2],cex=cex[2])
            else {
                points(x$B[,dimen],pch=20,col=col[2])
                tpos <- ifelse(x$B[,dimen[tposdim]]<0,tpos[1],tpos[2])
                text(x$B[,dimen],labels=rownames(x$B[,dimen]),pos=tpos,offset=0.2,col=col[2],cex=cex[2])
            }
        }
    else if(tmatch2==3) lines(x$B[,dimen],lty=lty[2],lwd=lwd[2],col=col[2])
    else if(tmatch2==4) lines(x$B[,dimen],type="b",lty=lty[2],pch=pch[2],lwd=lwd[2],cex=cex[2],col=col[2])
    else points(x$B[,dimen],pch=pch[2],cex=cex[2],col=col[2])

}

plot.unfolding<-function(x,y=NULL,dimen=1,
    discrete=attr(x,"plot_discrete"),
    use.rownames=discrete,
    xlab=paste("Dimension ",dimen),
    ...){

    if(!length(discrete)) discrete <- c(FALSE,FALSE)

    A.density <- density(x$A[,dimen])
    B.density <- density(x$B[,dimen])
    y.range <- range(
                    A.density$y,
                    B.density$y,
                    1/NROW(x$A[,dimen]),
                    1/NROW(x$B[,dimen])
                )
    x.range <- range(
                    A.density$x[!discrete[1]],
                    B.density$x[!discrete[2]],
                    x$A[,dimen],
                    x$B[,dimen]
                )
    plot(
            x=rbind(A.density$x,B.density$x),
            y=rbind(A.density$y,B.density$y),
            type="n",
            xlim=x.range,
            ylim=y.range,
            xlab=xlab,
            ylab="Density",
            ...
            )
    if(discrete[1]) {
            lines(  x=x$A[,dimen],
                    y=rep(1/NROW(x$A[,dimen]),NROW(x$A[,dimen])),
                    type="h")
            if(use.rownames[1])
                text(   x=x$A[,dimen],
                        y=rep(1/NROW(x$A[,dimen]),NROW(x$A[,dimen])),
                        labels=rownames(x$A),
                        adj=c(0,0.5),
                        srt=90
                        )
            }
    else {
            lines(  x=A.density$x,
                    y=A.density$y,
                    lty=1)
            rug(x$A[,dimen])
        }
    if(discrete[2]) {
            lines(  x=x$B[,dimen],
                    y=rep(1/NROW(x$B[,dimen]),NROW(x$B[,dimen])),
                    type="h")
            if(use.rownames[2])
                text(   x=x$B[,dimen],
                        y=rep(1/NROW(x$B[,dimen]),NROW(x$B[,dimen])),
                        labels=rownames(x$B),
                        adj=c(0,0.5),
                        srt=90
                )
            }
    else {
            lines(  x=B.density$x,
                    y=B.density$y,
                    lty=2)
            rug(x$B[,dimen])
        }
}

uapply <- function(x,FUN){
    x$A <- FUN(x$A)
    x$B <- FUN(x$B)
    x
}

print.unfolding <- function(x,...){
  cat("\n")
  cat(paste(x$dim,"dimensional",sep="-"),"unfolding solution")
  cat("\n\nGroup 1:",nrow(x$A),"points")
  cat("\nGroup 2:",nrow(x$B),"points")
  cat("\n\nStress:",x$stress)
  cat("\n")
}


relabel.unfolding <- function (x, ..., gsub = FALSE, fixed = TRUE, warn = FALSE){

  x$A <- rowrename(x$A,...,gsub=gsub,fixed=fixed,warn=warn)
  x$B <- rowrename(x$B,...,gsub=gsub,fixed=fixed,warn=warn)

  x

}