
svyfactanal<-function(formula, design, factors,n=c("none", "sample","degf","effective","min.effective"),...){
  v<-svyvar(formula,design)
  n<-match.arg(n)
  s2<-diag(v)
  ses2<-diag(matrix(SE(v), length(s2), length(s2)))
  neff<-2*(s2/ses2)^2
  n<-switch(n, sample=nrow(design)-1, degf=degf(design), effective=1/mean(1/neff),min.effective=min(neff), none=NA)+1
  f<-factanal(covmat=v, factors=factors, n.obs=n,...)
  f$call<-sys.call()
  f
}


svyprcomp<-function (formula, design, center = TRUE, scale. = FALSE, tol = NULL, scores=FALSE,
    ...) 
{
    tms<-terms(formula)
    attr(tms,"intercept")<-0
    mf<-model.frame(formula,model.frame(design))
    naa<-attr(mf,"na.action")
    x <- model.matrix(tms,mf)
    if(length(naa))
      w<-weights(design,"sampling")[-naa]
    else
      w<-weights(design,"sampling")
    
    x<-x*sqrt(w/mean(w))
    x <- scale(x, center = center, scale = scale.)
    cen <- attr(x, "scaled:center")
    sc <- attr(x, "scaled:scale")
    if (any(sc == 0)) 
        stop("cannot rescale a constant/zero column to unit variance")
    s <- svd(x, nu = 0)
    s$d <- s$d/sqrt(max(1, nrow(x) - 1))
    if (!is.null(tol)) {
        rank <- sum(s$d > (s$d[1L] * tol))
        if (rank < ncol(x)) {
            s$v <- s$v[, 1L:rank, drop = FALSE]
            s$d <- s$d[1L:rank]
        }
    }
    dimnames(s$v) <- list(colnames(x), paste("PC", seq_len(ncol(s$v)), 
        sep = ""))
    r <- list(sdev = s$d, rotation = s$v, center = if (is.null(cen)) FALSE else cen, 
        scale = if (is.null(sc)) FALSE else sc)
    r$weights<-w/mean(w)
    
    if (scores) 
        r$x <- (x %*% s$v)/sqrt(r$weights)
    
    r$naa<-naa
    r$design<-design
    class(r) <- c("svyprcomp","prcomp")
    r
}


biplot.svyprcomp<-function(x, cols=c("black","darkred"),xlabs=NULL,weight=c("transparent","scaled","none"),
                           max.alpha=0.5,max.cex=0.5,xlim=NULL,ylim=NULL,pc.biplot=FALSE,expand=1,xlab=NULL,ylab=NULL,
                           arrow.len=0.1,
                           ...){

   if(is.null(xlabs)){
     xlabs<-1:NROW(x$x)
   } else {
     if (inherits(xlabs,"formula")){
       mf<-model.frame(xlabs,model.frame(x$design),na.action=na.pass)
       if(length(x$na.action))
         mf<-mf[-x$na.action,]
       if(ncol(mf)>1) xlabs<-sapply(mf,paste,collapse=".") else xlabs<-as.character(mf[[1]])
     }
   }

   
   scores<-x$x

   lam <- x$sdev[1:2]
   n <- NROW(scores)
   lam <- lam * sqrt(n)
   
   if (pc.biplot) 
     lam <- lam/sqrt(n)

   xx<-t(t(scores[, 1:2])/lam)
   yy<-t(t(x$rotation[,1:2]) * lam)
   

   if (missing(xlabs)) {
     xlabs <- dimnames(x)[[1L]]
     if (is.null(xlabs)) 
       xlabs <- 1L:n
   }

   xlabs <- as.character(xlabs)
   dimnames(xx) <- list(xlabs, dimnames(xx)[[2L]])
   ylabs <- dimnames(yy)[[1L]]
   
   ylabs <- as.character(ylabs)
   dimnames(yy) <- list(ylabs, dimnames(yy)[[2L]])
   
   
   weight<-match.arg(weight)

   w<-weights(x$design)
   if (length(x$na.action)) w<-w[-x$na.action]
 
   if (weight=="transparent"){
     xcexs<-par("cex")*max.cex
     rgbcol<-col2rgb(rep(cols[1],length=length(w)))
     xcols<-rgb(rgbcol[1,],rgbcol[2,],rgbcol[3,],alpha=pmax(1,255*w*max.alpha/max(w)), maxColorValue=255)
   } else if (weight=="scaled"){
     xcexs<-par("cex")*pmax(0.2, max.cex*sqrt(w/max(w)))
     rgbcol<-col2rgb(cols[1])
     xcols<-rgb(rgbcol[1,],rgbcol[2,],rgbcol[3,],alpha=max.alpha*255, maxColorValue=255)
   } else if (weight=="none"){
     rgbcol<-col2rgb(cols[1])
     xcols<-rgb(rgbcol[1,],rgbcol[2,],rgbcol[3,],alpha=max.alpha*255, maxColorValue=255)
     xcexs<-par("cex")*max.cex
   }

   
    unsigned.range <- function(x) c(-abs(min(x, na.rm = TRUE)), abs(max(x, na.rm = TRUE)))
   
    rangx1 <- unsigned.range(xx[, 1L])
    rangx2 <- unsigned.range(xx[, 2L])
    rangy1 <- unsigned.range(yy[, 1L])
    rangy2 <- unsigned.range(yy[, 2L])
   
    if (is.null(xlim) && is.null(ylim)) 
        xlim <- ylim <- rangx1 <- rangx2 <- range(rangx1, rangx2)
    else if (is.null(xlim)) 
        xlim <- rangx1
    else if (is.null(ylim)) 
        ylim <- rangx2
   
    ratio <- max(rangy1/rangx1, rangy2/rangx2)/expand
    on.exit(par(op))
    op <- par(pty = "s")

    plot(xx, type = "n", xlim = xlim, ylim = ylim, col = cols[1], 
        xlab = xlab, ylab = ylab,  ...)
    text(xx, xlabs, cex = xcexs, col = xcols, ...)
    par(new = TRUE)
    plot(yy, axes = FALSE, type = "n", xlim = xlim * ratio, ylim = ylim * 
        ratio, xlab = "", ylab = "", col = xcols, ...)
    axis(3, col = cols[2L], ...)
    axis(4, col = cols[2L], ...)
    box(col = cols[1L])
    text(yy, labels = ylabs,  col = cols[2L], ...)
   arrows(0, 0, yy[, 1L] * 0.8, yy[, 2L] * 0.8, col = cols[2L], 
          length = arrow.len)
    invisible()   
}

