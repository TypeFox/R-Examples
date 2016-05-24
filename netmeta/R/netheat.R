netheat <- function(x, random=FALSE, tau.preset=NULL, ...){ 
  
  
  if (!inherits(x, "netmeta"))
    stop("Argument 'x' must be an object of class \"netmeta\"")
  
  
  if (random==FALSE & length(tau.preset)==0){ 
    nmak <- nma.krahn(x)
    decomp <- decomp.design(x) 
    residuals <- decomp$residuals.inc.detach 
    Q.inc.design <- decomp$Q.inc.design 
  } 
  ## 
  if (length(tau.preset)==1){ 
    nmak <- nma.krahn(x, tau.preset=tau.preset) 
    decomp <- decomp.design(x, tau.preset=tau.preset) 
    residuals <- decomp$residuals.inc.detach.random.preset 
    Q.inc.design <- decomp$Q.inc.design.random.preset 
  } 
  ## 
  if (random==TRUE & length(tau.preset)==0){ 
    tau.within <- tau.within(x)
    nmak <- nma.krahn(x, tau.preset=tau.within)
    decomp <- decomp.design(x, tau.preset=tau.within) 
    residuals <- decomp$residuals.inc.detach.random.preset 
    Q.inc.design <- decomp$Q.inc.design.random.preset 
  }
  
  
  if (nmak$d<=2){
    warning("Net heat plot not available due to small number of designs: ", nmak$d)
    return(invisible(NULL))
  }
  
  
  H <- nmak$H 
  V <- nmak$V 
  design <- nmak$design 
  
  
  Q.inc.design.typ <- apply(residuals, 2, function(x) t(x) %*% solve(V)*x)
  inc <- matrix(Q.inc.design,
                nrow=nrow(Q.inc.design.typ),
                ncol=ncol(Q.inc.design.typ))
  diff <- inc - Q.inc.design.typ
  colnames(diff) <- colnames(Q.inc.design.typ)
  rownames(diff) <- rownames(residuals)
  ##
  if (all(is.na(diff))){
    warning("Net heat plot not available as no between-design heterogeneity exists.")
    return(invisible(NULL))
  }
  
  
  t1 <- -1*diff
  wi <- which(apply(t1,2,function(x) sum(is.na(x))==nrow(t1)))
  if (length(wi)>0){
    t1 <- t1[-wi, -wi, drop=FALSE]
  }
  dmat <- t1
  d1 <- dist(dmat, method="manhattan")
  d1 <- d1 + dist(t(dmat), method="manhattan")
  h1 <- hclust(d1)
  t1 <- t1[h1$order, h1$order]
  
  
  tn <- matrix(NA, nrow=nrow(t1), ncol=nrow(t1))
  
  
  for (i in 1:(nrow(t1)))
    tn[i,] <- t1[nrow(t1)-(i-1),]
  
  
  tn[tn < -8] <- -8
  tn[tn > 8] <- 8
  
  
  ncolo <- 40
  sat <- min(max(abs(tn))/8, 1)
  nblue <- ifelse((max(max(tn), 0) - min(min(tn), 0))!=0,
                  round(max(max(tn),0)/ (max(max(tn),0) - min(min(tn),0))*ncolo),
                  0)
  bf <- min(max(max(tn), 0) / abs(min(min(tn), 0)), 1)*sat
  clev <- seq(from=1, to=1-bf, len=nblue)
  heat.vec <- heat.colors(round((ncolo-nblue-1) / max(sat, 0.00001)))
  mycol <- c(heat.vec[(length(heat.vec)-(ncolo-nblue-1)+1):length(heat.vec)],
             rgb(1,1,1),rgb(clev,clev,1))
  
  
  if (length(wi)>0)
    Hp <- H[(as.character(design$comparison)[-wi]), -wi]
  else
    Hp <- H[as.character(design$comparison),]
  ##
  Hp <- Hp[h1$order,h1$order]
  
  
  Hpn <- matrix(NA, nrow=nrow(Hp), ncol=ncol(Hp))
  ##
  for (i in 1:(nrow(Hp)))
    Hpn[i,] <- Hp[ncol(Hp)-(i-1),]
  ##
  Hpn <- t(Hpn)
  
  
  va.image <- function(x, size=abs(Hpn), labels=NULL,
                       cex=1, col=mycol, bg="grey", box.bg=TRUE,
                       min.size=0.10,
                       xlab=NULL, ylab=NULL, cex.axis=1.1,...){
    
    
    if (length(wi)>0)
      design2 <- design[-wi,]
    else
      design2 <- design
    ##
    design.comb <- rep(NA, length(as.character(design2$comparison)))
    ##
    for (i in 1:(length(as.character(design2$comparison)))){
      if((((design2$narms)[h1$order])[i])>2)
        design.comb[i] <- paste(as.character(design2$comparison)[h1$order][i],
                                as.character(design2$design)[h1$order][i], sep="_")
      else
        design.comb[i] <-as.character(design2$comparison)[h1$order][i]
    }
    ##
    sc <- max(nchar(design.comb))/4
    par(oma=c(1, 1, 1, 1) + c(0, sc, sc, 0))
    ##oldpar <- par(oma=c(1, 1, 1, 1) + c(0, sc, sc, 0))
    ##on.exit(par(oldpar))
    ##
    plot(0, type="n", xlim=c(0.036, 0.963), ylim=c(0.036, 0.963),
         bty="n", xlab="", ylab="", axes=FALSE)
    ##
    if (is.null(xlab) && !is.null(rownames(x))) xlab <- rownames(x)
    if (is.null(ylab) && !is.null(colnames(x))) ylab <- colnames(x)
    ##
    rect(0,0,1,1,col=bg,border=NA)
    ##
    runit <- 1/ncol(x)
    ##
    center.x <- seq(from=runit/2, to=1-runit/2, len=ncol(x))
    center.y <- seq(from=runit/2, to=1-runit/2, len=ncol(x))
    ##
    val.grid <- seq(from=min(x), to=max(x), len=length(col)+1)
    val.grid <- val.grid[-length(val.grid)]
    if (min(x) == max(x)) val.grid <- val.grid[1]
    ##
    col[unlist(lapply(x[,1], function(arg) sum(arg < val.grid)))]
    ##
    min.sv <- min(size)
    max.sv <- max(size)

    fg.col <- bg.col <- apply(col2rgb(col)/255, 2,
                              function(arg) rgb(arg[1], arg[2], arg[3]))
    if (box.bg)
      fg.col <- apply(col2rgb(col)/255*0.8, 2,
                      function(arg) rgb(arg[1], arg[2], arg[3]))
    
    
    for (i in 1:ncol(x)){
      actual.size <- (sqrt((size[,i] - min.sv)/(max.sv - min.sv))+min.size) /
        (1+4*min.size)*runit/2
      ##
      rect(center.x-runit/2, center.y[i]-runit/2,
           center.x+runit/2, center.y[i]+runit/2,
           col=bg.col[unlist(lapply(x[,i], function(arg) sum(arg >= val.grid)))],
           border=NA)
      ##
      if (box.bg)
        rect(center.x-actual.size, center.y[i]-actual.size,
             center.x+actual.size, center.y[i]+actual.size,
             col=fg.col[unlist(lapply(x[,i], function(arg) sum(arg >= val.grid)))],
             border=NA)
      if (!is.null(labels))
        text(center.x, center.y[i], labels[,i], cex=cex)
    }
    
    
    if (!is.null(ylab))
      axis(side=2, at=center.y, labels=ylab, tick=FALSE,
           line=NA, las=2, cex.axis=cex.axis)
    ##
    if (!is.null(xlab))
      axis(side=3, at=center.x, labels=xlab, tick=FALSE,
           line=NA, las=2, cex.axis=cex.axis)

    axis(2, at=center.x, lty=0, las=2,
         labels=design.comb[length(design.comb):1], cex.axis=1)
    axis(3, at=center.y, lty=0, las=3,
         labels=design.comb, cex.axis=1)
  }
  
  
  va.image(x=t(tn)+max(abs(tn)))
  
  
  legend.col <- function(col, lev){
    opar <- par
    n <- length(col)
    bx <- par("usr")
    box.cx <- c(bx[2] + (bx[2] - bx[1]) / 1000,
                bx[2] + (bx[2] - bx[1]) / 1000 + (bx[2] - bx[1]) / 50)
    box.cy <- c(bx[3], bx[3])
    box.sy <- (bx[4] - bx[3]) / n 
    xx <- rep(box.cx, each = 2)
    par(xpd = TRUE)
    ##
    for(i in 1:n){
      yy <- c(box.cy[1] + (box.sy * (i - 1)),
              box.cy[1] + (box.sy * (i)),
              box.cy[1] + (box.sy * (i)),
              box.cy[1] + (box.sy * (i - 1)))
      polygon(xx, yy, col = col[i], border = col[i]) 
    }
    ##
    par(new = TRUE)
    plot(0, 0, type = "n",
         ylim = c(min(lev), max(lev)),
         yaxt = "n", ylab = "",
         xaxt = "n", xlab = "",
         frame.plot = FALSE)
    axis(side = 4, las = 2, tick = FALSE, line = .25)
    par <- opar
  }
  
  
  if(min(round(t1,10)) != max(round(t1,10)))
    legend.col(rev(mycol),
               seq(from=max(-max(t1),-8), to=min(-min(t1),8), len=length(mycol)))
  
  
  box(lwd=1.1)

  invisible(NULL)
}
