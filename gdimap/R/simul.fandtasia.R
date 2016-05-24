simul.fandtasia <-
function(gdi="gqi", gridsz=32, b=4000, depth=3, sigma=0.01, clusterthr=0.6, showglyph=FALSE, savedir=tempdir(), ...)
{
  ## S2 shell grid
  s2 <- s2tessel.zorder(depth=depth)
  g0 <- s2$pc
  simul.fandtasiaSignal(g=s2$pc, gridsz=gridsz, b=b, sigma=sigma,
    savedir=savedir)
  sfield <- readniidata(fbase=savedir, filename="simfield.nii.gz")
  ## Estimate ODFs
  fielddata <- sfield@.Data
  odfs <- field.gqi(gdi=gdi, g0, fielddata, b=b, lambda=NULL,
                    savedir=savedir)
  ## with vMF clustering
  plotodfvmf(s2, odfs, clusterthr=clusterthr, showglyph=showglyph, 
             savedir=savedir, ...)
}

##----------------------

field.gqi <-
function(gdi="gqi", grad, fielddata, b=4000, lambda=NULL, savedir=tempdir())
{
  gdimethods <- c("gqi", "gqi2")
  gdimethod <- match(gdi, gdimethods)
  sz <- dim(fielddata)[4]-1
  dv <- dim(fielddata)
  odfs <- array(0, dim=c(dv[1:3], dv[4]-1))
  bn <- rep(b, dim(grad)[1])
  btable <- as.matrix(cbind(bn, grad))
  ##-----------------------------
  ## "gdimethod" process
  cat("Estimating slice odfs ...\n")
  switch(gdimethod,
      q2odf <- gqifn(odfvert=grad, btable=btable,
                     lambda=lambda),
      q2odf <- gqifn2(odfvert=grad, btable=btable,
                     lambda=lambda) )
  ##-----------------------------
  for(i in 1:dv[1]) { 
    for(j in 1:dv[2]) {
      S <- fielddata[i,j,1,]
      S <- S[2:dv[4]]
      odf <- q2odf%*%S
      odf <- odf - min(odf)
      odfs[i,j,1,] <- odf
    }
  }
  f <- file.path(savedir,"simtable.txt")
  write(t(btable), file=f, ncolumns=4)
  cat("wrote",f,"\n")
  invisible(odfs)
}

##----------------------
plotodfvmf <-
function(s2, odfsdata, clusterthr=0.6, showglyph=FALSE, savedir=tempdir(), ...)
{
  normvf <- function(x) { norm(matrix(x,length(x),1),"f") }
	## control parameters for movMF
  E <- list(...)[["E"]]
  if (is.null(E)) E <- "softmax"
  kappa <- list(...)[["kappa"]]
  if (is.null(kappa)) kappa <- "Newton_Fourier"
  minalpha <- list(...)[["minalpha"]]
  if (is.null(minalpha)) minalpha <- 8
  start <- list(...)[["start"]]
  if (is.null(start)) start <- "s"
  startctl=list(E=E, kappa=kappa, minalpha=minalpha, start=start) ## movMF inits
  ## 
  odfvertices <- s2$pc
  dm <- dim(odfsdata)
  odfs.mat <- array(odfsdata, dim=c(dm[1]*dm[2], dm[4]))
  odfs <- apply(odfs.mat, 1 , norm01)
  # gfas <- apply(odfs, 1, genfa)
  gfas <- apply(odfs, 2, genfa)
  # gfas <- norm01(gfas) ## 
  nn <- 8*dm[1]*dm[2]
  v <- matrix(0, nrow=nn, ncol=3)
  ck <- numeric(nn)
  m <- 1
  q <- 1
  npar1 <- 7
  npar2 <- 15
  volmask <- matrix(1, nrow=dm[1], ncol=dm[2])
  z2d <- which(volmask != 0, arr.ind=TRUE) 
  volgfa <- array(0, dim=c(dim(volmask),1)) ## gfas map
  V1 <- array(0, dim=c(dim(volmask),1, 3)) ## V1 direction
  V2 <- array(0, dim=c(dim(volmask),1, 3)) ## V2 direction
  lix <- dim(z2d)[1]
  v1perslice <- matrix(0, nrow=lix,ncol=3) # v1 directions 
  v2perslice <- matrix(0, nrow=lix,ncol=3) # v2 directions 
  tperc <- c(20, 40, 60, 80)
  tline <- floor(c(0.2,0.4,0.6,0.8)*lix) 
  cat("processing", lix,"voxels \n")
  for(m in 1:lix) {
    tt <- which(tline == m)
    if(length(tt) != 0) {
    cat(paste(tperc[tt],"% ", sep="")); cflush() }
    odf <- odfs[,m]
    gk <- gfas[m]
    ith <- which(odf < clusterthr) 
    ## odf <- odf - min(odf)
    vx <- odfvertices[-ith,]
    n <- dim(vx)[1]
    ##
    nc <- dim(vx)[2]
    kc <- 1:8
    npar <- nc*kc+kc-1
    bic <- -1.0e+10; nf <- 0; yy <- NULL
    for(k in seq(2,4,by=2)) {
      y2 <- movMF::movMF(vx, k=k, control=startctl) 
      par <- logb(n)*npar[k]
      bic2 <- 2*logLik(y2) - par
      if(bic2 > bic) {
        bic <- bic2
        nf <- k
        yy <- y2
      }  
    }
    np <- dim(yy$theta)[1]
    ##   pcoords <- yy$theta/max(yy$theta)
    pcoords <- yy$theta ## no scaling
    pk <- list(np=np , pcoords=t(pcoords))
    v1perslice[m,] <- pk$pcoords[,1]
    if(np == 4) {
      if(all.equal(abs(pk$pcoords[,1]), abs(pk$pcoords[,2]),
        toler=0.01) == TRUE) 
          v2perslice[m,] <- pk$pcoords[,3]
      else 
        v2perslice[m,] <- pk$pcoords[,2]
    }
    ## optional visualization of crossing-fiber glyphs
    if(showglyph) {
      if(np == 4) {
        if(rgl.cur() == 0) open3d()
        plotglyph(odf, odfvertices, pk, kdir=4)
        pp <- readline(
          "\ncontinue showing crossing-fiber glyphs ? ('n' to exit) ") 
        if(pp == "n" ) { showglyph <- FALSE; }
        else { rgl.clear( type = "shapes" ) }
      }
    }
    ## reorder 
    if(np == 4) {
      vref <- c(1,0,0)
      c1 <- crossprod(v1perslice[m,], vref)                
      c2 <- crossprod(v2perslice[m,], vref)                
      if(abs(c2) > abs(c1)) {
        tmp <- v1perslice[m,]
        v1perslice[m,] <- v2perslice[m,]
        v2perslice[m,] <- tmp
      }
    }
    ##----------------------------
    ## 1st direction of max odf values for  odfs
    pc <- odfvertices * as.vector(odf) 
    pc <- pc / (2*max(pc))
    pos <- c(z2d[m,],0)
    ## normalize for visualization
    mx <- apply(yy$theta, 1, normvf)
    coords <- t(yy$theta/mx)
    for(k in 1:min(pk$np, 4)) {
      zch <- coords[,k] * gk
        zch <- t(norm01(abs(zch)))
        ck[q] <- rgb(zch)
        ck[q+1] <- ck[q]
        pp <- coords[,k]/2
        # v[q,] <- -pp + pos
        v[q,] <- pos
        v[q+1,]  <-  pp + pos
        q <- q+2
    }
  }
  cat("100% completed\n")
  cat("\nplotting ... \n")
  ## open3d()
  segments3d(v[1:(q-1),], col=ck[1:(q-1)], lwd=2, alpha=1)
  rgl.viewpoint(0,0)
  par3d("windowRect"= c(100, 100, 700, 700), zoom=0.68)
  rgl.bringtotop()
  ##------------------
  ## storing nii volumes
  for(k in 1:3) {
    mx <- matrix(0, dm[1],dm[2])
    mx[z2d] <- v1perslice[,k]
    V1[,,1,k] <- mx
    mx <- matrix(0, dm[1],dm[2])
    mx[z2d] <- v2perslice[,k]
    V2[,,1,k] <- mx   # axial
  }
  ## gfas volume
  mx <- matrix(0, dm[1],dm[2])
  mx[z2d] <- gfas
  volgfa[,,1] <- mx # axial
  ##
  # fsave <- paste(savedir,"/vol",sl,".RData",sep="")
  # res <- list(gfa=volgfa,  v1=V1, v2=, file=fsave)
  # save(res, file=fsave)
  # cat("wrote", fsave,"\n")
  cat("\n")
  f <- file.path(savedir,"data_gfa")
  writeNIfTI(volgfa, filename=f, verbose=TRUE)
  cat("wrote",f,"\n")
  f <- file.path(savedir,"data_V1")
  writeNIfTI(V1, filename=f, verbose=TRUE)
  cat("wrote",f,"\n")
  f <- file.path(savedir,"data_V2")
  writeNIfTI(V2, filename=f, verbose=TRUE)
  cat("wrote",f,"\n")
  ##------------------
  rgl.bringtotop()
}


