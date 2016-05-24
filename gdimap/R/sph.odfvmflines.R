##
## SPH volume processing for visualization of ODF line maps
## using vonMises clustering 
##

sph.odfvmflines <-
function(run=TRUE, fbase=NULL, savedir=tempdir(), roi=NULL,  rg=c(1,1), swap=FALSE, btoption=2, threshold=0.4, kdir=4, zfactor=5, showglyph=FALSE, showimage="linesgfa", bview="coronal", bg="white", order=4, texture=NULL, clusterthr=0.6, aniso=NULL, ...)
{
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
  showimages <- c("none", "gfa", "lines", "linesgfa", "linesrgbmap", "linesdata") ## map types
  kshow <- match(showimage, showimages)
  stopifnot(is.na(kshow) != TRUE)
  bviews <- c("sagittal", "coronal", "axial")
  kv <- match(bview, bviews)
  stopifnot(is.na(kv) != TRUE)
  ##-----------
  ## Read data
  testfilexist(fbase=fbase, btoption=btoption)
  if(btoption == 1) { ## Option 1: S2-shell (DSI 203-point 3mm)
    btable <- as.matrix(
      readtable(fbase=fbase, filename="btable.txt"))
  }
  else {
    if(btoption == 2) { ## Option 2: 3D-dsi grid 
      bval <- scantable(fbase=fbase, filename="data.bval")
      # bvec <- readtable(fbase=fbase, filename="data.bvec")
      bvec <- scantable(fbase=fbase, filename="data.bvec")
      bvec <- matrix(bvec, ncol=3)
      btable <- cbind(bval,bvec)
      rm(bval, bvec)
    }
    else stop()
  }
  ##-----------
  b0 <- which(btable[,1] == 0)
  odfvertices <- btable[-b0,2:4]
  tc <-  geometry::delaunayn(odfvertices)
  tcsurf <- t( surf.tri(odfvertices,tc))  
  ##----------------------------
  gc()
  cat("Reading data ...")
  img.nifti  <- readniidata(fbase=fbase, filename="data.nii.gz")
  volimg <- img.nifti@.Data  
  if(is.null(roi))
    mask.nifti <-
      readniidata(fbase=fbase, filename="data_brain_mask.nii.gz")
  else 
    mask.nifti <- readniidata(fbase=fbase, filename=roi)
  volmask <- mask.nifti@.Data  
  rm(img.nifti, mask.nifti)
  gc()
  ##----------------------------
  d <- dim(volmask)
  volgfa <- array(0, dim=d)              ## gfas map
  V1 <- array(0, dim=c(d, 3)) ## V1 direction
  if(is.null(rg)) {
    switch(kv,
      { nslices <- d[1]}, # sagittal,
      { nslices <- d[2]}, # coronal
      { nslices <- d[3]})  # axial
    first <- 1; last <- nslices
  }
  else { first <- rg[1]; last <- rg[2] }
  cat("\n")
  ##-----------------------------
  ## SPH process preparation
  gradient <- t(odfvertices)
  z <- design.spheven(order,gradient,lambda=0.006)
  plz <- plzero(order)/2/pi
  ngrad <- dim(gradient)[2]
  ngrad0 <- ngrad
  lord <- rep(seq(0,order,2),2*seq(0,order,2)+1)
  while(length(lord)>=ngrad0){
    order <- order-2
    lord <- rep(seq(0,order,2),2*seq(0,order,2)+1)
     cat("Reduced order of spherical harmonics to",order,"\n")
   }
  cat("Using",length(lord),"spherical harmonics\n")
  L <- -diag(lord*(lord+1)) 
  ##-----------------------------
  ## storage of 1st vector directions
  nv1 <- length(first:last)
  v1list <- vector(mode="list", nv1)
  v1count <- 0
  npar1 <- 7
  npar2 <- 15
  # rglstart()
  for (sl in (first:last)) {
    cat("slice",sl,"\n")
    #-------------------
    slicedata <- read.slice(img=volimg, mask=volmask, slice=sl,
       swap=swap, bview=bview)
    ymaskdata <- premask(slicedata)
    if(ymaskdata$empty) next # empty mask
    maxslicedata <- max(slicedata$niislicets) ##????
    S <- ymaskdata$yn[-b0,]
    S <- S / maxslicedata
    s0 <- 1
    si <- apply(S, 2, datatrans, s0)
    sicoef <- z$matrix%*% si
    sphcoef <- plz%*%L%*%sicoef
    coef0 <- sphcoef[1,]
    sphcoef[1,] <- 1/2/sqrt(pi)
    sphcoef[-1,] <- sphcoef[-1,]/8/pi
    ## odfs
    odfs <- t(z$design) %*% sphcoef
    # odfs <- apply(odfs, 2, norm01) 
     odfs <- apply(odfs, 2, anisofn, aniso=aniso) 
    ## gfas
    gfas <- apply(odfs, 2, genfa)
    gfas <- norm01(gfas) ## ??
    z2d <- ymaskdata$kin
     ## mask out thresholded values
    zx <- which(gfas <= threshold)
    if(length(zx)) {
      z2d <- z2d[-zx,]
      gfas <- gfas[-zx]
      odfs <- odfs[,-zx]
    }
    if(is.null(dim(z2d))) next
    if(length(gfas) < 8) next # 8 elements as minimum number
    lix <- dim(z2d)[1]
    switch(kv,
      { nr <- d[2]; nc <- d[3]}, # sagittal,
      { nr <- d[1]; nc <- d[3]}, # coronal
      { nr <- d[1]; nc <- d[2]})  # axial
    # nn <- 8*nr*nc
    nn <- nr*nc
    ck <- numeric(nn)
    v <- matrix(0, nrow=nn, ncol=3)
    q <- 1
    v1perslice <- matrix(0, nrow=lix,ncol=3) # store v1 directions 
    nullvectors <- NULL
    if(run) {
      ## ptm <- proc.time()
      tperc <- c(20, 40, 60, 80)
      tline <- floor(c(0.2,0.4,0.6,0.8)*lix) 
      cat("vMF estimation for ", lix, "voxels, ...\n")
      for(m in 1:lix) {
        tt <- which(tline == m)
        if(length(tt) != 0) {
          cat(paste(tperc[tt],"% ", sep="")); cflush() }
        odf <- odfs[,m]
        ## Find peaks based on clusters
        ith <- which(odf < clusterthr) 
        vx <- odfvertices[-ith,]
        n <- dim(vx)[1]
        ## Fit a vMF mixture  with k=2
        y1 <- movMF::movMF(vx, k=2, control=startctl) 
        par1 <- logb(n)*npar1
        bic1 <- 2*logLik(y1) - par1
        ## Fit a vMF mixture  with k=4
        y2 <- movMF::movMF(vx, k=4, control=startctl) 
        par2 <- logb(n)*npar2
        bic2 <- 2*logLik(y2) - par2
        if(bic1 >= bic2) { yy <- y1 }
        else { yy <- y2 }
        np <- dim(yy$theta)[1]
        # reorder by alpha weight
        yys <- sort(yy$alpha, decreasing=TRUE, index.return=TRUE)
        yyn <- yy$theta[yys$ix,]
        yy <- list(theta=yyn, alpha= yys$x)
        # normalize for visualization
        if(length(yy$alpha) < 2)
          pn <- fnorm(yy$theta)
        else
          pn <- apply(yy$theta, 1, fnorm)
        pcoords=yy$theta/pn
        pk <- list(np=np , pcoords=t(pcoords))
        if(pk$np < 2) {
          nullvectors <- c(nullvectors, m)
          next
        }
        v1perslice[m,] <- pk$pcoords[,1]
        ## optional glyph visualization
        if(showglyph) {
          if(rgl.cur() == 0) rglinit()
          else rgl.clear()
          if(pk$np > 2) {
            plotglyph(odfs[,m], odfvertices, pk, kdir=kdir)
            pp <- readline(
              "\nmore crossing-fiber glyphs  ? ('n' to exit) ") 
            if(pp == "n" ) { rgl.close(); showglyph <- FALSE; }
            else { rgl.clear( type = "shapes" ) }
          }
        }
        ## directions of max odf values to define (colored) lines
        gk <- gfas[m] 
        # pos <- c(z2d[m,],sl) # use yy-swapped mask
        pos <- c(z2d[m,],0) # use yy-swapped mask
        for(k in 1:min(pk$np, kdir)) {
          coords <- pk$pcoords
          zch <- coords[,k] * gk
          zch <- t(norm01(abs(zch)))
          if(q+1 > nn) {
            ck <- append(ck , numeric(nn))
            v <- rbind(v,  matrix(0, nrow=nn, ncol=3))
            nn <- nn+nn
          }
          ck[q] <- rgb(zch)
          ck[q+1] <- ck[q]
          pp <- pk$pcoords[,k]/2
          v[q,] <- pos
          v[q+1,]  <-  pp + pos
          q <- q+2
        }
      }
      cat("100% completed\n")
      ## print(proc.time() - ptm)
      ## save slice data
      q <- q-1;  
      v <- v[1:q,]
      ck <- ck[1:q]
      res <- list( q=q, v=v, ck=ck, nullvectors=nullvectors,
        v1perslice=v1perslice)
      fsave <- paste(savedir,"/sl",sl,".RData",sep="")
      save(res, file=fsave)
      cat("wrote", fsave,"\n")
    }
    else {
      fsave <- paste(savedir,"/sl",sl,".RData",sep="")
      load(fsave)
      cat("loaded", fsave, "\n")
      q <- res$q; v <- res$v; ck <- res$ck;
      nullvectors <- res$nullvectors;
      v1perslice <- res$v1perslice
    }
    cat("\n")
    ## remove null pk vectors
    nvl <- lix
    nnv <- length(nullvectors)
    if(nnv > 0) {
      nvl <- nvl-nnv
      v1perslice <- v1perslice[-nullvectors,]
      z2d <- z2d[-nullvectors,]
      gfas <- gfas[-nullvectors]
    }
    if(is.null(dim(z2d))) next
    v1perslicelist <- list(n=nvl, v1=v1perslice)
    v1count <- v1count+1
    v1list[[v1count]] <- v1perslicelist 
    ##---
    if(kshow != 1){
      ## one image per slice
      # rgl.init()
      if(sl == first) {
        rglstart(bg=bg)
      }
      if(kshow > 2) {
        segments3d(v, col=ck, lwd=2, alpha=1)
        rgl.viewpoint(theta=0, phi=15)
        par3d('windowRect'=c(0,0,600,600), 'zoom'=0.6, skipRedraw=FALSE)
        rgl.bringtotop()
      }
      switch(kshow,
      { ovr <- FALSE },
      { ovr <- TRUE; imgfa <- matrix(0, nr, nc); imgfa[z2d ] <- gfas },
      { ovr <- FALSE },
      { ovr <- TRUE; imgfa <- matrix(0, nr, nc); imgfa[z2d ] <- gfas },
      { ovr<- TRUE; zfactor=0.1;
        imgfa <- matrix(0, nr, nc); imgfa[z2d ] <- gfas }, # linesrgb
      { ovr <- TRUE;
      imgfa <- slicedata$niislicets[,,1] * slicedata$mask;
      imgfa <- imgfa/max(imgfa) } )
      if(ovr) {
        bg3d(col=bg) 
        light3d()  
        gfasurf3d(imgfa, zfactor=zfactor, alpha=0.6, texture=texture, ...)
        # rgl.viewpoint(theta=0, phi=0)
        rgl.viewpoint(theta=0, phi=15)
         par3d('windowRect'=c(0,0,600,600), 'zoom'=0.6, skipRedraw=FALSE)
          rgl.bringtotop()
      }
    }
    ##---
    if(sl != last){
      pp <- readline("continue to next 'rg' slice ? ('n' to exit) ") 
      if(pp == "n" ) { break; }
      else { rgl.clear( type = "shapes" ) }
    }
  }
  cat("\n")
  rm(list=ls())
  gc()
  ## save V1 directions
  # if(run) {
  #   v1file <- paste(savedir,"/V1list.RData",sep="")
  #   save(v1list, file=v1file) # list of v1 vectors 
  #   cat("saved V1 directions ",v1file,"\n")
  # }
}

