gqi.odfvmf <-
function(gdi="gqi", run=TRUE, fbase=NULL, savedir=tempdir(), rg=NULL, swap=FALSE,
 lambda=NULL, depth=3, btoption=2, threshold=0.4, showglyph=FALSE, bview="coronal",
 clusterthr=0.6, aniso=NULL, ...)
{
  gdimethods <- c("gqi", "gqi2")
  gdimethod <- match(gdi, gdimethods)
  bviews <- c("sagittal", "coronal", "axial")
  kv <- match(bview, bviews)
  stopifnot(is.na(kv) != TRUE)
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
  ## generate S2 grid
  s2 <- s2tessel.zorder(depth=depth, viewgrid=FALSE)
  odfvertices <- s2$pc
  tcsurf <- s2$tcsurf
  ##----------------
  ## Read data
  testfilexist(fbase=fbase, btoption=btoption)
  if(btoption == 1){ ## Option 1: S2-shell (DSI 203-point 3mm)
    btable <- as.matrix(readtable(fbase=fbase, filename="btable.txt"))
  } else {
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
  gc()
  cat("Reading data ...\n") 
  ptm <- proc.time()
  img.nifti  <- readniidata(fbase=fbase, filename="data.nii.gz")
  volimg <- img.nifti@.Data  
  mask.nifti <- readniidata(fbase=fbase, filename="data_brain_mask.nii.gz")
  volmask <- mask.nifti@.Data  
  print(proc.time() - ptm)
  rm(img.nifti, mask.nifti)
  gc()
  ##----------------
  d <- dim(volmask)
  volgfa <- array(0, dim=d)   ## gfas map
  V1 <- array(0, dim=c(d, 3)) ## V1 direction
  V2 <- array(0, dim=c(d, 3)) ## V2 direction
  V3 <- array(0, dim=c(d, 3)) ## V3 direction
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
  ## "gdimethod" process
  cat("Estimating slice odfs ...\n")
  switch(gdimethod,
    q2odf <- gqifn(odfvert=odfvertices, btable=btable,
                   lambda=lambda),
    q2odf <- gqifn2(odfvert=odfvertices, btable=btable,
                   lambda=lambda) )
  ##-----------------------------
  ## store 1st vector directions for each non-thresholded voxel 
  ## v1list: vector of lists
  nv1 <- length(first:last)
  v1list <- vector(mode="list", nv1)
  v1count <- 0
  # rglinit()
  npar1 <- 7
  npar2 <- 15
  npar3 <- 23
  for (sl in (first:last)) {
    cat("slice",sl,"\n")
    if(run) {
      slicedata <- read.slice(img=volimg, mask=volmask, slice=sl,
       swap=swap, bview=bview)
      ymaskdata <- premask(slicedata)
      if(ymaskdata$empty) next # empty mask
      ## odfs
      odfs <- q2odf %*% (ymaskdata$yn)
      odfs <- apply(odfs, 2, anisofn, aniso=aniso) 
      ## gfas
      gfas <- apply(odfs, 2, genfa)
      gfas <- norm01(gfas) ##?
      z2d <- ymaskdata$kin
      ## mask out thresholded values
      zx <- which(gfas <= threshold)
      if(length(zx)) {
        z2d <- z2d[-zx,]
        gfas <- gfas[-zx]
        odfs <- odfs[,-zx]
      }
      if(is.null(dim(z2d))) next
      if(length(gfas) < 8) next 
      lix <- dim(z2d)[1]
      v1perslice <- matrix(0, nrow=lix,ncol=3) # v1 directions 
      v2perslice <- matrix(0, nrow=lix,ncol=3) # v2 directions 
      v3perslice <- matrix(0, nrow=lix,ncol=3) # v3 directions 
      nullvectors <- NULL
      tperc <- c(20, 40, 60, 80)
      tline <- floor(c(0.2,0.4,0.6,0.8)*lix) 
      cat("processing", lix,"voxels \n")
      ptm <- proc.time()
      for(m in 1:lix) {
        tt <- which(tline == m)
        if(length(tt) != 0) {
          cat(paste(tperc[tt],"% ", sep="")); cflush() }
        odf <- odfs[,m]
        ## Find peaks based on clusters
        ith <- which(odf < clusterthr) 
        vx <- odfvertices[-ith,]
        n <- dim(vx)[1]
        lbn <- logb(n)
        ## Fit a vMF mixture  with k=2
        y1 <- movMF::movMF(vx, k=2, control=startctl)
        ## Inspect the fitted parameters:
        par1 <- lbn*npar1
        bic1 <- 2*logLik(y1) - par1
        ## Fit a vMF mixture  with k=4
        y2 <- movMF::movMF(vx, k=4, control=startctl)
        par2 <- lbn*npar2
        bic2 <- 2*logLik(y2) - par2
        if(bic1 >= bic2) { yy <- y1 }
        else {
          ## Fit a vMF mixture  with k=6
          y3 <- movMF::movMF(vx, k=6, control=startctl) 
          par3 <- lbn*npar3
          bic3 <- 2*logLik(y3) - par3
          if(bic2 >= bic3) { yy <- y2 }
          else { yy <- y3 }
        }
        np <- dim(yy$theta)[1]
        if(np < 2) {
          nullvectors <- c(nullvectors, m)
          next
        }
        # reorder by alpha weight
        yys <- sort(yy$alpha, decreasing=TRUE, index.return=TRUE)
        yyn <- yy$theta[yys$ix,]
        yy <- list(theta=yyn, alpha= yys$x)
        pk <- list(np=np , pcoords=t(yy$theta))
        v1perslice[m,] <- pk$pcoords[,1]
        if(np >= 4) 
          v2perslice[m,] <- pk$pcoords[,3]
        if(np == 6) 
          v3perslice[m,] <- pk$pcoords[,5]
        if(showglyph) {
          if(pk$np > 2) { ## show crossing-fiber glyphs only
            # normalize for visualization
            if(length(yy$alpha) < 2)
              pn <- fnorm(yy$theta)
            else
              pn <- apply(yy$theta, 1, fnorm)
            pcoords=yy$theta/pn
            pkv <- list(np=np , pcoords=t(pcoords))
            plotglyph(odfs[,m], odfvertices, pkv, kdir=6)
            points3d(vx, col="violet")
            pp <- readline(
                "\nmore crossing-fiber glyphs  ? ('n' to exit) ") 
            if(pp == "n" ) { showglyph <- FALSE; }
            else { rgl.clear( type = "shapes" ) }
          }
        }
      }
      cat("100% completed\n")
      print(proc.time() - ptm)
      ## remove null pk vectors
      nvl <- lix
      nnv <- length(nullvectors)
      if(nnv > 0) {
        nvl <- nvl-nnv
        v1perslice <- v1perslice[-nullvectors,]
        v2perslice <- v2perslice[-nullvectors,]
        v3perslice <- v3perslice[-nullvectors,]
        z2d <- z2d[-nullvectors,]
        gfas <- gfas[-nullvectors]
      }
      if(is.null(dim(z2d))) next
      ## V volumes
      for(k in 1:3) {
        switch(kv,
        { mx <- matrix(0, d[2],d[3])
        mx[z2d] <- v1perslice[,k]
        V1[sl,,,k] <- mx
        mx <- matrix(0, d[2],d[3])
        mx[z2d] <- v2perslice[,k]
        V2[sl,,,k] <- mx 
        mx <- matrix(0, d[2],d[3])
        mx[z2d] <- v3perslice[,k]
        V3[sl,,,k] <- mx }, # sagittal
        { mx <- matrix(0, d[1],d[3])
        mx[z2d] <- v1perslice[,k]
        V1[,sl,,k] <- mx
        mx <- matrix(0, d[1],d[3])
        mx[z2d] <- v2perslice[,k]
        V2[,sl,,k] <- mx
        mx <- matrix(0, d[1],d[3])
        mx[z2d] <- v3perslice[,k]
        V3[,sl,,k] <- mx }, # coronal
        { mx <- matrix(0, d[1],d[2])
        mx[z2d] <- v1perslice[,k]
        V1[,,sl,k] <- mx
        mx <- matrix(0, d[1],d[2])
        mx[z2d] <- v2perslice[,k]
        V2[,,sl,k] <- mx 
        mx <- matrix(0, d[1],d[2])
        mx[z2d] <- v3perslice[,k]
        V3[,,sl,k] <- mx } ) # axial
      }
      ## gfas volume
      fsave <- paste(savedir,"/vol",sl,".RData",sep="")
      switch(kv,
      { mx <- matrix(0, d[2],d[3])
      mx[z2d] <- gfas
      volgfa[sl,,] <- mx 
      res <- list(kv=kv, gfa=volgfa[sl,,],  v1=V1[sl,,,], v2=V2[sl,,,], v3=V3[sl,,,],
        file=fsave) }, # sagittal
      { mx <- matrix(0, d[1],d[3])
      mx[z2d] <- gfas
      volgfa[,sl,] <- mx
      res <- list(kv=kv, gfa=volgfa[,sl,],  v1=V1[,sl,,], v2=V2[,sl,,], v3=V3[,sl,,],
        file=fsave) }, # coronal
      { mx <- matrix(0, d[1],d[2])
      mx[z2d] <- gfas
      volgfa[,,sl] <- mx
      res <- list(kv=kv, gfa=volgfa[,,sl],  v1=V1[,,sl,], v2=V2[,,sl,], v2=V3[,,sl,],
        file=fsave) } ) # axial
      ##
      save(res, file=fsave)
      cat("wrote", fsave,"\n")
    } else {
      fsave <- paste(savedir,"/vol",sl,".RData",sep="")
      load(fsave)
      cat("loaded", fsave, "\n")
      switch(res$kv,
      { V1[sl,,,] <- res$v1
        V2[sl,,,] <- res$v2
        V3[sl,,,] <- res$v3
        volgfa[sl,,] <- res$gfa },
      { V1[,sl,,] <- res$v1
        V2[,sl,,] <- res$v2
        V3[,sl,,] <- res$v3
        volgfa[,sl,] <- res$gfa },
      { V1[,,sl,] <- res$v1
        V2[,,sl,] <- res$v2
        V3[,,sl,] <- res$v3
        volgfa[,,sl] <- res$gfa } )
    }
  }
  f <- paste(savedir,"/data_gfa",sep="")
  writeNIfTI(volgfa, filename=f, verbose=TRUE)
  cat("wrote",f,"\n")
  f <- paste(savedir,"/data_V1",sep="")
  writeNIfTI(V1, filename=f, verbose=TRUE)
  cat("wrote",f,"\n")
  f <- paste(savedir,"/data_V2",sep="")
  writeNIfTI(V2, filename=f, verbose=TRUE)
  cat("wrote",f,"\n")
  f <- paste(savedir,"/data_V3",sep="")
  writeNIfTI(V3, filename=f, verbose=TRUE)
  cat("wrote",f,"\n")
  ##
  V123 <- array(0, dim=c(d, 9)) ## V1 direction
  V123[,,,1:3] <- V1
  V123[,,,4:6] <- V2
  V123[,,,7:9] <- V3
  f <- paste(savedir,"/data_V123",sep="")
  writeNIfTI(V123, filename=f, verbose=TRUE)
  cat("wrote",f,"\n")
}

