
sph.odfvxgrid <-
function(fbase=NULL, rg=c(1,1), swap=FALSE, btoption=2, threshold=0.4, kdir=4, zfactor=5, showimage="glyphgfa", bview="coronal", savedir=tempdir(), bg="white", order=4, texture=NULL, ...)
{
  showimages <- c("none", "gfa", "glyph", "glyphgfa", "glyphrgbmap", "glyphdata") ## map types
  kshow <- match(showimage, showimages)
  bviews <- c("sagittal", "coronal", "axial")
  kv <- match(bview, bviews)
  stopifnot(is.na(kv) != TRUE)
  ##-----------
  ## Read data
  testfilexist(fbase=fbase, btoption=btoption)
  if(btoption == 1) { ## Option 1: S2-shell (DSI 203-point 3mm)
    btable <- as.matrix(readtable(fbase=fbase, filename="btable.txt"))
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
  ##----------------------------
  b0 <- which(btable[,1] == 0)
  odfvertices <- btable[-b0,2:4]
  tc <-  geometry::delaunayn(odfvertices)
  tcsurf <- t( surf.tri(odfvertices,tc))  
  ##----------------------------
  gc()
  cat("Reading data ...")
  img.nifti  <- readniidata(fbase=fbase, filename="data.nii.gz")
  volimg <- img.nifti@.Data  
  mask.nifti <- readniidata(fbase=fbase, filename="data_brain_mask.nii.gz")
  volmask <- mask.nifti@.Data  
  rm(img.nifti, mask.nifti)
  gc()
  ##----------------------------
  d <- dim(volmask)
  # volgfa <- array(0, dim=dim(volmask)) ## gfas map
  # V1 <- array(0, dim=c(dim(volmask), 3)) ## V1 direction
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
  ## RGB channels - normalize for rgb color
  dt2 <- dim(odfvertices[tcsurf,])[1]
  for (sl in (first:last)) {
    cat("Processing slice",sl,", ")
    slicedata <- read.slice(img=volimg, mask=volmask, slice=sl, swap=swap, bview=bview)
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
    odfs <- apply(odfs, 2, norm01) 
    ## gfas
    gfas <- apply(odfs, 2, genfa)
    z2d <- ymaskdata$kin
    ## mask out thresholded values
    zx <- which(gfas <= threshold)
    if(length(zx)) {
      z2d <- z2d[-zx,]
      gfas <- gfas[-zx]
      odfs <- odfs[,-zx]
    }
    switch(kv,
      { nr <- d[2]; nc <- d[3]}, # sagittal,
      { nr <- d[1]; nc <- d[3]}, # coronal
      { nr <- d[1]; nc <- d[2]})  # axial
    imgfa <- matrix(0, nr, nc)
    imgfa[z2d ] <- gfas 
    nn <- d[1]*d[2]
    ck <- numeric(dim(odfvertices)[1])
    dm <- dim(odfs)[2]
    grid <- matrix(0, nrow=dm*dt2, ncol=3)
    vcolors <- matrix("#FFFFFF", nrow=dm, ncol=dt2)
    lix <- dim(z2d)[1]
    tperc <- c(20, 40, 60, 80)
    tline <- floor(c(0.2,0.4,0.6,0.8)*lix) 
    cat("processing", lix,"voxels \n")
    for(m in 1:lix) {
      tt <- which(tline == m)
      if(length(tt) != 0) {
        cat(paste(tperc[tt],"% ", sep="")); cflush() }
      odf <- odfs[,m]
      gk <- gfas[m] 
      ##-------------
      ## RGB channels: calibrated colors for each odf-vertice
      zch <- odfvertices * gk
      zch <- t(apply(abs(zch), 1, norm01)) 
      ck <- rgb(zch[,1],zch[,2],zch[,3])
      ## --------
      pc <- odfvertices * odf
      pc <- pc / (2*max(pc))
      # pc <- pc / max(pc)
      pos <- c(z2d[m,],0) # use yy-swapped mask
      pcsurf <- cbind(pc[tcsurf,1] + pos[1],
                pc[tcsurf,2] + pos[2] , pc[tcsurf,3])
      b <- (m-1)*dt2; e <- b+dt2
      grid[(b+1):e, ] <- pcsurf
      vcolors[m,] <- ck[tcsurf]
    }
    cat("100% completed\n")
    ##-----------
    # ptm <- proc.time()
    if(kshow != 1){
      if(kshow != 6) rm(slicedata)
      gc() 
      if(sl == first)
        rglstart(bg=bg)
      if(kshow > 2) {
        ## Display grid of glyphs
        if(dm <= 512) {
        # Plotting once all voxels: ! slow for large grid size and glyph dim !
          cat("Plotting all voxels ...\n")
          xx <- dt2*dm
          rgl.triangles(grid[1:xx,1], grid[1:xx,2], grid[1:xx,3],
            col=t(vcolors[1:dm,]))
        } else {
          ## plot in chunks of chunk voxels
          chunk <- 512
          cat("Plotting in chunks of ", chunk," voxels ...\n")
          mc <- ceiling(dm/chunk)
          ck <- dt2*chunk
          eg <- 0; ec <- 0
          i <- 1
          while(i < mc) { 
            cat("chunk",i,"\n")
            fg <- (i-1)*ck;    eg <- fg+ck 
            fc <- (i-1)*chunk; ec <- fc+chunk
            rgl.triangles(grid[(fg+1):eg,1], grid[(fg+1):eg,2],
              grid[(fg+1):eg,3], col=t(vcolors[((fc+1):ec),]))
            i <- i+1
          }
          cat("chunk",i,"\n")
          eg2 <- dt2*dm
          ec2 <- dm
          rgl.triangles(grid[(eg+1):eg2,1], grid[(eg+1):eg2,2], grid[(eg+1):eg2,3],
            col=t(vcolors[(ec+1):ec2,]))
        }
      }
      switch(kshow,
        { ovr <- FALSE },
        { ovr <- TRUE; imgfa <- matrix(0, nr, nc); imgfa[z2d ] <- gfas },
        { ovr <- FALSE },
        { ovr <- TRUE; imgfa <- matrix(0, nr, nc); imgfa[z2d ] <- gfas },
        { ovr<- TRUE; zfactor=0.1;
          imgfa <- matrix(0, nr, nc); imgfa[z2d ] <- gfas }, # +rgb
        { ovr <- TRUE;
          imgfa <- slicedata$niislicets[,,1] * slicedata$mask;
          imgfa <- imgfa/max(imgfa) } )
      if(ovr) {
        bg3d(col=bg) 
        light3d()  
        gfasurf3d(imgfa, zfactor=zfactor, alpha=0.6, texture=texture, ...)
        rgl.viewpoint(theta=0, phi=15)
         par3d('windowRect'=c(0,0,600,600), 'zoom'=0.6, skipRedraw=FALSE)
          rgl.bringtotop()
      }
    }
    # print(proc.time() - ptm)
    ##---
    if(sl != last) {
      pp <- readline("continue to next 'rg' slice ? ('n' to exit) ") 
      if(pp == "n" ) { break; }
      else { rgl.clear( type = "shapes" ) }
    }
  }
  cat("\n")
  rm(list=ls())
  gc()
}
