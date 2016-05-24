## GQI volume processing
## fslview-compatible gfa-map and V1 volumes 

gqi.odfpeaks <-
function(gdi="gqi", fbase=NULL, rg=NULL, swap=FALSE, lambda=NULL, depth=3, btoption=2, threshold=0.4, showglyph=FALSE, bview="coronal", savedir=tempdir(), aniso=NULL)
{
  gdimethods <- c("gqi", "gqi2")
  gdimethod <- match(gdi, gdimethods)
  bviews <- c("sagittal", "coronal", "axial")
  kv <- match(bview, bviews)
  stopifnot(is.na(kv) != TRUE)
  ##---------
  ## generate S2 grid
  s2 <- s2tessel.zorder(depth=depth, viewgrid=FALSE)
  odfvertices <- s2$pc
  tcsurf <- s2$tcsurf
  ##-----------
  ## Read data
  testfilexist(fbase=fbase, btoption=btoption)
  if(btoption == 1) { ## Option 1: S2-shell (DSI 203-point 3mm)
    btable <- as.matrix(readtable(fbase=fbase, filename="btable.txt"))
  }
  else {
    if(btoption == 2) {
      ## Option 2: using a 3D-dsi grid 
      bval <- scantable(fbase=fbase, filename="data.bval")
      # bvec <- readtable(fbase=fbase, filename="data.bvec")
      bvec <- scantable(fbase=fbase, filename="data.bvec")
      bvec <- matrix(bvec, ncol=3)
      btable <- cbind(bval,bvec)
      rm(bval, bvec)
    }
    else stop()
  }
  ##--------------------
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
  ##--------------------
  d <- dim(volmask)
  volgfa <- array(0, dim=d)   ## gfas map
  V1 <- array(0, dim=c(d, 3)) ## V1 direction
  if(is.null(rg)) {
    switch(kv,
      { nslices <- d[1]}, # sagittal,
      { nslices <- d[2]}, # coronal
      { nslices <- d[3]}) # axial
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
  for (sl in (first:last)) {
    cat(sl,"")
    ## slicedata <- read.slice(img=volimg, mask=volmask, slice=sl, swap=swap)
    slicedata <- read.slice(img=volimg, mask=volmask, slice=sl,
      swap=swap, bview=bview)
    ymaskdata <- premask(slicedata)
    if(ymaskdata$empty) next # empty mask
    ##------------------
    ## odfs
    odfs <- q2odf %*% (ymaskdata$yn)
    # odfs <- apply(odfs, 2, norm01) ## normalize 
   	odfs <- apply(odfs, 2, anisofn, aniso=aniso) 
    ##------------------
    ## gfas
    gfas <- apply(odfs, 2, genfa)
    gfas <- norm01(gfas) ##??
    z2d <- ymaskdata$kin
    zx <- which(gfas <= threshold)
    if(length(zx)) {
      z2d <- z2d[-zx,]
      gfas <- gfas[-zx]
      odfs <- odfs[,-zx]
    }
    if(is.null(dim(z2d))) next
    # if(length(gfas) < 2) next # 2 elements as minimum number
    lix <- dim(z2d)[1]
    v1perslice <- matrix(0, nrow=lix,ncol=3) # store v1 directions 
    nullvectors <- NULL
    for(m in 1:lix) {
      odf <- odfs[,m]
      ##-------------------
      ## find peaks
      odf <- odf[1:(length(odf)/2)] # use half sized odf in findpeak
      pk <- findpeak(odf, t(odfvertices), tcsurf)
      ## don't store eigenvector for cross-fiber voxels 
      if(length(pk$peaks) < 1 | (length(pk$peaks) > 2)) {
        nullvectors <- c(nullvectors, m)
        next
      }
      v1perslice[m,] <- pk$pcoords[,1]
      ## optional glyph visualization
      if(showglyph) {
        if(rgl.cur() == 0) rglinit()
        else rgl.clear()
        if(pk$np < 2) { # show 1st direction only
          plotglyph(odfs[,m], odfvertices, pk, kdir=2, vmfglyph=FALSE)
           pp <- readline(
            "\nmore glyphs  ? ('n' to exit) ")
          if(pp == "n" ) { rgl.close(); showglyph <- FALSE; }
          else { rgl.clear( type = "shapes" ) }
        }
      }
    }
    # remove null pk vectors
    nvl <- lix
    nnv <- length(nullvectors)
    if(nnv > 0) {
      nvl <- nvl-nnv
      v1perslice <- v1perslice[-nullvectors,]
      z2d <- z2d[-nullvectors,]
      gfas <- gfas[-nullvectors]
    }
    ## V1 volume
    if(is.null(dim(z2d))) next
    for(k in 1:3) {
      switch(kv,
        { mx <- matrix(0, d[2],d[3])
        mx[z2d] <- v1perslice[,k]
        V1[sl,,,k] <- mx }, # sagittal
        { mx <- matrix(0, d[1],d[3])
        mx[z2d] <- v1perslice[,k]
        V1[,sl,,k] <- mx }, # coronal
        { mx <- matrix(0, d[1],d[2])
        mx[z2d] <- v1perslice[,k]
        V1[,,sl,k] <- mx } ) # axial
    }
    ## gfas volume
    switch(kv,
      { mx <- matrix(0, d[2],d[3])
      mx[z2d] <- gfas
      volgfa[sl,,] <- mx }, # sagittal
      { mx <- matrix(0, d[1],d[3])
      mx[z2d] <- gfas
      volgfa[,sl,] <- mx }, # coronal
      { mx <- matrix(0, d[1],d[2])
      mx[z2d] <- gfas
      volgfa[,,sl] <- mx } ) # axial
  }
  print(proc.time() - ptm)
  cat("\n")
  ##-----------------------------
  f <- paste(savedir,"/data_gfa",sep="")
  writeNIfTI(volgfa, filename=f, verbose=TRUE)
  cat("wrote",f,"\n")
  f <- paste(savedir,"/data_V1",sep="")
  writeNIfTI(V1, filename=f, verbose=TRUE)
  cat("wrote",f,"\n")
}

