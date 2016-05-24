## SPH volume processing
## fslview-compatible gfa-map and V1 volumes 

sph.odfpeaks <-
function(fbase=NULL, rg=NULL, swap=FALSE, btoption=2, threshold=0.4, showglyph=FALSE, bview="coronal", savedir=tempdir(), order=4)
{
  bviews <- c("sagittal", "coronal", "axial")
  kv <- match(bview, bviews)
	stopifnot(is.na(kv) != TRUE)
  ##-----------
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
  b0 <- which(btable[,1] == 0)
  odfvertices <- matrix(btable[-b0,2:4], ncol=3)
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
  volgfa <- array(0, dim=d) ## gfas map
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

			if(length(pk$peaks) < 1) next 
      ## don't store eigenvector for cross-fiber voxels 
##      if(length(pk$peaks) < 1 | (length(pk$peaks) > 2)) {
##        nullvectors <- c(nullvectors, m)
##        next
##      }
      v1perslice[m,] <- pk$pcoords[,1]
      ## optional glyph visualization
      if(showglyph) {
        if(rgl.cur() == 0) rglinit()
        else rgl.clear()
        # if(pk$np > 2) {
          plotglyph(odfs[,m], odfvertices, pk, kdir=2, vmfglyph=FALSE)
           pp <- readline(
            "\nmore glyphs  ? ('n' to exit) ")
          if(pp == "n" ) { rgl.close(); showglyph <- FALSE; }
          else { rgl.clear( type = "shapes" ) }
        # }
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
  cat("\n")
  ##-----------------------------
  f <- paste(savedir,"/data_gfa",sep="")
  writeNIfTI(volgfa, filename=f, verbose=TRUE)
  cat("wrote",f,"\n")
  f <- paste(savedir,"/data_V1",sep="")
  writeNIfTI(V1, filename=f, verbose=TRUE)
  cat("wrote",f,"\n")
}

