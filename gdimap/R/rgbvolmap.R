rgbvolmap <-
function(fbase=NULL, rg=c(1,1), bview="coronal",
 texture=NULL, bg="black")
{
  bviews <- c("sagittal", "coronal", "axial")
  kv <- grep(bview, bviews)
  ##----------------------------
  img.nifti  <- readniidata(fbase=fbase, filename="data_gfa.nii.gz")
  gfavol <- img.nifti@.Data
  v1.nifti  <- readniidata(fbase=fbase, filename="data_V1.nii.gz")
  V1 <- v1.nifti@.Data
  ##----------------------------
  d <- dim(gfavol)
  switch(kv,
    { nr <- d[2]; nc <- d[3]}, # sagittal,
    { nr <- d[1]; nc <- d[3]}, # coronal
    { nr <- d[1]; nc <- d[2]})  # axial
  if(is.null(rg)) {
    switch(kv,
      { nslices <- d[1]}, # sagittal,
      { nslices <- d[2]}, # coronal
      { nslices <- d[3]})  # axial
    first <- 1; last <- nslices
    rg <- c(first,last)
  }
  else { first <- rg[1]; last <- rg[2] }
  displist <- vector("list",diff(rg)+1)
  j <- 0
  for (sl in (first:last)) {
    cat(sl,"")
    j <- j+1
    #---------------
    # RGB channels
    switch(kv,
      gfaslice <- gfavol[sl,,], # sagittal
      gfaslice <- gfavol[,sl,], # coronal
      gfaslice <- gfavol[,,sl]) # axial
    ix <- which(gfaslice != 0, arr.ind=TRUE)
    gfas <- gfaslice[ix]
    tmp <- matrix(0, nrow=nr, ncol=nc)  
    switch(kv,
      vtmp <- V1[sl,,,], # sagittal
      vtmp <- V1[,sl,,], # coronal
      vtmp <- V1[,,sl,]) # axial
    v <- matrix(0,nrow=dim(ix)[1],  ncol=3)
    for(k in 1:3) {
      va <- vtmp[,,k] 
      v[,k] <- va[ix]
    }
    cvm <- gfas*abs(v) ## choose ONE color as GFA*|Vpeak| 
    zp <- array(0, dim=c(nr,nc,3))
    tmp <- matrix(0, nrow=nr, ncol=nc)  
    for(k in 1:3) {
      tmp[ix] <- cvm[,k]
      zp[,,k] <- tmp 
      tmp[ix] <- 0
    }
    if(max(zp) > 0)
      zp <- zp/max(zp)
    zpr <- as.raster(zp)
    ##----------
    ## display
    mask <- gfaslice
    zpr[ mask == 0  ] <- bg
    zpr2 <- array(0, dim=c(nc,nr))
    if(!is.null(texture)){ 
      #nm <- max(nc,nr)
      #zpr2 <- array(0, dim=c(nm,nm)) # using square texture
      zpr2[nc:1, ] <- zpr[ ,nc:1] 
      ww <- 7; hh <- 7
    }
    else {
      zpr2 <- array(0, dim=c(nc,nr))
      zpr2[1:nc, ] <- zpr[ ,nc:1]
      mm <- range(nr,nc)
      # ww <- 7; hh <- 7*mm[1]/mm[2]
      nn <- as.integer(sqrt(last-first)+1)
      w <- 10; h <- 10
      # ww <- 7*1/nn; hh <- 7*mm[1]/(mm[2] * nn)
      ww <- w*1/nn; hh <- h*mm[1]/(mm[2] * nn)
    }
    r1 <- grid::rasterGrob(zpr2, interpolate=TRUE, width=ww,
															 	height=hh, default.units="inches")
    displist[j] <- list(r1)
  }
  cat("\n")
  ## Initialize Cairo display for texture
  if(!is.null(texture)){ 
    dev.new(type="cairo", bg=bg)
    do.call(gridExtra::grid.arrange,  displist) 
    savePlot(texture)
    cat("saved",texture,"\n")
    Sys.sleep(0.25)
		dev.off()
  }
  else 
    do.call(gridExtra::grid.arrange,  displist) 
}

