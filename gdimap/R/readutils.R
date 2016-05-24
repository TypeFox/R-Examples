## $Id: read.niftivol.R,v 1.2 2012/09/15 17:48:50 arfs Exp $
##
## Default data: see DSI_Studio
## data may be prefiltered by FSL 
##   niislicets: slice time-series
##  mask: slice mask


read.slice <-
function (img, mask, slice, swap=FALSE, reorient=FALSE, bview="axial") 
{
  bviews <- c("sagittal", "coronal", "axial")
  kv <- match(bview, bviews)
  stopifnot(is.na(kv) != TRUE)
  X <- nrow(img)
  Xm <- nrow(mask)
  if (swap) { ## swap=T to be represented as in fslview
    switch(kv, ##?
    { niislicets <- img[slice, Xm:1, , ]
      mask <- mask[slice, Xm:1, ]}, 
    { niislicets <- img[Xm:1, slice, , ]  
      mask <- mask[Xm:1, slice, ]},
    { niislicets <- img[X:1, , slice, ]
      mask <- mask[Xm:1, , slice]})
  }
  else {
    switch(kv,
    { niislicets <- img[slice, , , ]
      mask <- mask[slice, , ]}, 
    { niislicets <- img[, slice , , ]  
      mask <- mask[, slice , ]},
     {  niislicets <- img[, , slice , ]
      mask <- mask[, , slice]} )
  }
  invisible(list(slice=slice, niislicets=niislicets, mask=mask, swap=swap))
}

##----------------------------------------------------------
readniidata <-
function(fbase=NULL, filename)
{
  if(is.null(fbase)) {
    file <- system.file(file.path("extdata", filename), package = "gdimap")
  }
  else {
	  file <- file.path(fbase,filename)
  }
  options(niftiAuditTrail = FALSE)
  vol.nifti <- readNIfTI(file, reorient=FALSE)
}


##--------------------------------------
#
# Mask out slice times series and keep indices
#
premask <-
function (slicedata) 
{
  slice <- slicedata$slice
  niislicets <- slicedata$niislicets
  mask <- slicedata$mask
  # kin <- which(mask == 1, arr.ind = T) # indices of pixels in mask 
  kin <- which(mask >= 1, arr.ind = T) # indices of pixels in mask 
  kin <- matrix(kin, ncol=2)
  d <- dim(kin)
  if (d[1] < 2 ) { # minimum 2 dots in mask
    return(list(empty=TRUE)) }
  yn <- matrix(0, nrow=dim(niislicets)[3], ncol=d[1])
  rx <- numeric(d[1])
  ## do not include null time series even if mask is 1 
  for (i in 1:d[1]) {
    yx <- niislicets[kin[i, 1], kin[i, 2], ]
    if (sd(yx)) 
      yn[,i] <- yx
    else 
      rx[i] <- 1
  }
  ri <- which(rx == 1)
  if(length(ri) != 0) {
    mask[kin[ri, 1], kin[ri, 2]] <- 0
    kin <- kin[-ri,]
    yn <- yn[,-ri]
  }
  # stdf <- function(y) { return((y - mean(y))/sd(y)); }
  # yn <- apply(yn, 2, stdf)
  nobs <- slicedata$nobs
  stopifnot(nobs == nrow(yn))
  nreg <- ncol(yn)
  invisible(list(yn = yn, kin = kin, nreg = nreg, empty=FALSE))
}

##----------------------------------------------------------
scantable <-
function(fbase=NULL, filename)
{
  if(is.null(fbase)) 
    file <- system.file(file.path("extdata", filename), package = "gdimap")
  else 
    file <- file.path(fbase,filename)
  invisible(scan(file))
}

##----------------------------------------------------------
readtable <-
function(fbase=NULL, filename)
{
  if(is.null(fbase)) 
    file <- system.file(file.path("extdata",filename), package = "gdimap")
  else 
    file <- file.path(fbase,filename)
  invisible(as.matrix(read.table(file)))
}


#-------------------------
testfilexist <-
function(fbase=getwd(), btoption=2)
{
  getfilename <- function(filename, fbase=getwd()) {
    if(is.null(fbase)) 
      file <- system.file(file.path("extdata",filename), package = "gdimap")
    else 
      file <- file.path(fbase,filename)
    invisible(file)
  }
  data.nii.gz <- getfilename(filename="data.nii.gz", fbase=fbase)
  data.nii <- getfilename(filename="data.nii", fbase=fbase)
  stopifnot(file.exists(data.nii.gz) | file.exists(data.nii))
  data_brain_mask.nii.gz <- getfilename(filename="data_brain_mask.nii.gz", fbase=fbase)
  data_brain_mask.nii <- getfilename(filename="data_brain_mask.nii", fbase=fbase)
  stopifnot(file.exists(data_brain_mask.nii.gz) | file.exists(data_brain_mask.nii))
  ##
  if(btoption == 1) {
  btable.txt <- getfilename(filename="btable.txt", fbase=fbase)
  stopifnot(file.exists(btable.txt))
  }
  if(btoption == 2) {
    data.bvec <- getfilename(filename="data.bvec", fbase=fbase)
    stopifnot(file.exists(data.bvec))
    data.bval <- getfilename(filename="data.bval", fbase=fbase)
    stopifnot(file.exists(data.bval))
  }
}

#-------------------------
rglstart <- 
function(bg="white")
{
  rgl.open()
  rgl.clear("all")
  # colorlut <- terrain.colors(12) # height color lookup table
  # rgl.bg(color=c(colorlut[1],"white"))
  # rgl.light()
  colorlut <- colorspace::terrain_hcl(12, h = c(0, -100), c. = c(40, 80), l = c(75, 40), power = 1)
  bg3d(col=bg)
  light3d()  
}

#-------------------------
rglinit <- 
function()
{
  rgl.open()
  rgl.clear("all")
  rgl.bg(color="white")
  rgl.light()
}

#-------------------------
cflush <-
function() 
{
  if (Sys.info()[1] == "Windows") flush.console()
  return()
}

#-------------------------
# hacked from "geometry" to elimate stack imbalace warning in R-3.* !!!

surf.tri <-
function(p,t){
  # original by Per-Olof Persson (c) 2005 for MATLAB
  # ported to R and modified for efficiency by Raoul Grasman (c) 2005

  # construct all faces
  faces = rbind(t[,-4], t[,-3], t[,-2], t[,-1]);
  node4 = rbind(t[, 4], t[, 3], t[, 2], t[, 1]);

#  #original translated from MATLAB:
#  # select the faces that occur only once --> these are the surface boundary faces
#  faces = t(apply(faces,1,sort));                     # sort each row
#  foo   = apply(faces,1,function(x) do.call("paste",as.list(x,sep=" "))); # makes a string from each row
#  vec   = table(foo);                           # tabulates the number of occurences of each string
#  ix  = sapply(names(vec[vec==1]),function(b) which(b==foo))      # obtain indices of faces with single occurence
#  tri   = faces[ix,];
#  node4 = node4[ix];
  # we wish to achieve
  #   > faces = t(apply(faces,1,sort));
  # but this is much too slow, we therefore use max.col and the fact
  # that there are only 3 columns in faces

  faces <- round(faces) ## hack
  i.max = 3*(1:nrow(faces)-1) + max.col(faces)
  i.min = 3*(1:nrow(faces)-1) + max.col(-faces)
  faces = t(faces)
  faces = cbind(faces[i.min], faces[-c(i.max,i.min)], faces[i.max])
  ix = order(faces[,1], faces[,2], faces[,3])

  # Next, we wish to detect duplicated rows in faces, that is,
  #   > qx = duplicated(faces[ix,],MARGIN=1)        # logical indicating duplicates
  # but this is also much to slow, we therefore use the fact that
  # faces[ix,] has the duplicate rows ordered beneath each other
  # and the fact that each row occurs exactly once or twice
  fo = apply(faces[ix,],2,diff)
  dup = (abs(fo) %*% rep(1,3)) == 0    # a row of only zeros indicates duplicate
  dup = c(FALSE,dup)             # first is never a duplicate
  qx = diff(dup)==0            # only zero if two consecutive elems are not duplicates
  qx = c(qx, !dup[length(dup)])      # last row is either non-duplicate or should not be selected
  tri = faces[ix[qx],]           # ix[qx] are indices of singly occuring faces
  node4 = node4[ix[qx]]

  # compute face orientations
  v1 = p[tri[,2],] - p[tri[,1],]; # edge vectors
  v2 = p[tri[,3],] - p[tri[,1],];
  v3 = p[node4,]   - p[tri[,1],];
  ix = which( apply(geometry::extprod3d(v1,v2) * v3, 1, sum) > 0 )
  tri[ix,c(2,3)] = tri[ix,c(3,2)]
  rownames(tri) = NULL
  tri
}

