##    DATforDCEMRI: a Deconvolution Analysis Tool for DCE-MRI
##    Copyright 2013 Genentech, Inc.
##
##    This package is distributed under the
##    Creative Commons Attribution-NonCommercial-ShareAlike 3.0 License
##    at http://creativecommons.org/licenses/by-nc-sa/3.0/us/
##
##    For questions or comments, please contact
##    Gregory Z. Ferl, Ph.D.
##    Genentech, Inc.
##    Development Sciences
##    1 DNA Way, Mail stop 463A
##    South San Francisco, CA, United States of America
##    E-mail: ferl.gregory@gene.com

DAT.checkData <-
function(file.name, vector.times, map.CC, mask.ROI, vector.AIF, slice.start=1, slice.stop="not.specified"){

require(grid)
require(graphics)
require(locfit)
require(matlab)
require(xtable)
require(akima)
require(R.oo)
require(R.methodsS3)
require(lattice)

if(length(dim(map.CC))==4 && slice.stop=="not.specified")
 stop("there are more than one slices of data in these files; please specify a value for slice.stop=n, so that a range of slices from 1 to n is extracted")

st <- slice.start
sp <- slice.stop

if(slice.stop=="not.specified"){
  time <- as.vector(vector.times)
  cc <- map.CC
  roi <- mask.ROI
  aif <- as.vector(vector.AIF)
}

if(slice.stop!="not.specified"){
  time <- as.vector(vector.times)
  cc <- map.CC[,,st:sp,]
  roi <- mask.ROI[,,st:sp]
  aif <- as.vector(vector.AIF)
}



cat("\n")
cat("checking dimensions of vectors and arrays...", "\n")
cat("\n")
cat("length of vector.times is", length(time), "\n")
cat("length of vector.AIF is", length(aif), "\n")

if(length(time) != length(aif))
  stop("vector.times and vector.AIF must be the same length")

if(length(dim(cc)) != (length(dim(roi))+1))
  stop("map.CC array must have n+1 dimensions when mask.ROI array has n dimensions (one of these arrays may have data corresponding to a single slice while the other has multiple slices)")

if(length(dim(cc))==4){
cat("dimensions of map.CC array are", length(cc[,1,1,1]), "x", length(cc[1,,1,1]), "x", length(cc[1,1,,1]), "slices x", length(cc[1,1,1,]), "time points", "\n")

if(length(time) != length(cc[1,1,1,]))
  stop("length of vector.times and nt dimension of the map.CC array must be equal")

if(length(aif) != length(cc[1,1,1,]))
  stop("length of vector.AIF and nt dimension of the map.CC array must be equal")

cat("dimensions of mask.ROI array are", length(roi[,1,1]), "x", length(roi[1,,1]), "x", length(roi[1,1,]), "slices", "\n")

if(length(roi[,1,1]) != length(cc[,1,1,1]))
  stop("nx dimension of the mask.ROI and map.CC arrays must be equal")

if(length(roi[1,,1]) != length(cc[1,,1,1]))
  stop("ny dimension of the mask.ROI and map.CC arrays must be equal")

if(length(roi[1,1,]) != length(cc[1,1,,1]))
  stop("number of slices within mask.ROI and map.CC arrays must be equal")
}

if(length(dim(cc))==3){
cat("dimensions of map.CC array are", length(cc[,1,1]), "x", length(cc[1,,1]), "x", length(cc[1,1,]), "time points", "\n")

if(length(time) != length(cc[1,1,]))
  stop("length of vector.times and nt dimension of the map.CC array must be equal")

if(length(aif) != length(cc[1,1,]))
  stop("length of vector.AIF and nt dimension of the map.CC array must be equal")

cat("dimensions of mask.ROI array are", length(roi[,1]), "x", length(roi[1,]), "\n")

if(length(roi[,1]) != length(cc[,1,1]))
  stop("nx dimension of the mask.ROI and map.CC arrays must be equal")

if(length(roi[1,]) != length(cc[1,,1]))
  stop("ny dimension of the mask.ROI and map.CC arrays must be equal")
}

cat("\n")
cat("...vector and array dimensions are okay.", "\n")
cat("\n")
cat("Saving data in a single R file...", "\n")

file.name <- paste(file.name, "_s", st, "-s", sp, ".RData", sep="")
dcemri.data <- list(time, cc, roi, aif, slice.start, slice.stop)
names(dcemri.data) <- c("vectorTimes", "mapCC", "maskROI", "vectorAIF", "slice.start", "slice.stop")
save(dcemri.data, file=file.name)

cat("...file saved as", file.name, "...")
cat("\n")
cat("...use the DAT() function to analyze data within this file.")
cat("\n")
cat("\n")
}

