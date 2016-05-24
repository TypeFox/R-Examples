#'@import waveslim
WISA <-
function(data, wf1 = "la8", n.levels = 11, J0 = 7, J1 = 11, boundary = "reflection", quantile = 0.9)
{
# check if data is in the acceptable format
if(is.matrix(data)==F & is.array(data)==F& is.data.frame(data)==F) stop("wrong data format: input data should be matrix or array") else {

#N<-length(data)
n.station<-dim(data)[2] #each column of the data matrix are records for each individual station
N<-dim(data)[1]         # number of observations for each station
#row-1 col-2

# if length of data >= one year, we have to
# decompose up to last level to remove the trend
if(N>=(365*1440)) n.levels<-floor(log(N,2))        #floor: result can't exceed the value of x; log(n,base) base 10\ base 2\ base exp(1)

# create the matrix to store data
preindex <- matrix(data = 0, ncol = n.station, nrow = N)

# MRA leveles used to extract periodic component
mra.sq <- matrix(data = 0, ncol = n.station, nrow = N)

# The following part constructs the part of the index, all except daily component part
#perform analysis for each line/station
for(i in 1:n.station){

### apply wavelet transform  #############
### Maximal overlap discrete wavelet transform:############
### This function performs a level J decomposition of the input vector using discrete wavelete transform.
### Boundary: Periodic-- the vector is assumed to be periodic on its defined interval.
###           Reflection-- the vector beyond its boundaries is assumed to be a symmetric reflection of itself.


data.wt <- modwt(data[, i], wf=wf1, n.levels=n.levels, boundary=boundary)

# levels in thresholding
which.levels<-1:J0

# Step II:

# threshold and change the name of the wavelet coeff. to
# emphasis that which.levels were shrinked, while others unchanged

data.wt.th<-quantile.manual.thresh.scalewise(data.wt,
which.levels=which.levels, hard=FALSE, quantile=quantile)  #From routine r script

# obtain MRA
data.wt.th.mra<-mra.wt(data.wt.th)

# prepare to reconstruct

data.recon<-numeric(N)

# add shrinked levels to data.recon
for(f in 1:J0) data.recon<-data.recon+data.wt.th.mra[[f]]

# MRA levels used to remove daily component
for(f in (J0+1):J1) mra.sq[,i] <- mra.sq[,i] + data.wt.th.mra[[f]]

# add remaining details
for(j in (J1+1):n.levels) data.recon<-data.recon+data.wt.th.mra[[j]]  #look pca r script

# compute the smooth
smooth<-data.wt.th.mra[[n.levels+1]]

# if length of data >= one year, replace smooth by its average
if(N>=(365*1440)) data.recon<-data.recon+mean(smooth)

# otherwise add unchanged smooth
if(N<(365*1440)) data.recon<-data.recon+smooth

preindex[,i]<-preindex[,i]+data.recon
}

# add deseasonalized levels, and compute Sq
# remove the daily component (pseudo-Sq)
# use the data from all stations
# vector to store pseudo-Sq component

deco<-rem.daily(mra.sq)
#data2fd  is changed to Data2fd

preindex<-preindex+deco$recon

# recon - filtered data (matrix)
# pseudo.sq - estimated Sq
# smooth - smooth
multi_return <- function() {
   my_list <- list(preindex = preindex, pseudo.sq = deco$SQ, smooth = smooth)
   return(my_list)
}
multi_return()
#return(preindex = preindex, pseudo.sq = deco$SQ, smooth = smooth)
}
}
