# mQTL package: Adapted by Jean-Baptiste Cazier and Lyamine Hedjazi

if(getRversion() >= "2.15.1")  utils::globalVariables(c("recursion","peakParam","MAX_DIST_FACTOR","MIN_RC", 
"maxiL","ppm","res","best","permo","peaksValidated","ref","p","top","maxi","Xuv","Uvscaled","indicesinf","n.mar","met","extendedSegment"))

###############################################################Code RSPA###################################################################"
align_mQTL<-function(datafile,outdat )
{
# Recursive Segment-Wise Peak Alignment (RSPA) for accounting peak position
# variation across multiple 1H NMR biological spectra
# Output: Aligned spectra
# Author, L. Hedjazi, ICAN, Paris 2013

data<-read.csv(datafile,header=T)
np=dim(data)[2]
Sp=data[,-c(np-2,np-1,np)]
ppm=as.numeric(gsub("ppm_","",colnames(Sp)))

setupRSPA(ppm)

if (!is.matrix(Sp)){Sp<-as.matrix(Sp)}

obs<-dim(Sp)[1]
dim<-dim(Sp)[2]

donormalise<-TRUE

#Quotient probabilistic normalisation
normD<-normalise(abs(Sp),'prob')

print('...Automatic selection of a reference spectrum...')
idx<-selectRefSp(normD$Sp,recursion$step)
refSp<-normD$Sp[idx,]


#segmentate a reference spectrum
refSegments<- segmentateSp(refSp, peakParam)

for (index in 1:obs)
{
    # segmentate a test spectrum 
    testSegments<-segmentateSp(normD$Sp[index,], peakParam)

    # match test and reference segments
    attachedSegs<-attachSegments(refSegments,testSegments)

    refSegments<-attachedSegs$refSegmentsNew
    testSegments<-attachedSegs$testSegmentsNew

    Segs<-matchSegments(refSp,normD$Sp[index,], testSegments,refSegments,MAX_DIST_FACTOR, MIN_RC)

    # align a test spectrum 

Sp[index,]<- alignSp(refSp,Segs$refSegs,normD$Sp[index,],Segs$testSegs,recursion,MAX_DIST_FACTOR, MIN_RC)


    #undo normalisation
    if (donormalise==FALSE)
    {
       Sp[index,]<-Sp[index,] * normD$factors[index]
    }
}

if (donormalise==TRUE)
{
    print('...')
    print('...Returning normalised and aligned spectra...')
}else{
    print('')
    print('...Returning aligned spectra...')
}

data[,-c(np-2,np-1,np)]<-Sp
write.table(data, outdat, quote=FALSE, row.names=FALSE,col.names=TRUE,sep=",")

}

setupRSPA<-function(ppm){

# setup of alignment parameters
# input: ppm - chemical shift scale
# L. Hedjazi, ICAN, 2013

# Configuration of the algorithm invariant parameters 

configureRSPA(ppm)

### These parameters can be changed to improve the algorithm performance
#################### Recursive alignment parameters ########################

recursion$resamblance<-0.95 # Stop criterium of the recursion indicating
#the complete alignment of segment peaks
# default 98#

######################## Segmentation parameters############################
peakParam$ppmDist <- 0.03# (ppm)  #distance to concatenate adjacent peaks 
#to represent a multiplet in a single segment ### default 0.03# 
peakParam$ampThr <- 0.3 # amplitude value to threshold small peaks # 
#[] - automatic determination based on the 5# of the most intensive peaks

############Prevention of Misalignment on a local scale##############
recursion$segShift<-0.02#(ppm)  max peak shift for large peaks
recursion$inbetweenShift<-0.02 #(ppm) max shift for small peaks
recursion$acceptance<-0.5 # if resemblance after the alignment between modified test 
#and reference segments is less than the acceptance value, the alignment is not accepted. 
#Higher value is more stringent, and align only highly resemble peaks.

if (!missing(ppm)){
    peakParam$ppmDist<-ppmToPt(peakParam$ppmDist,0,ppm[2]-ppm[1])
    recursion$segShift<-ppmToPt(recursion$segShift,0,ppm[2]-ppm[1])
    recursion$inbetweenShift<-ppmToPt(recursion$inbetweenShift,0,ppm[2]-ppm[1])}

assign("peakParam", peakParam,envir = parent.frame())
assign("recursion", recursion,envir = parent.frame())
assign("MAX_DIST_FACTOR", MAX_DIST_FACTOR,envir = parent.frame())
assign("MIN_RC", MIN_RC,envir = parent.frame())
}



configureRSPA<-function(ppm){

recursion<-list()
peakParam<-list()

################Recursion minimum segment size##################
recursion$minSegWidth<-0.01 #(ppm) Stop criteria of the recursion - the size of the smallest peak
recursion$step<-0.02 #(ppm) - used for calculation of variance-scaled CC
#################### Peak peaking parameters ###############################
peakParam$minPeakWidth <- 0.005 #min peak width in ppm scale
peakParam$iFrameLen<-0.005 #Savitzky-Golay frame length in ppm scale
peakParam$iOrder<-3 #polynomial order of Savitzky - Golay filter
peakParam$peakEdgeMax<-0.2 
########### Matching parameters ############################################
MAX_DIST_FACTOR<-0.5 # The distance matching parameter (0.5*peak_width)
MIN_RC<-0.25 # Minimum resamblance coefficient

if (!missing(ppm)) {
    peakParam$minPeakWidth<-ppmToPt(peakParam$minPeakWidth,0,ppm[2]-ppm[1])
    peakParam$iFrameLen<-ppmToPt(peakParam$iFrameLen,0,ppm[2]-ppm[1])
    recursion$minSegWidth<-ppmToPt(recursion$minSegWidth,0,ppm[2]-ppm[1])
    recursion$step<-ppmToPt(recursion$step,0,ppm[2]-ppm[1])}

assign("peakParam", peakParam,envir= parent.frame())
assign("recursion", recursion,envir= parent.frame())
assign("MAX_DIST_FACTOR", MAX_DIST_FACTOR,envir = parent.frame())
assign("MIN_RC", MIN_RC,envir = parent.frame())
}




ppmToPt<- function(ppmValues, firstPtPpm, resolution){

if(nargs() < 2 | missing(firstPtPpm))
    {stop('nargs < 2 || missing(firstPtPpm)')}


if(nargs() < 3 | missing(resolution))
    {resolution <- ppmValues[2] - ppmValues[1]}

if(length(firstPtPpm)>1)
    {stop(paste('First ppm should be a number, got non-scalar value:', firstPtPpm))}

if(length(resolution)>1)
    {stop(paste('Resolution ppm should be a number, got non-scalar value:', resolution))}

ppmShift <- ppmValues - firstPtPpm

pt <- round(ppmShift / resolution) + 1


return(pt)

}


normalise<-function(X,method)
{
# Removing dilutions between biofluid samples (normalisation of spectra)

# Reference: "Probabilistic quotient normalization as robust method...
# to account for dilution of complex biological mixtures". 

# X [observation dimension]
# method - total area (method<-"total")
# or qoutient probabistic method (method<-"prob")
# L. Hedjazi, ICAN, 2013

if (nargs()<2)
{
    stop('incorrect number of input parameters')
}

obs<-dim(X)[1]
dimm<-dim(X)[2]


# Preliminary normalisation to a total area
factors<-rep(NaN,obs)
for (i in 1:obs) {
    factors[i]<-sum(X[i,])
    X[i,]<-X[i,]/factors[i]}

switch(method,
       total=return(),
       prob={X[0==X]<-0.00000001
           if (nargs()<3) {
            #normRef<-median(X)
            normRef<-X[1,]}

        F<-X/(matrix(rep(normRef, each=obs), ncol=length(normRef)))

        for (i in 1:obs) {
            X[i,]<-10000*X[i,]/median(F[i,])
            factors[i]<-(factors[i]*median(F[i,]))/10000}
        })

return(list(Sp=X, factors=factors, NormSp=normRef))
}

selectRefSp<-function(X,step)
{
# Automated selection of a reference spectrum based on the highest similarity to
# all other spectra
# Input: X - spectra
#        step - used to scale spectral regions down to specific bin_size

obs<-dim(X)[1]
splength<-dim(X)[2]

if (step >= splength)
{
    CC<-cor(p,ref)
    CC<-CC[1,2]
    return()

}

bin_count<-ceiling(splength/step)
bin_width<-ceiling(splength/bin_count)
bins<-seq(1,splength, bin_width)

if (bins[length(bins)]!=splength) {
    bins<-c(bins,splength)
    bin_count<-bin_count+1
}

for (i in 1:(length(bins)-1))
{
    istart<-bins[i]
    iend<-bins[i+1]-1
    seglength<-iend-istart+1
    X[,istart:iend]<-X[,istart:iend]-apply(t(X[,istart:iend]),2,mean)*rep(1,seglength)
    stdX<-apply(t(X[,istart:iend]),1,sd)
    stdX[stdX==0]<-1
    X[,istart:iend]<-X[,istart:iend]/(stdX*rep(1,seglength))
}

CC<-abs(cor(t(X)))
index<-which(apply(CC,2,prod)==max(apply(CC,2,prod)))
return(index[1])
}

sgolay<-function(k,F,W)
{
# Check if the input arguments are valid
if (round(F) != F) {stop('Frame length must be an integer.')}
if (F%%2!= 1) {stop('Frame length must be odd.')}
if (round(k)!= k) {stop('Polynomial degree must be an integer.')}
if (k > F-1) {stop('The degree must be less than the frame length.')}
if (nargs() < 3) {
   # No weighting matrix, make W an identity
   W <- diag(F)
}else{
   #Check for right length of W
   if (length(W) != F) { stop('The weight vector must be of the same length as the frame length.')}
   #Check to see if all elements are positive
   if (min(W) <= 0) {stop('All the elements of the weight vector must be greater than zero.')}
   #Diagonalize the vector to form the weighting matrix
   W <- diag(W)
}

#Compute the projection matrix B

S<-outer((-(F-1)/2):((F-1)/2), 0:k, FUN="^") # Compute the Vandermonde matrix  

qrm<- qr(sqrt(W)%*%S)

R<-upper.tri(qrm$qr,diag=TRUE)*qrm$qr

G<-S %*% pinv(R)%*% t(pinv(R))# Find the matrix of differentiators

B <- G%*%t(S)%*%W

return(G)
}

pinv <- function (A)
{
    s <- svd(A)
    # D <- diag(s$d) Dinv <- diag(1/s$d)
    # U <- s$u V <- s$v
    # A <- U D V'
    # X <- V Dinv U'
    s$v %*% diag(1/s$d) %*% t(s$u)
}


sgolayDeriv <- function(dpSpectr,iOrder,iFrameLen,j) 
{
# Calculate smoothed derivates using Savitzky - Golay filter
# iFrameLen- the length of frame window

if (nargs()<1) {stop('Incorrect number of input arguments')}

if (nargs()<2) {iOrder <- 3}

if (nargs()<3) {iFramLen<-11}

if (nargs()<4) {j<-2} #Derivative

iFrameLen<-(floor(iFrameLen/2))*2+1 # iFramLen must be odd

iSpecLen <- length(dpSpectr)

g<- sgolay(iOrder,iFrameLen)

#dpDerivs[1:iFrameLen]<- 0
dpDerivs<-as.vector(rep(0,iFrameLen))
dpDerivs[(iSpecLen-((iFrameLen+1)/2)):iSpecLen]<-0

for (n in ((iFrameLen+1)/2):(iSpecLen-(iFrameLen+1)/2)){
    #calculate first order derivate
    dpDerivs[n]<-(t(g[,j]) %*%(dpSpectr[(n - ((iFrameLen+1)/2)+ 1): (n + ((iFrameLen+1)/2) - 1)]))   
    }
 return(dpDerivs)
}

peakPeaks<-function(SpSmooth,dpDerivs,Sp)
{
# Peak peaking:
# Input: SpSmooth - smoothed spectrum
#        dpDerivs - smoothed derivative of the spectrum

# - the peak is identified if derivative crosses zero,
# i.e. sign(X'(i))>sing(X'(i+1))

# Author Lyamine Hedjazi, ICAN 2013


iSpecLen<-length(SpSmooth)
iPeakInd <-1

# Pre-location
peaks<-list()
#peaks$maxPos<-rep(NaN,iSpecLen)

for (i in (1:(iSpecLen-1)))
{
    # coarse peak maximum position
    if ((dpDerivs[i]>=0)&& (dpDerivs[i+1]<0)) 
     {
        peaks$maxPos[iPeakInd]<-i+1
        # Temporary starting and ending peak positions
        iPeakInd<-iPeakInd+1
    }
}

peakCount<-iPeakInd-1
peaks$maxPos<-peaks$maxPos[1:peakCount]

targetPkIdx<-1

for (srcPkIdx in 1:peakCount)
{
    maxPos<-peaks$maxPos[srcPkIdx]
    
    while ((maxPos > 2) && (maxPos < (iSpecLen-2)))
     {
        if (SpSmooth[maxPos-1]<=SpSmooth[maxPos]&&SpSmooth[maxPos]>=SpSmooth[maxPos+1]) 
        {
            
           if (targetPkIdx > 1 && peaks$maxPos[targetPkIdx-1]==maxPos) 
           {
                # the same maximum value - just skip it
                break
           }
            # save the new index:
            peaks$maxPos[targetPkIdx]<- maxPos
            targetPkIdx <- targetPkIdx + 1
            break
         }
        if (SpSmooth[maxPos]<=SpSmooth[maxPos+1])
         {
            maxPos<-maxPos+1
        }else { 
            if(SpSmooth[maxPos]<=SpSmooth[maxPos-1])
            {maxPos<-maxPos-1}
        }
    }
}

peakCount<-targetPkIdx-1
peaks$maxPos<-peaks$maxPos[1:peakCount]


for (i in 1:peakCount)
  {
    j<-peaks$maxPos[i]
    k<-peaks$maxPos[i]

    # left boundary
    while ((SpSmooth[j]>=SpSmooth[j-1]) && (j-1!=1)) #first index
        {j<-j-1}

    # right boundary
    while ((SpSmooth[k]>=SpSmooth[k+1]) && (k+1 != iSpecLen)) #last index
        {k<-k+1}


    peaks$startPos[i]<-j
    peaks$endPos[i]<-k
    peaks$centre[i]<-ceiling((k+j)/2)
    peaks$startVal[i]<-SpSmooth[j]
    peaks$endVal[i]<-SpSmooth[k]
    peaks$index[i]<-i
    #Use peak maximum position from original spectrum
    #instead of smoothed one.
    peaks$maxVal[i]<-max(Sp[j:k])
    maxInd<-which.max(Sp[j:k])
    peaks$maxPos[i]<-j+maxInd-1

    #estimate the baseline as minimum value:
    peaks$basl[i] <- min(cbind(SpSmooth[k], SpSmooth[j]))
}

peaks<-as.data.frame(peaks)
return(peaks)
}


validatePeaks<-function(SpSmooth,peaks,peakParam) 
{
# input:          Peak peaking details
#                 peaks:  maxPos - peak maxium position
#                         startPos - start position
#                         endPos - end position
#                         maxVal - maximum value
#                         startVal - start value
#                         endVal - end value
#                         basl - baseline value
#                         index - peak index
########################################################################
#             Peak validation parameters
#             peakParam:  minPeakWidth - minimum peak width
#             ampThr - amplitude threshold automatically determined if it
#             is zero
# Author, L. Hedjazi, ICAN 2013

peakCount<-nrow(peaks)
# Matrix pre-location
validatedPeaks<-peaks
minPeakWidth<-peakParam$minPeakWidth

if (("ampThr" %in% names(peakParam))==FALSE|| peakParam$ampThr==FALSE){
    ampThr<-getAmpThr(peaks)
    peakParam$ampThr<-ampThr
}else{
    ampThr<-peakParam$ampThr}

if ((peakParam$ampThr>(0.9*max(SpSmooth)))||(peakParam$ampThr<(1.1*min(SpSmooth)))) {
    stop('Peak validation threshold exceeds spectrum maximum and minimum values')}

index<-1

for (i in 1:peakCount) 
{
    if (((peaks$endPos[i]-peaks$startPos[i]) > minPeakWidth) && (peaks$maxVal[i]-peaks$basl[i] > ampThr)) {
 
        validatedPeaks[index,]<-peaks[i,]
        index<-index+1}
}

if (index> 1)
 {
    PeakCount<-index-1
        validatedPeaks<-validatedPeaks[1:PeakCount,]
}else{
    stop('wrong peak peaking parameters: No Validated peaks')
}


minsegwidth<-1e10
for (i in 1:PeakCount)
{
    startPos<-validatedPeaks$startPos[i]
    maxPos<-validatedPeaks$maxPos[i]
    endPos<-validatedPeaks$endPos[i]
    segwidth<-endPos-startPos
    # Determine the peak boundaries
    edgeVal<-validatedPeaks$maxVal[i]*peakParam$peakEdgeMax

    tstleft<-which((SpSmooth[startPos:maxPos]-validatedPeaks$basl[i])>= edgeVal)
     if (length(tstleft)==0)
      {
        validatedPeaks$LeftEdge[i]<-startPos
    }else{
        LeftEdge<-which((SpSmooth[startPos:maxPos]-validatedPeaks$basl[i])>= edgeVal)[[1]]
        validatedPeaks$LeftEdge[i]<-startPos+LeftEdge-1
    }
    
    tstright<-which(SpSmooth[maxPos:endPos]-validatedPeaks$basl[i]>= edgeVal)
    
    if (length(tstright)==0)
     {
        validatedPeaks$RightEdge[i]<-endPos
    }else{
        RightEdge<-tail(which(SpSmooth[maxPos:endPos]-validatedPeaks$basl[i]>= edgeVal),n=1)
        validatedPeaks$RightEdge[i]<-maxPos+RightEdge-1
    }

    if (minsegwidth>segwidth) 
    {minsegwidth<-segwidth}
}

assign("peakParam", peakParam,envir= parent.frame())
assign("peaksValidated",validatedPeaks,envir=parent.frame())
assign("minsegwidth", minsegwidth,envir= parent.frame())
}


getAmpThr<-function(peaks)
{
# Automatic determination of amplitude threshold for peak peaking
# based on the 5% of the most intensive peaks 
# Author, L. Hedjazi, ICAN Paris, 2013

PeakCount<-nrow(peaks)
peakMaxValues <- rep(NaN,PeakCount)

for (i in 1:PeakCount)
{
      peakMaxValues[i]<-peaks$maxVal[i]-peaks$basl[i]
}
### Select threshold based on 5% of the most intensive peaks

index<-floor(PeakCount*0.95)
peakSortedValuess<-sort(peakMaxValues)
ampThr<-peakSortedValuess[index]
return(ampThr)

}


segmentate<-function(Sp, peaks, peakParam)
{
# Combination of adjacent peaks into larger segments
#  Input: SpSmooth - spectrum of interest
#
#                 peaks.  maxPos - peak maxium position
#                         startPos - starting position
#                         endPos - ending position
#                         maxVal - maximum value
#                         startVal - starting value
#                         endVal - ending value
#                         basl - baseline value
#                         index - peak index
#
#                 peakParam.ppmDist - distance to combine adjacent peaks
#
# Author: Lyamine Hedjazi, ICAN Paris, 2013

peakCount<-nrow(peaks)
segments<-list()
ppmDist<-peakParam$ppmDist

segmentIndex<-1
peakIndex<-1

while (peakIndex<=peakCount)
{

    segments$start[segmentIndex]<-peaks$startPos[peakIndex]
    segments$PeakLeftBoundary[segmentIndex]<-as.data.frame(peaks$LeftEdge[peakIndex])
    segments$PeakRightBoundary[segmentIndex]<-as.data.frame(peaks$RightEdge[peakIndex])
    #segments$Peaks<-rbind(segments$Peaks,peaks[peakIndex,])

    segments$Peaks[[segmentIndex]]<-peaks[peakIndex,]
    
    while (peakIndex<=peakCount)
     {
         # check whether the next peak is part of the same segment
        #TODO: optimise no matter to store PeakLeft(Right)Boundary if we
        #store segment peaks themeselves.
        includePeak <- (peakIndex<peakCount) && ((peaks$maxPos[peakIndex+1]-peaks$maxPos[peakIndex])<ppmDist)

        if (includePeak)
        {
          peakIndex<-peakIndex+1
          segments$PeakLeftBoundary[[segmentIndex]]<-c(segments$PeakLeftBoundary[[segmentIndex]], peaks$LeftEdge[peakIndex])
          segments$PeakRightBoundary[[segmentIndex]]<-c(segments$PeakRightBoundary[[segmentIndex]],peaks$RightEdge[peakIndex])
          segments$Peaks[[segmentIndex]]<-rbind(segments$Peaks[[segmentIndex]], peaks[peakIndex,])
          #segments$end[segmentIndex]<-NULL
        }else {
            segments$end[segmentIndex]<-peaks$endPos[peakIndex]
            segments$centre[segmentIndex]<-ceiling((segments$start[segmentIndex]+segments$end[segmentIndex])/2)
            segmentIndex<-segmentIndex+1
            peakIndex<-peakIndex+1
            break
        }
    }
}

segmentCount<-segmentIndex-1

segments$PeakLeftBoundary <-segments$PeakLeftBoundary[1:segmentCount]
segments$PeakRightBoundary<-segments$PeakRightBoundary[1:segmentCount]
segments$Peaks<-segments$Peaks[1:segmentCount]
segments$end<-segments$end[1:segmentCount]
segments$centre<-segments$centre[1:segmentCount]

segmentVld<-segments

 for (i in 1:segmentCount)
  {
        segmentVld$start[i]<-min(segmentVld$start[i], segmentVld$PeakLeftBoundary[[i]])
        segmentVld$end[i]<-max(segmentVld$end[i], segmentVld$PeakRightBoundary[[i]])
  }

return(segmentVld)

}

segmentateSp<-function(Sp,peakParam)
{
# Determination of highly intensive peaks in the spectrum of interest and 
# subsequent concatenation of closely located peaks into 
# larger segments 
# Algorithm: - smooth spectrum X using SG
#            - locate peak maxima when the first 
#              derivative crosses zero, i.e. sign(X'(i))>sing(X'(i+1))
#            - validate peaks (/or eliminate noisy peaks)
#            - concatenate closely located peaks into larger segments
# Input: 
#                             Sp - spectrum
# Peak parameters:  peakParam$                             
#                             ampThr - amplitude threshold   [default 2*median(peaksMaxValues)] 
#                             iFrameLen - Savitzky-Golay frame length
#                             iOrder - polynomial order of Savitzky - Golay filter
#                             minPeakWidth - min peak size
#                             ppmDist - distance to concatenate adjacent peaks
#Output:                     segments
#Author: L. Hedjazi, ICAN Paris 2013 

#perform Savitzkiy Golay smoothing
SpDerivs <- sgolayDeriv(Sp,peakParam$iOrder,peakParam$iFrameLen,2)
SpSmooth <- sgolayDeriv(Sp,peakParam$iOrder,peakParam$iFrameLen,1)

#indentify peaks
peaks<-peakPeaks(SpSmooth,SpDerivs,Sp)
#validate peaks
validatePeaks(SpSmooth,peaks,peakParam)
#locate segments
segments<-segmentate(Sp,peaksValidated,peakParam)

assign("peakParam", peakParam,envir = parent.frame())
return(segments)
}


attachSegments<-function(refSegments,testSegments)
{
#Concatenation of test and reference segments to their ensure one-to-one
#correspondence
#Algorithm
#    - For each reference segment within segment boundaries, i.e. between
#      initial and final positions, find all centre (middle) positions of test segments 
#      and merge those segments, if more than one centre position is found
#    - Apply the same procedure for each test segment
#    Input: refSegments
#           testSegments

testSegmentsNew<-validateSegments(refSegments,testSegments)
refSegmentsNew<-validateSegments(testSegmentsNew,refSegments)

return(list(testSegmentsNew=testSegmentsNew,refSegmentsNew=refSegmentsNew))
}


validateSegments<-function(refSegments,testSegments)
{
# Combine segments 
ref_length<-length(refSegments$Peaks)
int_length<-length(testSegments$Peaks)

i_centre<- get_central_pos(testSegments)
index<-1
segments<-list()
concatination<-list()
concatination$segments<-list()

for (i in 1:ref_length)
{
    left_bnd<-refSegments$start[i]
    right_bnd<-refSegments$end[i]
    ind<-which((i_centre>left_bnd) & (i_centre<right_bnd))

    if (length(ind)>1)
     {
        concatination$segments[[index]]<-ind
        index<-index+1
    }
}

if (index==1) 
{
segments<-testSegments
}else{

    conc_index<-1
    seg_index<-1
    i<-1
    while (i<=int_length)
     {
        if (conc_index<index)
         {
            if (concatination$segments[[conc_index]][1]==i)
             {
                indexes<-concatination$segments[[conc_index]]

                tstSg<-list()
                tstSg$start<-testSegments$start[indexes]
                tstSg$PeakLeftBoundary<-testSegments$PeakLeftBoundary[indexes]
                tstSg$PeakRightBoundary<-testSegments$PeakRightBoundary[indexes]
                tstSg$Peaks<-testSegments$Peaks[indexes]
                tstSg$end<-testSegments$end[indexes]
                tstSg$centre<-testSegments$centre[indexes]


                segment<-unite_segments(tstSg)

                #segments(seg_index)<-segment

                segments$start[seg_index]<-segment$start
                segments$PeakLeftBoundary[[seg_index]]<-segment$PeakLeftBoundary
                segments$PeakRightBoundary[[seg_index]]<-segment$PeakRightBoundary
                segments$Peaks[[seg_index]]<-segment$Peaks
                segments$end[seg_index]<-segment$end
                segments$centre[seg_index]<-segment$centre

                i<-tail(concatination$segments[[conc_index]],n=1)+1
                seg_index<-seg_index+1
                conc_index<-conc_index+1
                next
            }
       }
        #segments[seg_index]<-testSegments[i]
 
                segments$start[seg_index]<-testSegments$start[i]
                segments$PeakLeftBoundary[seg_index]<-testSegments$PeakLeftBoundary[i]
                segments$PeakRightBoundary[seg_index]<-testSegments$PeakRightBoundary[i]
                segments$Peaks[seg_index]<-testSegments$Peaks[i]
                segments$end[seg_index]<-testSegments$end[i]
                segments$centre[seg_index]<-testSegments$centre[i]

        i<-i+1
        seg_index<-seg_index+1
    }
}
return(segments)
}



get_central_pos<-function(segments)
{
# segments 
ilength<-length(segments$start)
i_centre<-segments$centre[1:ilength]

return(i_centre)
}

unite_segments<-function(segments)
{
# concatination of segments
ilength<-length(segments$start)
i_centre<-NULL
segment<-list()

segment$start<-segments$start[1]

for (i in 1:ilength)
{
    segment$PeakLeftBoundary<-c(segment$PeakLeftBoundary, segments$PeakLeftBoundary[[i]])
    segment$PeakRightBoundary<-c(segment$PeakRightBoundary,segments$PeakRightBoundary[[i]])
    segment$Peaks<-rbind(segment$Peaks,segments$Peaks[[i]])
}

segment$end<-tail(segments$end,n=1)

segment$centre<-ceiling((segment$start+segment$end)/2)
return(segment)
}



matchSegments<-function(refSp,intSp,intSegments,refSegments,MAX_DIST_FACTOR, MIN_RC)
{
# Matching of the segment of interest to the corresponding reference using
# Fuzzy logic approach
# Algorithm: - take segment of interest
#            - take reference segments
#            - calculate relative distance between them
#            - calculate relative resamblance between them
#            - find min value of relative distance and resamblance 
#            - use it as representative of similiarity between target...
#              and reference segments 
#            - find the segment that has the highest value of both relative
#            distance and resamblance
# Input:  intSegments - segments of spectrum of interest
#         refSegments- segments of reference spectrum
#         intSp - spectrum of interest
#         refSp - reference spectrum 
# Output: intsegment$refInd - reference segment or []
# Author: L. Hedjazi, ICAN Paris 2013
# 
intSegLength<-length(intSegments$Peaks)
rC1<-vector()

for (index in 1:intSegLength)
 {
        testSeg<-list()
        testSeg$start<- intSegments$start[index]
        testSeg$PeakLeftBoundary<- intSegments$PeakLeftBoundary[[index]]
        testSeg$PeakRightBoundary<- intSegments$PeakRightBoundary[[index]]
        testSeg$Peaks<- intSegments$Peaks[[index]]
        testSeg$end<- intSegments$end[index]
        testSeg$centre<- intSegments$centre[index]

    peaksCompared<-comparePeaks(refSp, refSegments, intSp, testSeg, MAX_DIST_FACTOR, FALSE)

    if ((peaksCompared$rC>MIN_RC) && (peaksCompared$rC!=0))
    {
        intSegments$refIndex[index]<-peaksCompared$iSimilarPeakInd
        rC1<-cbind(rC1,peaksCompared$rC)
        refSegments$used[peaksCompared$iSimilarPeakInd]<-index
    }else{
        intSegments$refIndex[index]<-NA
    }
}

startPos<-numeric(0)
endPos<-numeric(0)

return(list(testSegs=intSegments,refSegs=refSegments))

}



comparePeaks<-function(dpReference,refPeaks, dpSpectr, curIPeak, MAX_DIST_FACTOR, tryDoubleReference)
{
# function iSimilarPeakInd <- comparePeaks(SPeakCurrent, refPeaks)
# compare current peak SPeakCurrent with all the peaks in refPeaks array and
# return index of similar peak
# dpReference - reference spectrum (whole!)
# dpSpectr - 'input' spectrum (whole!)
# curPeak - structure (.start, .end, .centre) of the current peak
# which we're aligning to our reference (one peak)
# refPeaks cell array of reference peak structures. We're looking for the
# best matching peak in these reference peak structures.
# MAX_DIST_FACTOR - 'zero' of distance
# Author: L. Hedjazi, ICAN Paris 2013

if (nargs()< 5)
{
    debug<-TRUE
    print('comparePeaks in debug mode!!!')
    MAX_DIST_FACTOR <- 50
    #load('../data/testPeaks')
    dpReference <- dpSpectr
}

iPeaksCount <- length(refPeaks$Peaks)
if (iPeaksCount<1)
{
    stop('Not enought peaks to compare')
}

#evluate comparison parameters

maxDistDeviation <- (curIPeak$end - curIPeak$start) * MAX_DIST_FACTOR
dMinParam   <- rep(0, iPeaksCount)
corCoeffs   <- rep(0, iPeaksCount)
distCoeffs  <- rep(0, iPeaksCount)

for (peakInd in 1:iPeaksCount)
{
refPeak<-list()
    if("used" %in% names(refPeaks)&& !(is.na(refPeaks$used[peakInd])))
    {
        #no matter to check this reference sicne it has been used already
        next
    }

    if((refPeaks$centre[peakInd] - curIPeak$centre) < -maxDistDeviation)
     {  next
     }
  
    if((refPeaks$centre[peakInd] - curIPeak$centre) > maxDistDeviation)
     {   break
     }

    if(tryDoubleReference && (peakInd < iPeaksCount))
     {
        refPeak$start   <- refPeaks$start[peakInd]
        refPeak$end     <- refPeaks$end[peakInd+1]
        refPeak$centre <- mean(refPeak$start, refPeak$end)
    }else{
        refPeak$start<- refPeaks$start[peakInd]
        refPeak$PeakLeftBoundary<- refPeaks$PeakLeftBoundary[[peakInd]]
        refPeak$PeakRightBoundary<- refPeaks$PeakRightBoundary[[peakInd]]
        refPeak$Peaks<- refPeaks$Peaks[[peakInd]]
        refPeak$end<- refPeaks$end[peakInd]
        refPeak$centre<- refPeaks$centre[peakInd]
    }
    #evaluate relative distance

    dDistCurr <- 1 - abs(refPeak$centre - curIPeak$centre) / maxDistDeviation
    distCoeffs[peakInd] <- dDistCurr
    
    ##-- little optimisation: if we got negative 'dDistCurr' value (our
    ##peaks are TOO far from each other, no matter to do FFT
    ##cross-correlation for them since we're interested only in positive rC
    ##values
    if(dDistCurr < 0)
     {
        dMinParam[peakInd] <- dDistCurr
        next
    }
    
    #evaluate cross-correlation
    refPeakShape <- dpReference[refPeak$start : refPeak$end]

    if (curIPeak$start<=0||curIPeak$end>length(dpSpectr))
     {
        dMinParam[peakInd]<-0
    }else{
        targetPeak <- dpSpectr[curIPeak$start : curIPeak$end]

        maxLen <- max(length(refPeakShape), length(targetPeak))
        fCorrCurr <- getCorellation(zeroPad(refPeakShape, maxLen), zeroPad(targetPeak, maxLen), maxDistDeviation)

        #get minimal parameter
        dMinParam[peakInd] <- min(dDistCurr,fCorrCurr)

        corCoeffs[peakInd]  <- fCorrCurr
    }
    }

#Get simillar peak index
rC <- max(dMinParam)
iSimilarPeakInd<-which.max(dMinParam)
distCoeff <- distCoeffs[iSimilarPeakInd]
corrCoeff <- corCoeffs[iSimilarPeakInd]


return(list(rC=rC, iSimilarPeakInd=iSimilarPeakInd, distCoeff=distCoeff, corrCoeff=corrCoeff))

}

getCorellation<-function(dpReferencePeak, dpInputPeak, maxShift)
{
if (nargs()<3)
 {
    stop('Invalid input parameters count')
 }
if (length(dpReferencePeak)!=length(dpInputPeak))
 {
    stop('dpReferencePeak and dpInputPeak sizes must agree')
 }

#evaluate lag by Wong
lag <- FFTcorr( dpInputPeak,dpReferencePeak, maxShift)
aligned <- shift(dpInputPeak,lag)

#get corellation coefficients
corellation <- cor(dpReferencePeak, aligned)

return(corellation)
}

########
zeroPad<-function(peak, maxLen)
{
if(!is.vector(peak))
 {
    stop('!isvector(peak)')
 }

if(!is.null(dim(peak)))
{
    stop('peak should be a single row vector!')
}

zerosCnt <- maxLen - length(peak)
if(zerosCnt > 0)
 {
    leftPaddCnt <- floor(zerosCnt /2 )
    rightPaddCnt <- zerosCnt - leftPaddCnt
    peak <- c(rep(0,leftPaddCnt),peak,rep(0,rightPaddCnt))
 }
return(peak)

}

####### FFT cross-correlation###########
FFTcorr<-function(testSegment, target,shift)
{
#padding
M<-length(testSegment)
diff <- 2^(ceiling(log2(M))) - M
testSegment<-testSegment-min(testSegment)
target<-target-min(target)

# append our ref & test segments with zeros at the end.
target<-c(target,rep(0,diff))
testSegment<-c(testSegment,rep(0,diff))
M<- M+diff
X<-fft(target)
Y<-fft(testSegment)
R<-X*Conj(Y)
R<-R/(M)
rev<-fft(R, inverse = TRUE)/length(R)
vals<-Re(rev)
maxpos <- 1
maxi <- -1
if (M<shift)
 {
    shift <- M
 }

for (i in 1:shift)
{
    if (vals[i]>maxi)
    {
        maxi <- vals[i]
        maxpos <- i
    }
    if (vals[length(vals)-i+1] > maxi)
     {
        maxi <- vals[length(vals)-i+1]
        maxpos <- length(vals)-i+1
    }
}

if (maxpos > (length(vals)/2))
 {
    lag <- maxpos-length(vals)-1
}else{
    lag <-maxpos-1
}

return(lag)
}

########shift segments#######################

shift<-function(seg, lag)
{
if ((lag == 0) || (lag >= length(seg)))
{
    shiftedSeg <- seg
    return(shiftedSeg)
}

if (lag > 0)
{
    ins <- rep(1,lag) * seg[1]
    shiftedSeg <-c(ins,seg[1:(length(seg) - lag)])
}else{ if (lag < 0)
    {
    lag <- abs(lag)
    ins <- rep(1,lag) * seg[length(seg)]
    shiftedSeg <-c(seg[(lag+1):length(seg)],ins)
    }
}

return(shiftedSeg)
}

corrcoef_aligned<-function(sp,ref,step)
{
# Calculation of correlation coefficient:
# Filter scales the spectral regions through specified step to unit variance 
# so as the high intensive and low intensive peaks would contribute 
# equally in the similarity


ilength<-length(sp)
if (ilength<20) {CC<-0}

if (step>=ilength)
 {
    CC<-cor(sp,ref)
    CC<-CC[1]
    return(CC)
}

bin_count<-ceiling(ilength/step)
bin_width<-ceiling(ilength/bin_count)
bins<-seq(1,ilength, bin_width)

if (bins[length(bins)]!=ilength)
{
    bins<-c(bins,ilength)
    bin_count<-bin_count+1
}

for (i in 1:(bin_count-1))

{
    istart<-bins[i]
    iend<-bins[i+1]-1
    sp[istart:iend]<-sp[istart:iend]-mean(sp[istart:iend])
    ref[istart:iend]<-ref[istart:iend]-mean(ref[istart:iend])

 if (istart!= iend)
  {
    if (var(sp[istart:iend])!=0)
    {
        sp[istart:iend]<-sp[istart:iend]/sd(sp[istart:iend])
    }
    if (var(ref[istart:iend])!=0)
     {
        ref[istart:iend]<-ref[istart:iend]/sd(ref[istart:iend])
    }
   }
}

CC<-cor(sp[1:(length(sp)-1)],ref[1:(length(ref)-1)])
CC<-CC[1]

return(CC)

}


#########find position to divide segment#########################
findMid<-function(testSeg,refSeg)
{
specLn<-length(testSeg)
M<-ceiling(length(testSeg)/2)
specM<-testSeg[(M-floor(M/4)):(M+floor(M/4))]
refM<-refSeg[(M-floor(M/4)):(M+floor(M/4))]

I<-which.min(specM*refM)
#[C,I]<-min(specM)
mid <- I[1]+M-floor(M/4)


#move to a point of a local minima
index<-1

while (TRUE)
{
    if (((mid-1) <= 1)||((mid+1) >= specLn))
     {
        mid<-NA
        break
    }

    if ((testSeg[mid] <= testSeg[mid+1])&&(testSeg[mid] <= testSeg[mid-1]))
       { 
        break}
   
    if (testSeg[mid] >= testSeg[mid+1])
     {
        mid<-mid+1
    }else{ if(testSeg[mid]>=testSeg[mid-1])
        {
         mid<-mid-1}
    }
}

return(mid)

}


######
localRecurAlign<-function(testSegment, refSegment,recursion,isSegment,lookahead)
{
# Input: recursion$minSegWidth
#                  minInbetweenWidth
#                  resamblance
#                  acceptance
#                  segShift
#                  inbetweenShift
#       isSegment == true takes segment parameters
# Author: L.Hedjazi, ICAN Paris 2013

if (!is.vector(testSegment)) {stop('!is.vector(testSegment)')}

if (!is.vector(refSegment)) {stop('!is.vector(refSegment)')}


if (length(refSegment) != length(testSegment))
 {
    stop('Reference and testSegment of unequal lengths')

}else {if (length(refSegment)== 1)
    {
    stop('Reference cannot be of length 1')
    }
}

recursion$minWidth<-recursion$minSegWidth

if (isSegment==TRUE)
{
    recursion$shift<-recursion$segShift
}else{
    recursion$shift<-recursion$inbetweenShift
}


alignedSegment<-recurAlign(testSegment,refSegment,recursion,lookahead)

return(alignedSegment)

}

#########Recursive segmentation################

recurAlign<-function(testSegment, refSegment, recursion, lookahead)
{

if (length(testSegment) < recursion$minWidth)
{
    alignedSegment <- testSegment
    return(alignedSegment)
}

if (var(testSegment)==0 | var(refSegment)==0)
{
    alignedSegment<-testSegment
    return(alignedSegment)
}


lag <- FFTcorr(testSegment,refSegment,recursion$shift)


## stop if the segment is perfectly aligned and there is no need to lookahead

alignedSegment <- testSegment


if (abs(lag) < length(testSegment))
{
    alignedTemp <- shift(testSegment,lag)

    if (var(alignedTemp)<0.000001) 
   {
       return(alignedSegment)
    }

    CorrCoef<-corrcoef_aligned(refSegment,alignedTemp,recursion$step)

    if (CorrCoef>=recursion$acceptance)
     {
  
        alignedSegment <-alignedTemp
    }else{
        if (var(testSegment)==0 || var(refSegment)==0)
         {

            alignedSegment<-testSegment
            return(alignedSegment)
        }
        CorrCoef<-corrcoef_aligned(refSegment,alignedSegment,recursion$step)
    }
}

CorrCoef<-corrcoef_aligned(refSegment,alignedSegment,recursion$step)


# Can be adjusted the recursion stops if the resemblance between the
# referebce and the segment of interest is e.g. 98%

if (CorrCoef>=recursion$resamblance) 
{return(alignedSegment)}

# If the middle point is not the local min then divide


mid <- findMid(alignedSegment,refSegment)



if (is.na(mid)){
return(alignedSegment)
}

firstSH<- alignedSegment[1 : mid]
firstRH<- refSegment[1 : mid]
secSH <- alignedSegment[(mid+1):(length(alignedSegment))]
secRH <- refSegment[(mid+1):(length(refSegment))]
alignedSeg1 <- recurAlign(firstSH,firstRH,recursion, lookahead)
alignedSeg2 <- recurAlign(secSH,secRH,recursion, lookahead)
alignedSegment <- c(alignedSeg1,alignedSeg2)


return(alignedSegment)

}

###### Spectrum Alignment#######
 alignSp<-function(refSp, refSegments, intSp,intSegments,recursion,MAX_DIST_FACTOR, MIN_RC)
{
# Input: refSp - reference spectrum
#        intSp - spectrum of interest
#        refSegments$used - reference
#        intSegments$refInd - matched reference segments
#        recursion - parameters for recursive alignment
# Output: alignedSpectrum
#         extendedSegments
#Author: L. Hedjazi, ICAN Paris 2013

if(length(refSp) != length(intSp))
{
    stop('length(refSp) != length(intSp)')
}

specLen <- length(refSp)

alignedSpectrum <- rep(NaN,specLen)
prevGeneralEnd <- 0

iSegmentInd <- 1
intSegLength<-length(intSegments$Peaks)
refSegLength<-length(refSegments$Peaks)
extensionCount<-0
extendedSegments<-NA
extensionInfo<-list()

###########################################
while (iSegmentInd <= intSegLength)
{
    ## Equi. to iSegment = intSegments[iSegmentInd]
    iSegment<-list()
     
    iSegment$start <- intSegments$start[iSegmentInd]
    iSegment$PeakLeftBoundary <- intSegments$PeakLeftBoundary[[iSegmentInd]]
    iSegment$PeakRightBoundary <- intSegments$PeakRightBoundary[[iSegmentInd]]
    iSegment$Peaks <- intSegments$Peaks[[iSegmentInd]]
    iSegment$end <- intSegments$end[iSegmentInd]
    iSegment$centre <- intSegments$centre[iSegmentInd]
    iSegment$refIndex <- intSegments$refIndex[iSegmentInd]

    if(is.na(iSegment$refIndex))
    {
      
        iSegmentInd <- iSegmentInd + 1
        next
    }
    ###### Segment of interest #######
    iLeftB<-iSegment$PeakLeftBoundary
    iRightB<-iSegment$PeakRightBoundary
    iPeaks<-iSegment$Peaks
    
    ###### Corresponding Reference segment ######
    referenceInd<-iSegment$refIndex



    #refSegment <- refSegments[referenceInd]


    refStart <- refSegments$start[referenceInd]
    refLeftB<-refSegments$PeakLeftBoundary[[referenceInd]]
    refRightB<-refSegments$PeakRightBoundary[[referenceInd]]
    rPeaks<-refSegments$Peaks[[referenceInd]]
    refCentre <- refSegments$centre[referenceInd]
    refUsed <- refSegments$used[referenceInd]
    iStart <- iSegment$start

    ######## Find joint starting position #########
    iSegmentIndex<-iSegmentInd
    refIndex<-referenceInd

    generalStart <- min(iStart, refStart)
    
    #search for the general start preventing overlapping with the previous segments
    while (TRUE)
     {
        #no segments before
        if ((iSegmentIndex<=1) && (refIndex<=1))
        { 
          break
        }

        # the segment of interest is first
        if (iSegmentIndex<=1)
         {
            if (generalStart<refSegments$end[refIndex-1])
             {
                generalStart<-min(generalStart,refSegments$start[refIndex-1])
                refLeftB<-c(refSegments$PeakLeftBoundary[[refIndex-1]], refLeftB)
                refRightB<-c(refSegments$PeakRightBoundary[[refIndex-1]],refRightB)
                rPeaks<-rbind(refSegments$Peaks[[refIndex-1]],rPeaks)
                extensionCount<-extensionCount+1
            }

            break
        }

        #the reference segment is first

        if (refIndex<=1)
        {
            if (generalStart<intSegments$end[iSegmentIndex-1])
             {
                generalStart<-min(generalStart,intSegments$start[iSegmentIndex-1])
                iLeftB<-c(intSegments$PeakLeftBoundary[[iSegmentIndex-1]],iLeftB)
                iRightB<-c(intSegments$PeakRightBoundary[[iSegmentIndex-1]],iRightB)
                iPeaks<-rbind(intSegments$Peaks[[iSegmentIndex-1]],iPeaks)
                extensionCount<-extensionCount+1
            }

            break
        }

        #both segments end before the general start
        if ((intSegments$end[iSegmentIndex-1]<=generalStart)&&(refSegments$end[refIndex-1]<=generalStart))
         {
           break
        }

        #both segments end after the general start (in fact impossible)
        if ((generalStart<intSegments$end[iSegmentIndex-1])&&(generalStart<=refSegments$end[refIndex-1]))
         {
            generalStart<-min(c(generalStart,refSegments$start[refIndex-1],intSegments$start[iSegmentIndex-1]))

            iLeftB<-c(intSegments$PeakLeftBoundary[[iSegmentIndex-1]], iLeftB)
            iRightB<-c(intSegments$PeakRightBoundary[[(iSegmentIndex-1)]],iRightB)
            refLeftB<-c(refSegments$PeakLeftBoundary[[refIndex-1]],refLeftB)
            refRightB<-c(refSegments$PeakRightBoundary[[refIndex-1]],refRightB)
            iPeaks<-rbind(intSegments$Peaks[[iSegmentIndex-1]],iPeaks)
            rPeaks<-rbind(refSegments$Peaks[[refIndex-1]],rPeaks)
            iSegmentIndex<-iSegmentIndex-1
            refIndex<-refIndex-1
            extensionCount<-extensionCount+1
            next
         }

        #the segment of interest ends after the general start
        if (generalStart<intSegments$end[iSegmentIndex-1])
         {
            generalStart<-min(generalStart,intSegments$start[iSegmentIndex-1])
            
            iLeftB<-c(intSegments$PeakLeftBoundary[[iSegmentIndex-1]],iLeftB)
            iRightB<-c(intSegments$PeakRightBoundary[[iSegmentIndex-1]],iRightB)
            iPeaks<-rbind(intSegments$Peaks[[iSegmentIndex-1]],iPeaks)
            iSegmentIndex<-iSegmentIndex-1
            extensionCount<-extensionCount+1
            next
        }

        #the reference segment ends after the general start
        if (generalStart<refSegments$end[refIndex-1])
         {
            generalStart<-min(generalStart,refSegments$start[refIndex-1])
            
            extensionCount<-extensionCount+1
            refLeftB<-c(refSegments$PeakLeftBoundary[[refIndex-1]],refLeftB)
            refRightB<-c(refSegments$PeakRightBoundary[[refIndex-1]],refRightB)
            rPeaks<-rbind(refSegments$Peaks[[(refIndex-1)]],rPeaks)
            refIndex<-refIndex-1
            next
        }
    }

    #search for 'generalEnd' preventing overlapping with the following segments
    iEnd <- iSegment$end
    refEnd <- refSegments$end[referenceInd]
    generalEnd <- max(iEnd,refEnd)


    while(TRUE)
     {
        # No segments ahead
        if ((iSegmentInd>=intSegLength)&&(referenceInd>=refSegLength))
         {
            break
         }

        # No segment ahead in spectrum of interest
        if (iSegmentInd>=intSegLength)
         {
            if (generalEnd>refSegments$start[referenceInd+1])
             {
                generalEnd<- max(generalEnd,refSegments$end[referenceInd+1])
                refLeftB<-c(refSegments$PeakLeftBoundary[[referenceInd+1]],refLeftB)
                refRightB<-c(refSegments$PeakRightBoundary[[referenceInd+1]],refRightB)
                rPeaks<-rbind(rPeaks,refSegments$Peaks[[(referenceInd+1)]])
                extensionCount<-extensionCount+1
                break
            }
            break
        }

        # No segment ahead in reference spectrum
        if (referenceInd>=refSegLength)
        {
            if (generalEnd>intSegments$start[iSegmentInd+1])
             {
                generalEnd<-max(generalEnd,intSegments$end[iSegmentInd+1])
                iLeftB<-c(iLeftB,intSegments$PeakLeftBoundary[[iSegmentInd+1]])
                iRightB<-c(iRightB,intSegments$PeakRightBoundary[[iSegmentInd+1]])
                iPeaks<-rbind(iPeaks,intSegments$Peaks[[iSegmentInd+1]])
                extensionCount<-extensionCount+1
                break
            }
            break
        }

        #Both subsequent segments start after the current general end

        if ((generalEnd <= intSegments$start[iSegmentInd+1]) && (generalEnd <= refSegments$start[referenceInd+1]))
         {
            break
        }

        #Both segments starts before the General End
        if ((generalEnd>intSegments$start[iSegmentInd+1])&&(generalEnd>refSegments$start[referenceInd+1]))
         {
            generalEnd<-max(c(generalEnd,intSegments$end[iSegmentInd+1],refSegments$end[iSegmentInd+1]))
            iLeftB<-c(iLeftB, intSegments$PeakLeftBoundary[[iSegmentInd+1]])
            iRightB<-c(iRightB, intSegments$PeakRightBoundary[[iSegmentInd+1]])
            refLeftB<-c(refLeftB,refSegments$PeakLeftBoundary[[referenceInd+1]])
            refRightB<-c(refRightB, refSegments$PeakRightBoundary[[referenceInd+1]])
            iPeaks<-rbind(iPeaks,intSegments$Peaks[[iSegmentInd+1]])
            rPeaks<-rbind(rPeaks,refSegments$Peaks[[referenceInd+1]])
            referenceInd<-referenceInd+1
            iSegmentInd<-iSegmentInd+1
            extensionCount<-extensionCount+1
            next
        }

        #If the next segment in intSp starts before the general end
        if ((generalEnd>intSegments$start[iSegmentInd+1])&&(is.na(intSegments$refIndex[iSegmentInd+1])))
        {
            generalEnd<-max(generalEnd,intSegments$end[iSegmentInd+1])
            iLeftB<-c(iLeftB,intSegments$PeakLeftBoundary[[iSegmentInd+1]])
            iRightB<-c(iRightB,intSegments$PeakRightBoundary[[iSegmentInd+1]])
            iPeaks<-rbind(iPeaks,intSegments$Peaks[[iSegmentInd+1]])
            iSegmentInd<-iSegmentInd+1
            extensionCount<-extensionCount+1
            next
        }else{  if (generalEnd>intSegments$start[iSegmentInd+1])
           {
            refInd<-referenceInd+1
            referenceInd<-intSegments$refIndex[iSegmentInd+1]
            generalEnd<-max(c(generalEnd,intSegments$end[iSegmentInd+1],refSegments$end[referenceInd]))
            iLeftB<-c(iLeftB,intSegments$PeakLeftBoundary[[iSegmentInd+1]])
            iRightB<-c(iRightB,intSegments$PeakRightBoundary[[iSegmentInd+1]])
            iPeaks<-rbind(iPeaks,intSegments$Peaks[[iSegmentInd+1]])

            for (i in (refInd:referenceInd))
             {
                refLeftB<-c(refLeftB,refSegments$PeakLeftBoundary[[i]])
                refRightB<-c(refRightB,refSegments$PeakRightBoundary[[i]])
                rPeaks<-rbind(rPeaks,refSegments$Peaks[[i]])

            }
            iSegmentInd<-iSegmentInd+1
            extensionCount<-extensionCount+1
            next
           }
          }

        #If the next segment in refSp starts before the general end

        if ((generalEnd>refSegments$start[referenceInd+1])&&(is.na(refSegments$used[referenceInd+1])))
         {  
            generalEnd<-max(generalEnd,refSegments$end[referenceInd+1])
            refLeftB<-c(refLeftB,refSegments$PeakLeftBoundary[[referenceInd+1]])
            refRightB<-c(refRightB,refSegments$PeakRightBoundary[[referenceInd+1]])
            rPeaks<-rbind(rPeaks,refSegments$Peaks[[referenceInd+1]])
            referenceInd<-referenceInd+1
            extensionCount<-extensionCount+1
            next
        }else{ if (generalEnd>refSegments$start[referenceInd+1])
           {
            iSegIndex<-iSegmentInd+1
            iSegmentInd<-refSegments$used[referenceInd+1]
            generalEnd<-max(c(generalEnd,intSegments$end[iSegmentInd],refSegments$end[referenceInd+1]))
            for (i in (iSegIndex:iSegmentInd))
             {
                iLeftB<-c(iLeftB,intSegments$PeakLeftBoundary[[i]])
                iRightB<-c(iRightB,intSegments$PeakRightBoundary[[i]])
                iPeaks<-rbind(iPeaks,intSegments$Peaks[[i]])

             }
            refLeftB<-c(refLeftB,refSegments$PeakLeftBoundary[[referenceInd+1]])
            refRightB<-c(refRightB,refSegments$PeakRightBoundary[[referenceInd+1]])
            rPeaks<-rbind(rPeaks,refSegments$Peaks[[referenceInd+1]])
            referenceInd<-referenceInd+1
            extensionCount<-extensionCount+1
            next
            }
          }
      }

    refSegment    <- refSp[generalStart : generalEnd]
    testSegment   <- intSp[generalStart : generalEnd]

    Bnd<-list()
    Bnd$refLeftB<-refLeftB-generalStart+1
    Bnd$refRightB<-refRightB-generalStart+1
    Bnd$iLeftB<-iLeftB-generalStart+1
    Bnd$iRightB<-iRightB-generalStart+1

    alignedSegment <- localRecurAlign(testSegment, refSegment, recursion,TRUE,TRUE)  


    if (any(is.nan(alignedSegment)))
     {
        readline("any(is.nan(alignedSegment))")
    }
    alignedSpectrum[generalStart : generalEnd] <- alignedSegment

    #############################################################

    ##-- align 'grass':

    grassStart  <- prevGeneralEnd + 1
    grassEnd    <- generalStart - 1

    if( grassEnd > grassStart )
    {
        refSegment <- refSp[ grassStart : grassEnd ]
        testSegment  <- intSp[ grassStart : grassEnd ]

        # do not want to visualize grass

        alignedSegment <- localRecurAlign(testSegment, refSegment, recursion,FALSE,TRUE)
        alignedSpectrum[grassStart : grassEnd] <- alignedSegment
    }
    prevGeneralEnd <- generalEnd

    # don't forget to increase the counter!!!
    iSegmentInd <- iSegmentInd + 1
}

if(extensionCount > 0)
 {
    extensionInfo$extensionCount <- extensionCount
 } 


##########################################################
if(!(is.na(extendedSegments)))
 {
    maxExtensionCount <- -1
    if (extendedSegment == extendedSegments)
     {
        if(extendedSegment$extensionCount > maxExtensionCount)
         {
            maxExtensionCount <- extendedSegment$extensionCount
            maxExtInd <- extendedSegment$extensionSegmentInd
         }
     }
    maxExtensionInfo$extensionSegmentInd <- maxExtInd
    maxExtensionInfo$extensionCount <- maxExtensionCount
    extendedSegments <- c(extendedSegments, maxExtensionInfo)
  }

grassStart<- prevGeneralEnd + 1
grassEnd <- specLen

if( grassEnd > grassStart )
{
    refSegment<- refSp[grassStart : grassEnd]
    testSegment<- intSp[grassStart : grassEnd]

    alignedSegment <- localRecurAlign(testSegment, refSegment, recursion, FALSE,TRUE)
    alignedSpectrum[grassStart : grassEnd]<- alignedSegment
}
#return(list(alignedSpectrum=alignedSpectrum, extendedSegments=extendedSegments))
 return(alignedSpectrum)
}

drop_outliers=function(data,ppm,n){
  #par(mfrow=c(2,1))
  nppm=length(ppm)
  for(p in 1:nppm){
  #  hist(lip.SRVlog10_1[p,],br=50,main=paste(p,": Before"))
    for (i in 1:n){
      data[p,]<-rm.outlier(data[p,] ,fill=T,median=T)
    }
   # hist(lip.SRVlog10_1[p,],br=50,main="After")
  }
 return(data)
}
format_mQTL=function(datafile, genofile, physdat,outdat, outgeno){
  # Routine to reformat the data into the required csvs used by the R/QTL package
  print(paste("Start formatting the datafile",datafile,"and the genotype file",genofile,"into the cvsv files:",outdat,outgeno))

  ur.dat<-read.table(datafile,as.is=T,header=T,sep="\t")
  ur.geno<-read.table(genofile,as.is=T,header=T,na.strings="ND",sep="\t")
  ur.mrks=colnames(ur.geno[,-1])
  ur.genmouse=paste("X",ur.geno[-c(1,2,3),1],sep="")

  #ur.gen<-matrix(9,nrow=dim(ur.geno.CHD)[[1]]+dim(ur.geno.HFD)[[1]]-5,ncol=(dim(ur.geno.CHD)[[2]]+dim(ur.geno.HFD)[[2]]-1)
  ur.gen<-matrix("NA",nrow=dim(ur.geno)[[1]]-3,ncol=dim(ur.geno)[[2]]-1)
  ur.gen[ur.geno[-c(1,2,3),-1]=="a"]="A"
  ur.gen[ur.geno[-c(1,2,3),-1]=="h"]="H"
  ur.gen[ur.geno[-c(1,2,3),-1]=="b"]="B"
  dim(ur.gen)<-c(dim(ur.geno)[[1]]-3,dim(ur.geno)[[2]]-1)
  row.names(ur.gen)<-ur.genmouse
  colnames(ur.gen)<-ur.mrks
  # Define the set of samples to use with genotype, phenotype and data, in the genotype order
  ur.set<-intersect(ur.genmouse,colnames(ur.dat))
  ur.dat.names<-grep("X",unlist(strsplit(ur.set,"_")),value=T)
  ur.gen2<-ur.gen[ur.set,]
  #write.table(ur.gen2,file="ur_geno.ehk.dat",sep="\t",row.names=F,col.names=F)
  
  ur.data<-ur.dat[,ur.set]
  ur.nm<-dim(ur.data)[2]
  
  ## PREPARE THE DATA
  ur.chr<-ur.geno[1,-1]
  ur.loc<-ur.geno[2,-1]
  
  # Export the cleaned data in the proper csvs format
  ur.dat2<-t(ur.data)
  ur.p<-as.numeric(colnames(ur.dat2))


   ur.extr<-read.table(physdat,as.is=T,header=T,sep="\t")
   sex=ur.extr[1,ur.set]
   pgm=ur.extr[2,ur.set]

   ur.dat3<-cbind(ur.dat2, t(sex), t(pgm), rownames(ur.gen2))



  colnames(ur.dat3)<-c(paste("ppm",ur.p,sep="_"),"sex","pgm","id")
  write.csv(ur.dat3,outdat,quote=F,row.names=F)
  
  ur.gen3<-rbind(ur.chr, ur.loc, ur.gen2)
  ur.gen3<-cbind(c("","",rownames(ur.gen3[-c(1,2),])), ur.gen3)
  colnames(ur.gen3)[1]<-"id"
  write.csv(ur.gen3, outgeno, quote=F, row.names=F)
}

pre_mQTL=function(infile, outfile, met,corrT=0.9 ){
  # Routine to preprocess the NMR file into SRV 
 
  dat<-read.csv(infile)
  ind=row.names(dat)
  n_ind=length(ind)

  np<-dim(dat)[[2]]-3
  dat2<-as.matrix(dat[,1:np])
  data<-matrix(nrow=n_ind, ncol=np)
  p<-as.numeric(sub("e\\.","e-",sub("ppm_","",sub("ppm_\\.","-", colnames(dat2)))))
  dat2[which(is.na(dat2))]=0
  if(met!="running"){
    #Perform the SRV analysis to reduce the number of dimension
    meth<-get(met)
  
    SRV<-SRV(dat2, 10, corrT,clustf=meth)
    data<-SRV[[3]]
    #data<-log10(SRV[[3]]+abs(floor(min(SRV[[3]]))))
    ppm<-(p[SRV[[1]] ]+p[SRV[[2]] ])/2
    nppm=length(ppm)
    
    write.table(t(rbind(p[SRV[[1]] ] ,p[SRV[[2]] ],ppm)),paste(met,"SRV.ppm",sep="_"),col.names=c("Min","Max","Mean"),sep="\t",row.names=F,quote=F)
  
  }else{ # Use the original Running function rather than SRV
    mymid=25
    mybase=5
    for(r in 1:n_ind){
      data[r,] = running(dat2[r,],width=0.00025,mid=mymid,base=mybase)
    }
    ppm=p
    nppm=length(ppm)
    # Drop the border effect
    start=mymid+mybase+1
    end=nppm-mymid-mybase-1
    
    data[,1:start]=0 # Set borders to 0 
    data[,end:nppm]=0
  }

  SRVlog10_1 = log10(data+abs(floor(min(data)))+1)
  dim(SRVlog10_1)= c(n_ind,nppm)
  
  #Drop outliers
  #drop_outliers(SRVlog10_1,nppm,2)
  
  dat3<-cbind(SRVlog10_1,dat[1:3+np])
  colnames(dat3)<-c(paste("ppm",ppm,sep="_"),"sex","pgm","id")
  if(!exists("outfile")){outfile<-paste("ur",met,"dat",sep=".")}
  write.csv(dat3,outfile,quote=F,row.names=F)
}

process_mQTL=function(datfile, genfile, nperm=0){
  # Script to process the tissue extract of the individuals
  # Part of the main script
  # Input:  genotype and phenotype
  # Output: 2D LOD score table
  
  ## PROCESS THE DATA
   
  st=5
  err=0.001
  cross<-read.cross("csvs", genfile=genfile,phefile=datfile)
  np<-nphe(cross)-3
  ppm<-as.numeric(sub("e\\.","e-",sub("ppm_","",sub("ppm_\\.","-", names(cross$pheno)[1:np]))))
  k=1:np
  cross=calc.genoprob(cross, step=st, error.prob=err)
  cross=sim.geno(cross, step=st, error.prob=err,n.draws=64)
  TotM<-0
  for(i in 1:nchr(cross)){ TotM<-TotM+dim(cross$geno[[i]]$prob)[[2]] }
  print(paste("TotM is ",TotM,"and np is",np))
  res = matrix(nrow=TotM,ncol=np)
  ehk = matrix(nrow=TotM,ncol=np)
  permo = list()
  for (i in 1:np){
    best<-scanone(cross, pheno.col=i, method ="ehk")
        if (nperm>0){ 
       permo[[i]]<-scanone(cross,pheno.col=i,method ="ehk",n.perm=nperm,verbose=F) 
           }
    ehk[,i]<-best[,"lod"]
    if(i%%min(round(np/10),100)==0){
      print(c(i,max(ehk,na.rm=T)))
      summary(best)
    }
  }
  res[,k]=ehk
  top<-ceiling(which.max(res)/TotM)
  best[,"lod"]<-res[,top]
  maxi<-rownames(best)[which.max(res)-(top-1)*TotM]
  maxiL<-round(max(res),2)
  print(paste("SRV ehk in this area was",maxiL,"for",ppm[top],"ppm at",maxi))
  if(nperm>0){
    summary(best,perms=permo[[top]],pvalues=T)
  }else{
    summary(best,maxiL-1)
  }

  return(list(res=res,ppm=ppm,best=best,top=top,maxi=maxi,maxiL=maxiL,permo=permo))
}

post_mQTL=function(results, probs=c(0.95,0.99,0.999,0.9999)){
# Function to plot the results of a given run
# Input:  2D results of LOD scores
# Output: Graphs and summaries

  for(met in names(results)){
 
    print(paste("Post processing",met))
    for(i in 1:length(results[[met]])){assign(names(results[[met]])[i],results[[met]][[i]])}
 
    ## PLOT THE RESULTS
    
 #   x11(width=11,height=11,pointsize=10)
    split.screen(c(2,2))
    split.screen(c(2,1),screen=4)
    split.screen(c(2,1),screen=1)
    
    title<-paste("for All using",met,":\n",maxiL,"for",ppm[top],"ppm at",maxi)
    screen(7)
    qt<-quantile(res,probs=probs)
    h<-hist(res,br=50,freq=F,main="LOD Distribution",xlab="LOD")
    abline(v=qt,col=c(1,4,3,2),lt=3)
    text(qt,9:6*0.1*max(h$intensities),paste(format(1-probs),"%",sep=""),col=c(1,4,3,2),pos=3)
    #text(qt,0.9*max(h$intensities),c("1%","0.1%","0.01%","0.001%"),col=c(1,4,3,2),pos=4)
    
    screen(8)
    plot(ppm,res[which(res==max(res),arr.ind=T)[1] ,],lty=1,main=paste("LOD for top locus",maxi),ty="b",col="red",ylab="LOD",xlab="Resonance, ppm",pch=3)
    rug(ppm,ticksize=-0.03)
    
    screen(6)
    plot(best,lty=1,main=paste("LOD for top Shift",ppm[top],"ppm"))
    
    screen(5)
    topc<-summary(best)[which.max(summary(best)[,3]),1]
    plot(best,lty=1,main=paste("LOD for top Shift",ppm[top],"ppm on chr",topc),chr=topc,col="blue")
    
    screen(3)
    pplot(res,"Full 2D Profile", ppm, best, quantile(res,probs=probs))
    
    screen(2)
    ppersp(res, ppm, paste("Top LOD",title))
    
    close.screen(all.screens=T)
  }
}

summarize=function(resu, met, Th=5){
# Function to summarize a specific result called by summary_mQTL
# Input:  2D results of LOD scores
# Output: Summary

    for(i in 1:length(resu)){assign(names(resu)[i],resu[[i]])}
 
    overT<-unique(which(res>=Th,arr.ind=T)[,2])
    if(length(overT)==0){
       overT<-unique(which(res>=max(res)*4/5,arr.ind=T)[,2])
    }
    summar<-data.frame(check.names=F)
    min_max<-read.table(paste(met,"SRV.ppm",sep="_"),header=T)
    for(i in overT){
      sc<-best
      sc[,"lod"]<-res[,i]
      if(length(permo)>0){
        summar<-rbind(summar,cbind(min_max[i,1],min_max[i,2],ppm[i],summary(sc, perms=permo[[i]], alpha=0.1, pvalues=T)))
      }else{
        summar<-rbind(summar,(cbind(min_max[i,1],min_max[i,2],ppm[i],summary(sc,Th))))
      }
    }
    colnames(summar)[1:3]<-c("Min","Max","ppm")
    return(summar)
}

summary_mQTL=function(results, Th=5){
# Function to summarize the results of a all the runs and their differences
# Input:  2D results of LOD scores
# Output: Summaries

  for(met in names(results)){
 
    print(paste("Summarize single runs",met))
    summa<-summarize(results[[met]], met, Th)
    print(summa)
    write.table(summa, file=paste("signif",met,"sum",sep="."),quote=F,sep="\t")
  }
#  print("Prepare the pairs")

}

run_QTL=function(genofile,data,group,pgm,ppm,wd,nm,st=5,err=0.01,nperms=0){

  np=length(ppm)
  res = matrix(nrow=nm,ncol=np)
  k=1:np
  cross=read.cross(genfile=genofile,phefile=cbind(t(data[k,group]),pgm), format="gary",pnamesfile=c(paste(ppm[k],rep("ppm",np)),"pgm"),dir=wd)
  cross=calc.genoprob(cross, step=st, error.prob=err)
  ehk = matrix(nrow=nm,ncol=np)
  for (i in 1:np){
    #lip.male.ehk[,i]<-scanone(lip.male.cross,pheno.col=i,method ="ehk",chr=c('-X'))[,"lod"]
    best<-scanone(cross,pheno.col=i,method ="ehk",addcovar=pgm,maxit=100000)
    ehk[,i]<-best[,"lod"]
  }
  res[,k]=ehk
  top<-ceiling(which.max(res)/nm)
  best[,"lod"]<-ehk[,top]
  maxi<-rownames(best)[which.max(res)-(top-1)*nm]
  maxiL<-round(max(res),2)
  print(paste("SRV Sex corrected ehk in this area was",maxiL,"for",ppm[top],"ppm at",maxi))
  if (nperms>0){
    permo<-scanone(cross,pheno.col=top,method ="ehk",addcovar=pgm,n.perm=nperms)
  }else{
    permo=NULL
  }
  summary(best,maxiL-1)
  return( list(res=res,best=best,top=top,maxi=maxi,maxiL=maxiL,permo=permo))
}

panel.cor <- function(x, y, digits=2, prefix="", cex.cor)
{
    usr <- par("usr"); on.exit(par(usr))
    par(usr = c(0, 1, 0, 1))
    r <- abs(cor(x, y))
    txt <- format(c(r, 0.123456789), digits=digits)[1]
    txt <- paste(prefix, txt, sep="")
    if(missing(cex.cor)) cex.cor <- 0.8/strwidth(txt)
    text(0.5, 0.5, txt, cex = cex.cor * r)
}

panel.hist <- function(x, ...)
{
    usr <- par("usr"); on.exit(par(usr))
    par(usr = c(usr[1:2], 0, 1.5) )
    h <- hist(x, plot = FALSE)
    breaks <- h$breaks; nB <- length(breaks)
    y <- h$counts; y <- y/max(y)
    rect(breaks[-nB], 0, breaks[-1], y, col="cyan", ...)
}

Plot_summary=function(res,title){
  
  pairs(res ,upper.panel= panel.cor, lower.panel=panel.smooth,cex = 1.5,diag.panel=panel.hist, main=title)
}

Centered=function(A){
  #Centered Centre les variables d'une matrice.
  #   Class support for input A:
  #      float: double, single
  
  no<-dim(A)[[1]] 
  s <-  sum(A);
  m <- s/no;
  
  Ca<-A-m;   

  return(Ca)
}

UVscaled=function(A){
  #UVscaled Centre et reduit les variables d'une matrice.
  #   Class support for input A:
  #      float: double, single

  no<-dim(A)[[1]] 
  nv<-dim(A)[[2]] 
  s <-  sum(A);
  m <- s/no;
  Ua<-matrix(0,nrow=no,ncol=nv);
  m<-rep(0,nv)
  ET<-rep(0,nv)
  for (j in 1:nv){
    ET[j]<-sd(A[,j]);
  }
 
  for (j in 1:nv){
    M<-m[j]*rep(1,no);
    Ua[,j]<-(A[,j]-M)/ET[j];
  }

  return(Ua)
}

SRV=function(X, minsize, correl, clustf=median){
  # Statistical Recoupling of Variables.
  # input: data matrix X, singlet size, bucketting resolution, correlation threshold
  # output: supercluster borders: indicesdebf and indicesfinf; 
  # superclusters.

  dX<-dim(X);
  X[is.nan(X)]<-0
  dim(X)<-dX
  Xct<-Centered(X)
  Xuv<-UVscaled(X)
  nr<-dim(X)[[1]]
  nc<-dim(X)[[2]]
  B<-rep(0,nc-1)

  print(paste("Processing a",nr,"x",nc,"matrix with minsize of",minsize," and a correlation Threshold of",correl,"using the",met))

  # Calculation of the covariance/correlation profile between consecutive variables
  for (j in 1:(nc-1)){
    B[j]<-var(Xct[,j],Xct[,j+1])/abs(cor(Xuv[,j],Xuv[,j+1]));
  }

  #Identification of the residual water area
  wat<-which(apply(X,2,sum)==0)
  if (length(wat)>0){
    lsup<-min(wat)-1
    linf<-max(wat)

    B[lsup:linf]<-0;
  }

  #Scan of the profile for the identification of local minima
  #starting and ending variables of each cluster are stored in indicesdeb and
  #indicesfin according to the covariance/correlation profile.
  #the first variable belongs to the first cluster and the last variabel to
  #the last cluster.  
  indicesdeb<-rep(0,nc);
  indicesfin<-rep(0,nc);
  ndeb<-1;
  nfin<-0;

  indicesdeb[ndeb]<-1;
  for (j in 2:(nc-2)){
    if ((B[j]<B[j-1]) && (B[j]<B[j+1])){
      nfin<-nfin+1;
      indicesfin[nfin]<-j;
      ndeb<-ndeb+1;
      indicesdeb[ndeb]<-j+1;
    }
  }
#  length(indicesdeb)<-ndeb
#  length(indicesfin)<-nfin

  if(indicesfin[nfin]!=nc){
  	nfin=nfin+1
    indicesfin[nfin]<-nc;
  }
  length(indicesdeb)<-ndeb
  length(indicesfin)<-nfin

  #Measurement of the size of each cluster.
  #Comparison of the size with a resolution criterium:
  #floor(singletsize/resolution), where singletsize is the size of a resolved singlet in a 700MHz spectrum.
  #Clusters that do not contain at least a number of variables equals to "limit" are discarded
  indices<-indicesfin-indicesdeb+1; 
  si<-length(indices);
  destruction<-rep(0,si);
  ndestruction<-0;
  limit<-minsize
  for (j in 1:si){
    if(indices[j]<limit){
      ndestruction<-ndestruction+1;
      destruction[ndestruction]<-j;
    }
  }
  length(destruction)<-ndestruction
  if (ndestruction>0){
    indicesdeb=indicesdeb[-destruction];
    indicesfin=indicesfin[-destruction];
  }


  #Measurement of the mean value of the signal in each cluster stored in Xcluster.

  s1<-length(indicesdeb)
  print(paste("length of events:",s1))
  Xcluster<-matrix(0,nrow=nr,ncol=s1);
  jbstore<-matrix(0,nrow=nr*s1,ncol=4)
  colnames(jbstore)<-c("Max","Sum","Mean","Median")

  for (j in 1:s1){
    Xcluster[,j]<-apply(X[,indicesdeb[j]:indicesfin[j]],1,clustf)
    ind<-(j-1)*nr+1
    jbstore[ind:(ind+nr-1),1]<-as.vector(apply(X[,indicesdeb[j]:indicesfin[j]],1,max))
    jbstore[ind:(ind+nr-1),2]<-as.vector(apply(X[,indicesdeb[j]:indicesfin[j]],1,sum))
    jbstore[ind:(ind+nr-1),3]<-as.vector(apply(X[,indicesdeb[j]:indicesfin[j]],1,mean))
    jbstore[ind:(ind+nr-1),4]<-as.vector(apply(X[,indicesdeb[j]:indicesfin[j]],1,median))
  }

  print(c("Xcluster",dim(Xcluster),s1))

  #Correlation of neighboring clusters stored in Clustercorrelation.
  #Identification of highly correlated clusters.
  #We limit the agregation of clusters to 3 clusters according to a
  #sufficient level of correlation (0.9) that allow the identification of
  #chemically relevant superclusters (singlet, doublet, triplet, mulatiplet) on a 1D spectrum.
  Clustercorrelation<-rep(0,s1-1); 
  for (j in 1:(s1-1)){
    Clustercorrelation[j]<-abs(cor(c(Xcluster[,j],0.123456789),c(Xcluster[,j+1],0.1234566789)));
      # Add a hopefully unlikely constant to both vector to avoid a constant vector which result in NA correlation
  }
  
  n2<-0;
  balises<-rep(0,s1);
  for (j in 1:(s1-1)){
    if(Clustercorrelation[j]>correl){
      n2<-n2+1;
      balises[n2]<-j;
    }
  }
  length(balises)<-n2

  print(c("balise length:",n2))
  
  s2<-length(balises)
  choixbalises<-rep(0,s2);
  if(s2>1){ # If there are more than 1 supercluster
    for (k in 2:(s2-1)){
      if((balises[k]+1)==(balises[k+1])){
        choixbalises[k]<-1;
      }else{
#        if((balises[k]+1)!=(balises[k+1])){
          choixbalises[k]<-0;
#        }
      }
    }
    if (balises[1]+1== balises[2]){
      choixbalises[1]<-1;
    }else{ 
  	choixbalises[1]<-0;
    }

    if (balises[s2]-1== balises[s2-1]){
      choixbalises[s2]<-1;
    }else{ 
    	choixbalises[s2]<-0;
    }

    for (j in 2:(s2-1)){
      if ((choixbalises[j-1]==choixbalises[j]) &&(choixbalises[j+1]==choixbalises[j])){
        choixbalises[j]<-0;
      }
    }

    debcluster2<-rep(0,s2);
    fincluster2<-rep(0,s2);
    ndeb2<-0;
    nfin2<-0;
    for (k in 2:(s2-1)){
      if (choixbalises[k]==0 && choixbalises[k-1]!=1){
        ndeb2<-ndeb2+1;
        nfin2<-nfin2+1;
        debcluster2[ndeb2]<-balises[k];
        fincluster2[nfin2]<-balises[k];
      }else{
        if (choixbalises[k]>choixbalises[k-1]){
          ndeb2<-ndeb2+1;
          debcluster2[ndeb2]<-balises[k];
        }else{  
          if (choixbalises[k]==0 && choixbalises[k-1]==1){
            nfin2<-nfin2+1;
            fincluster2[nfin2]<-balises[k];
          }
        }
      } 
    }
    length(debcluster2)<-ndeb2
    length(fincluster2)<-nfin2

    if(choixbalises[s2-1]==1){
      nfin2<-nfin2+1;
      fincluster2[nfin2]<-balises[s2];
    }else{
      nfin2<-nfin2+1;
      ndeb2<-ndeb2+1;
      fincluster2[nfin2]<-balises[s2];
      debcluster2[ndeb2]<-balises[s2];
    }

    if(choixbalises[1]==1){
      debcluster2<-c(balises[1],debcluster2);
    }else{
      debcluster2<-c(balises[1],debcluster2);
      fincluster2<-c(balises[1],fincluster2);
    }

    #Measurement of the mean value of the clusters involved in a supercluster,
    #stored in Xcluster2.
    s3<-length(debcluster2);
    Xcluster2<-matrix(0,nrow=nr,ncol=s3);
    jbstore2<-matrix(0,nrow=nr*s3,ncol=4)
    colnames(jbstore2)<-c("Max","Sum","Mean","Median")
    for (j in 1:s3){
      Xcluster2[,j]<-apply(Xcluster[,debcluster2[j]:(fincluster2[j]+1)],1,clustf)
      ind<-(j-1)*nr+1
      jbstore2[ind:(ind+nr-1),1]<-as.vector(apply(Xcluster[,debcluster2[j]:(fincluster2[j]+1)],1,max))
      jbstore2[ind:(ind+nr-1),2]<-as.vector(apply(Xcluster[,debcluster2[j]:(fincluster2[j]+1)],1,sum))
      jbstore2[ind:(ind+nr-1),3]<-as.vector(apply(Xcluster[,debcluster2[j]:(fincluster2[j]+1)],1,mean))
      jbstore2[ind:(ind+nr-1),4]<-as.vector(apply(Xcluster[,debcluster2[j]:(fincluster2[j]+1)],1,median))
    }
  }else{ #If there is only 1 correlation
    if(s2==1){
      debcluster2<-balises[1]
      fincluster2<-balises[1]
      s3<-1
      Xcluster2<-matrix(0,nrow=nr,ncol=s3);
      jbstore2<-matrix(0,nrow=nr*s3,ncol=4)
      colnames(jbstore2)<-c("Max","Sum","Mean","Median")
      Xcluster2[,1]<-Xcluster[,debcluster2[1]]
      ind<-1
      print("jbstore2 single")
      jbstore2[,1]<-Xcluster[,debcluster2[1]]
      jbstore2[,2]<-Xcluster[,debcluster2[1]]
      jbstore2[,3]<-Xcluster[,debcluster2[1]]
      jbstore2[,4]<-Xcluster[,debcluster2[1]]
    }
  }

  #Clusters and superclusters are finally stored in Xclusterf, after 
  #adjustment of the cluster borders in case of overlapping after agregation.
  #Xclusterf<-matrix(0,nrow=nr,ncol=s3);
  Xclusterf<-matrix(0,nrow=nr,ncol=0);
  indicesdebf<-NULL;
  indicesfinf<-NULL;
  if (s2>1){ # If there are more than one superclusters Xcluster2
    if (debcluster2[1]>1){
      Xclusterf<-Xcluster[,1:(debcluster2[1]-1)];
      indicesdebf<-indicesdeb[1:(debcluster2[1]-1)];
      indicesfinf<-indicesfin[1:(debcluster2[1]-1)];
    }
    for (j in 1:(s3-1)){
      if (fincluster2[j]+2<=debcluster2[j+1]-1){
        Xclusterf<-cbind(Xclusterf,Xcluster2[,j],Xcluster[,(fincluster2[j]+2):(debcluster2[j+1]-1)]);
        indicesdebf<-c(indicesdebf,indicesdeb[debcluster2[j]],indicesdeb[(fincluster2[j]+2):(debcluster2[j+1]-1)]);
        indicesfinf<-c(indicesfinf,indicesfin[fincluster2[j]+1],indicesfin[(fincluster2[j]+2):(debcluster2[j+1]-1)]);
      }else{
        Xclusterf<-cbind(Xclusterf,Xcluster2[,j]);
        indicesdebf<-c(indicesdebf,indicesdeb[debcluster2[j]]);
        indicesfinf<-c(indicesfinf,indicesfin[fincluster2[j]+1]);
      }
    }
  
    if (fincluster2[s3]+2<=s1){
      Xclusterf<-cbind(Xclusterf,Xcluster2[,s3], Xcluster[,(fincluster2[s3]+2):s1]);
      indicesdebf<-c(indicesdebf,indicesdeb[debcluster2[s3]],indicesdeb[(fincluster2[s3]+2):s1]);
      indicesfinf<-c(indicesfinf,indicesfin[fincluster2[s3]+1],indicesfin[(fincluster2[s3]+2):s1]);
    }else{
      Xclusterf<-cbind(Xclusterf,Xcluster2[,s3])
      indicesdebf<-c(indicesdebf,indicesdeb[debcluster2[s3]]);
      indicesfinf<-c(indicesfinf,indicesfin[fincluster2[s3]+1]);
    }
  }else{
     if (s2==1){ # If there one supercluster Xcluster2
  #      print(paste("s1,s2,s3",s1,s2,s3))
  #      print("debcluster2, fincluster2")
  #      print(c(debcluster2, fincluster2))
  #      print("indicesdeb,indicesfin")
  #      print(c(indicesdeb,indicesfin))
        if (debcluster2[1]>1){
          Xclusterf<-Xcluster[,1:(debcluster2[1]-1)];
          indicesdebf<-indicesdeb[1:(debcluster2[1]-1)];
          indicesfinf<-indicesfin[1:(debcluster2[1]-1)];
        }
        Xclusterf<-cbind(Xclusterf,Xcluster2[,1]);
        indicesdebf<-c(indicesdebf,indicesdeb[debcluster2[1]]);
        indicesfinf<-c(indicesfinf,indicesfin[fincluster2[1]+1]);
        if (fincluster2[1]+2<=s1){
          Xclusterf<-cbind(Xclusterf,Xcluster[,(fincluster2[1]+2):s1]);
          indicesdebf<-c(indicesdebf,indicesdeb[(fincluster2[1]+2):s1]);
          indicesfinf<-c(indicesfinf,indicesfin[(fincluster2[1]+2):s1]);
        }
     }else{ # If there are no supercluser Xcluster2
       Xclusterf=Xcluster
       indicesdebf<-indicesdeb;
       indicesfinf<-indicesinf;
    }
  }

  nrt<-dim(Xclusterf)[[1]];
  nct<-dim(Xclusterf)[[2]];

  for (j in 1:(nct-1)){
    if (indicesdebf[j+1]<indicesfinf[j]){
      indicesfinf[j]<-indicesdebf[j+1]-1;
    }
  }

  print(paste("Found",length(indicesdebf)," clusters"))
  return(list(indicesdebf,indicesfinf, Xclusterf,jbstore,jbstore2,indicesdeb,indicesfin))
 # return(list(indicesdebf,indicesfinf,Xcluster,Xcluster2, Xclusterf,jbstore,jbstore2,indicesdeb,indicesfin))

}

SRV_Corr=function(X, Y, minsize, correl){
  # Correlation generation for Statistical Recoupling of Variables.
  # input: data matrix X, class matrix Y, bucketting resolution
  # output: Pfinal:pvalue vector, supercluster boreders: indicesdebf and
  #    indicesfinf; Correlation of superclusters/clusters and residual variables with the Y matrix.

  #Analysis of variance by the one-way anova function of matlab modified to
  #remove all printing options.
  nc<-dim(X)[[2]]
  #Identification of the residual water area
  wat<-which(apply(X,2,sum)==0)
  if (length(wat)>0){
    lsup<-min(wat)-1
    linf<-max(wat)
  }
  
  res<-SRV(X, minsize, correl)
  indicesdebf<-res[[1]]
  indicesfinf<-res[[2]]
  Xclusterf<-res[[3]]

  nct<-dim(Xclusterf)[[2]];

  Pvalue<-rep(0,nct);
  for (j in 1:nct){
  	g<-lm(Xclusterf[,j]~Y[,2])
    Pvalue[j]<-anova(g)$"Pr(>F)"[1];
  }

  #Computation of the Benjamini-Yekutieli correction on the nct-1 clusters (a
  #cluster corresponds to the residual water signal and is removed before analysis).
  I<-order(Pvalue)
  PBY<-rep(0,nct);
  for (j in 1:nct){
    PBY[j]<-Pvalue[j]*(nct-1)*(log(nct-1)+0.5772156649)/I[j];
  }
  if (length(wat)>0){
    PBY[lsup:linf]<-1;
  }
  #Comparison of the computed q-value with the statisitcal threshold of 0.05.
  PvBY<-rep(0,nct);
  for (i in 1:nct){
    if (PBY[i]<0.05){
      PvBY[i]<-2;
    }else{
    	if (PBY[i]>=0.05){
          PvBY[i]<-1;
        }
    }
  }

  #Representaion of the final p-value vector for the initial nc variables
  Pfinal<-rep(0,nc);

  for (i in 1:(indicesdebf[1]-1)){
    Pfinal[i]<-1;
  }

  for (i in 1:(nct-1)){
    for (j in (indicesfinf[i]+1):(indicesdebf[i+1]-1)){
      Pfinal[j]<-1;
    }
  }

  for (i in 1:nct){
    if (PvBY[i]==2){
      for (j in indicesdebf[i]:indicesfinf[i]){
        Pfinal[j]<-2;
      }
    }else{
      for (j in indicesdebf[i]:indicesfinf[i]){
        Pfinal[j]<-1;
      }
    }
  }

  if (indicesfinf[nct]<nc){
  	for (i in (indicesfinf[nct]+1):nc){
      Pfinal[i]<-1;
    }
  }

  #measurement of the correlation of the clusters/superclusters/residual
  #variables with the Y matrix.
  CorrelationB<-rep(0,nc);
  for (j in 1:nc){
    CorrelationB[j]=abs(cor(Xuv[,j],Y[,2]));
  }
  
  Xclusterfuv<-Uvscaled(Xclusterf);
  CorrelationC<-rep(0,nct);
  for (j in 1:nct){
      CorrelationC[j]<-abs(cor(Xclusterfuv[,j],Y[,2]));
  }

  Correlation<-rep(0,nc);
  if (indicesdebf[1]>1){
    for (i in 1:(indicesdebf[1]-1)){
      Correlation[i]<-CorrelationB[i];
    }
  }

  for (i in 1:(nct-1)){
    for (j in (indicesfinf[i]+1):(indicesdebf[i+1]-1)){
      Correlation[i]<-CorrelationB[i];
    }
  }

  for (i in 1:nct){
    for (j in indicesdebf[i]:indicesfinf[i]){
      Correlation[j]<-CorrelationC[i];
    }
  }

  if (indicesfinf[nct]<nc){
    for (i in (indicesfinf[nct]+1):nc){
      Correlation[i]<-CorrelationB[i];
    }
  }
  
  return(list(Pfinal,indicesdebf,indicesfinf,Xclusterf,Correlation))

}

rectangle=function(data){

  n=length(data)
  s=sum(data)
  t=min(data)*n
  if(s-t<0){print(paste(s,t,n,s-t))}
  #s=s-(data[1]+data[n])*n/2
 # print(c(min(data),max(data),s,t))
  return(s-t)
}

trapeze=function(data){
  trapeze <- 0
  n=length(data)
  s=sum(data)
  a=(data[n]-data[1])/n
  b=data[1]-a
  t=sum(pmin(data,1:n*a+b))
  #s=s-(data[1]+data[n])*n/2
 # print(c(min(data),max(data),s,t))
  return(s-t)
}

# Simulation scripts:
add_peak=function(data,gen,loc,amp,eff,s=1,sloc=1,samp=0.1,seff=0.1){
# Script to add a simulated peak
  nr=dim(data)[[1]]
  nc=dim(data)[[2]]
  for(r in 1:nr){
    if(gen[r]==9){
      g=sample(c(0,1,2))[1]
    }else{
      g=gen[r]
    }
   # data[r,]=data[r,]+abs(rnorm(1,amp,samp)+rnorm(1,eff,seff)*g)*dnorm(1:nc,rnorm(1,loc,sloc),s)
    sdamp<-abs(rnorm(1,amp,samp*amp)*(1+rnorm(1,eff,seff*eff))*g)
    data[r,]=data[r,]+sdamp*dnorm(1:nc,rnorm(1,loc,sloc),s)*s*sqrt(2*pi)
  }
#  print(paste(nr,nc,loc,amp,eff))
  return(data)
}

set_baseline=function(nr,nc,b=1){
# Script to create a simulated baseline
  print(c(nr,nc))
  out<-runif(nr*nc,min=0,max=b)
  print(length(out))
  print(dim(out))
  dim(out)<-c(nr,nc)
  print(dim(out))
  return(out)
}

NMR.plot=function(data,ppm,x=5000,t=2000,k=1,ylim=range(data[k,x:(x+t)]),title="NMR profile"){
 # Plot the region of size t, starting at x


  y=-min(data[k,x:(x+t)])/10;
  if (length(k)==1){
    plot(ppm[x:(x+t)],data[k,x:(x+t)],type="l",col=4,ylim=ylim,xlab="NMR resonance",ylab="Intensity",main=title) ; 
  }else{
    plot(ppm[x:(x+t)],data[1,x:(x+t)],type="l",col=1,ylim=ylim,xlab="NMR resonance",ylab="Intensity",main=title) ; 
    j=2
    for (i in k[-1]){
      points(ppm[x:(x+t)],data[i,x:(x+t)],type="l",col=j) ; 
      j=j+1 %% 8
    }
  }
}

SRV.plot=function(data,SRV,x=5000,t=2000,k=1,ylim=range(data[k,x:(x+t)]),title="Cluster plot"){
 # Plot the region of size t, starting at x
#Plot arrows defined by SRV on data


  y=-min(data[k,x:(x+t)])/10;
  if (length(k)==1){
    plot(data[k,x:(x+t)],type="l",col=4,ylim=ylim,xlab="NMR resonance",ylab="Intensity",main=title) ; 
  }else{
    plot(data[1,x:(x+t)],type="l",col=1,ylim=ylim,xlab="NMR resonance",ylab="Intensity",main=title) ; 
    j=2
    for (i in k[-1]){
      points(data[i,x:(x+t)],type="l",col=j) ; 
      j=j+1 %% 8
    }
  }
  abline(v=SRV[[1]]-x,col=2); 
  abline(v=SRV[[2]]-x,col=3) ;
  arrows(SRV[[1]]-x,y,SRV[[2]]-x,y,0.1)

}

arrow.plot=function(file="ur.alig.dat",lo=4.00,hi=4.10,k=1:168,title="Peak calling consensus")
{
  # Plot NMR profile plus SRV regions and consensus across the various statistics


  d1<-read.csv(file)
  f<-dim(d1)[[2]]
  dat<-d1[,-c(f-2,f-1,f)]
  ppm<-as.numeric(sub("e\\.","e-",sub("ppm_","",sub("ppm_\\.","-", colnames(dat)))))
  st<-which(ppm>=lo & ppm<=hi)
  start=min(st)
  size=length(st)
  mm<-max(dat[start:(start+size)])
  NMR.plot(dat,ppm,start,size,k=k,title=title,ylim=c(-mm/3,mm))
  par(new=T)
  n=start:(start+size)*3

  consensus<-read.table("consensus.dat")
  me<-read.table("mean_SRV.ppm",header=T)
  med<-read.table("median_SRV.ppm",header=T)
  ma<-read.table("max_SRV.ppm",header=T)
  su<-read.table("sum_SRV.ppm",header=T)
  tr<-read.table("trapeze_SRV.ppm",header=T)
  re<-read.table("rectangle_SRV.ppm",header=T)
  plot(consensus[n,2],-consensus[n,3],ty="b",pch=3,yaxt='n',ylim=c(-70,210),ylab="",xlab="",axes=F,cex=0.5);
  axis(4,labels=c(64,32,16,8,4,2,1,0),at=c(-64,-32,-16,-8,-4,-2,-1,0))
  mtext("Consensus Score",4)
  arrows(me[,1],rep(-55,dim(me)[[1]]),me[,2],rep(-55,dim(me)[[1]]),col="red",pch=4,length=0.1,angle=15);
  arrows(ma[,1],rep(-50,dim(ma)[[1]]),ma[,2],rep(-50,dim(ma)[[1]]),col="green",pch=3,length=0.1,angle=15);
  arrows(su[,1],rep(-45,dim(su)[[1]]),su[,2],rep(-45,dim(su)[[1]]),col="blue",pch=2,length=0.1,angle=15);
  arrows(med[,1],rep(-40,dim(med)[[1]]),med[,2],rep(-40,dim(med)[[1]]),col=6,pch=6,length=0.1,angle=15);
  arrows(tr[,1],rep(-35,dim(tr)[[1]]),tr[,2],rep(-35,dim(tr)[[1]]),col=7,pch=7,length=0.1,angle=15);
  arrows(re[,1],rep(-30,dim(re)[[1]]),re[,2],rep(-30,dim(re)[[1]]),col=8,length=0.1,angle=15);
}

simple.plot=function(file="ur.alig.dat",lo=4.00,hi=4.10,k=1:168,title="Peak calling consensus")
{
  # Plot NMR profile plus SRV regions 


  d1<-read.csv(file)
  f<-dim(d1)[[2]]
  dat<-d1[,-c(f-2,f-1,f)]
  ppm<-as.numeric(sub("e\\.","e-",sub("ppm_","",sub("ppm_\\.","-", colnames(dat)))))
  st<-which(ppm>=lo & ppm<=hi)
  start=min(st)
  size=length(st)
  NMR.plot(dat,ppm,start,size,k=k,title=title)

  mtext("Overall Profile")
}

ppersp = function (z,ppm,titre,theta=-15, phi=15,r=50)
{
#  z<-lip.SRVlog10_1.female.ehk
  #map<-1:TotM
  map<-1:dim(z)[1]
  jet.colors <- colorRampPalette( c("blue","green","red"))
  color <- jet.colors(100)
  nrz <- nrow(z)
  ncz <- ncol(z)
  zfacet <- z[-1, -1] + z[-1, -ncz] + z[-nrz, -1] + z[-nrz, -ncz]
  facetcol <- cut(zfacet, 100)
  persp(map, as.numeric(ppm), z, theta = theta, phi = phi, r=r, col = color[facetcol], main=titre,ticktype="detailed",nticks=10,  border = NA,xlab="Location",ylab="Shift",zlab="LOD")
}


psave = function(z,p,titre,ppm,res,LT)
{# Save in a file scores of z higher than LT for this permutation p

  zgt<-which(z>LT,arr.ind=T)
  if (nrow(zgt)>0){
    results<-cbind(p,z[zgt],rownames(res)[zgt[,1]],ppm[zgt[,2]])
    colnames(results)<-c("Perm","MaxScore","Marker","Shift")
    write.table(results,paste(titre,"dat",sep="."),row.names=F,col.names=F,quote=F,append=T,sep="\t")
  }
}


running = function (spectrum,width = 0.00025,mid=30,base=10) 
{
        n = length(spectrum);
        out = spectrum;
        for(i in (mid+1):(mid+base)) {
                out <- out - (c(spectrum[(i+1):n],rep(0,i)) + c(rep(0,i),spectrum[1:(n-i)])) * mid/base
        }
        out <- out + 0.5*(c(spectrum[(mid+1):n],rep(0,mid)) + c(rep(0,mid),spectrum[1:(n-mid)]))
        for(i in 1:(mid-1)) {
                out <- out + c(spectrum[(i+1):n],rep(0,i)) + c(rep(0,i),spectrum[1:(n-i)])
        }
        out * width
}

pplot = function(z,titre,ppm,res,LT=c(5,10,15,20))
{# Plot the results with a color scale y layer over 3 in 2D

  nm=dim(res)[1]
  which(z>LT[1],arr.ind=T)-> zgt4
  which(z>LT[2],arr.ind=T)-> zgt5
  which(z>LT[3],arr.ind=T)-> zgt6
  which(z>LT[4],arr.ind=T)-> zgt7
  plot(y=100, x= -1,xlim=range(ppm),ylim=c(1,nm),xlab="Resonance",main=titre,type="n",yaxt="n",ylab="Chrom postion")
  points(x=ppm[zgt4[,2]], y= zgt4[,1],col="green",pch=22,bg="green",cex=0.2)
  points(x=ppm[zgt5[,2]], y= zgt5[,1],col="blue",pch=22,bg="blue",cex=0.2)
  points(x=ppm[zgt6[,2]], y= zgt6[,1],col="red",pch=22,bg="red",cex=0.2)
  points(x=ppm[zgt7[,2]], y= zgt7[,1],col="yellow",pch=22,bg="yellow",cex=0.2)
  midpoints=1:length(levels(res[,1]))
  chrchanges=1:length(levels(res[,1]))+1
  prev=0
  step=0
  start=0
  chr="1"
  k=1
  i=1
  chrchanges[1]=0
  while(i <= nm){
    if(res[i,1]!=chr){
      midpoints[k]=step+(prev-start)/2
      step=prev
      k=k+1
      start=i
      chrchanges[k]=start-0.5
      chr=res[i,1]
    }
    prev=i
    i=i+1
  }
  midpoints[k]=(prev-start)/2+step
  chrchanges[k+1]=nm

  axis(side=2,at=midpoints,labels=levels(res[,1]), tick=F, line=-0.2, lty=6);
  axis(side=2,at=1+chrchanges,labels=F)
  rug(ppm,ticksize=0.01)

}


diff_res = function (res1, res2, nr=1000, nc=200, met=max){
# Function to compare 2 arrays by first binning them into a nr x nc matrix using the method met to summarize the region 

  dr1<-dim(res1)[1]
  dc1<-dim(res1)[2]
  dr2<-dim(res2)[1]
  dc2<-dim(res2)[2]
  if (nr>dr1 || nr>dr2){nr=min(dr1,dr2)}
  if (nc>dc1 || nc>dc2){nc=min(dc1,dc2)}
  res=matrix(0,nrow=nr,ncol=nc)

  print(paste("Setting up",nr,"x",nc,"from",dr1,"x",dc1,"and",dr2,"x",dc2))
  nr1<-c(1,floor(1:nr/nr*dr1))
  nc1<-c(1,floor(1:nc/nc*dc1))
  nr2<-c(1,floor(1:nr/nr*dr2))
  nc2<-c(1,floor(1:nc/nc*dc2))
  for(i in 1:nr){
    for( j in 1:nc){
      res[i,j]<-met(res1[nr1[i]:nr1[i+1],nc1[j]:nc1[j+1]])-met(res2[nr2[i]:nr2[i+1],nc2[j]:nc2[j+1]])
    }
  }

  return(res)
}

# Author : Dr. Vincent Navratil
# Institution :  ENS Lyon
# Date of creation : 28/07/2009
# MSEA function Args
# total: two column matrix or table - first column = compound_id (KEGG, HMDB ...); second column = pathway class (KEGG, HMDB ...)
# sample: a vector of sample's compound_id (KEGG, HMDB ...)
# description: a two column matrix or table - first column = pathway class (KEGG, HMDB ...); second column = pathway description
# method: see p.adjust multiple testing adjustment method

msea = function(total,sample,description,method){

total=unique(total[,c(1,2)])
names(total)=c("id","class")
names(description)=c("class","description")

#prepare total class summary
total=unique(total)
class=total$class
class_summary=as.data.frame.table(sort(table(class)))
names(class_summary)=c("class","n")

total_id_number=length(unique(total$id))

#prepare sample class summary
sample_total=total[which(total$id %in% sample),]
sample_class=sample_total$class
sample_class_summary=as.data.frame.table(sort(table(sample_class)))
names(sample_class_summary)=c("class","n")

sample_id_number=length(unique(sample_total$id))

#prepare contingency table

contingency=merge(class_summary,sample_class_summary,by.x="class",by.y="class")
names(contingency)=c("class","total","sample")

contingency_p_value=unlist(lapply(1:length(contingency$class),function(x){fisher.test(matrix(c(contingency$total[x],total_id_number-contingency$total[x],contingency$sample[x],sample_id_number-contingency$sample[x]),nrow=2),alternative="less")$p.value}))
contingency_p_value_adj=p.adjust(contingency_p_value,method=method)

odds=(contingency$sample/sample_id_number)/(contingency$total/total_id_number)

result=cbind(as.data.frame(contingency),odds,contingency_p_value,contingency_p_value_adj)
names(result)=c("class","total","sample","odds","pval","pval_adjust_BH")

result=merge(result,description,by.x="class",by.y="class")

result
}
# End of MSEA by V. Navratil
####################################################

read.cross.gary = function (dir, genfile, mnamesfile, chridfile, phefile, pnamesfile, mapfile, estimate.map, na.strings) 
		{
		    if (missing(genfile)) genfile <- "geno.dat"
		    if (missing(mnamesfile)) mnamesfile <- "mnames.txt"
		    if (missing(chridfile)) chridfile <- "chrid.dat"
		    if (missing(phefile)) phefile <- "pheno.dat"
		    if (missing(pnamesfile)) pnamesfile <- "pnames.txt"
	            if (missing(mapfile)) mapfile <- "markerpos.txt"
		    if (!missing(dir) && dir != "") {
	              genfile <- file.path(dir, genfile)
	              mnamesfile <- file.path(dir, mnamesfile)
		      chridfile <- file.path(dir, chridfile)
	              if (!is.null(mapfile)) mapfile <- file.path(dir, mapfile)
	          }
	          allgeno <- as.matrix(read.table(genfile, na.strings = "9")) + 1
	          		  pheno <- phefile
	          chr <- scan(chridfile, what = character(), quiet = TRUE)
	          mnames <- scan(mnamesfile, what = character(), quiet = TRUE)
	          if (!is.null(mapfile)) {
		     map <- read.table(mapfile, row.names = 1)
		     map <- map[mnames, 1]
		     map.included <- TRUE
		  }else{
		     map <- seq(0, by = 5, len = length(mnames))
	             map.included <- FALSE
	          }
	          if (!is.null(pnamesfile)) pnames <- pnamesfile
		    else pnames <- paste("pheno", 1:ncol(pheno), sep = "")
	          uchr <- unique(chr)
	          n.chr <- length(uchr)
		  geno <- vector("list", n.chr)
		  names(geno) <- uchr
		  min.mar <- 1
		  for (i in 1:n.chr) {
		     temp.map <- map[chr == uchr[i]]
		     if (any(is.na(temp.map))) {
		        o <- (seq(along = temp.map))[is.na(temp.map)]
		        for (j in o) {
		           if (j == 1 || all(is.na(temp.map[1:(j - 1)]))) {
		               z <- min((seq(along = temp.map))[-o])
		               temp.map[j] <- min(temp.map, na.rm = TRUE) - (z - j + 1)
		           } else if (j == length(temp.map) || all(is.na(temp.map[-(1:j)]))) {
			               z <- max((seq(along = temp.map))[-o])
			               temp.map[j] <- max(temp.map, na.rm = TRUE) + (j - z + 1)
		                  } else {
				       temp.map[j] <- (min(temp.map[-(1:j)], na.rm = TRUE) + max(temp.map[1:(j - 1)], na.rm = TRUE))/2
	                          }
			}
		    }
		    names(temp.map) <- mnames[chr == uchr[i]]
	            data <- allgeno[, min.mar:(length(temp.map) + min.mar - 1), drop = FALSE]
                    min.mar <- min.mar + length(temp.map)
                    colnames(data) <- names(temp.map)
		    geno[[i]] <- list(data = data, map = temp.map)
		    if (uchr[i] == "X" || uchr[i] == "x") class(geno[[i]]) <- "X" else class(geno[[i]]) <- "A"
            }
	    colnames(pheno) <- pnames
            sw2numeric <- function(x) {
	      pattern <- "^[ \t]*-*[0-9]*[.]*[0-9]*[ \t]*$"
              n <- sum(!is.na(x))
              if (length(grep(pattern, as.character(x[!is.na(x)]))) == n) return(as.numeric(as.character(x)))
                else return(x)
            }
	    pheno <- data.frame(lapply(as.data.frame(pheno), sw2numeric))
	    n.mar1 <- sapply(geno, function(a) ncol(a$data)) 
	    n.mar2 <- sapply(geno, function(a) length(a$map))
	    n.phe <- ncol(pheno)
	    n.ind1 <- nrow(pheno)
	    n.ind2 <- sapply(geno, function(a) nrow(a$data))
	    if (any(n.ind1 != n.ind2)) {
	      print(c(n.ind1, n.ind2)) 
	      stop("Number of individuals in genotypes and phenotypes do not match.")
	    }
            if (any(n.mar1 != n.mar2)) {
  	      print(c(n.mar, n.mar2))
	      stop("Numbers of markers in genotypes and marker names files do not match.")
	    }
	    cat(" --Read the following data:\n") 
	    cat("\t", n.ind1, " individuals\n") 
	    cat("\t", sum(n.mar1), " markers\n") 
	    cat("\t", n.phe, " phenotypes\n")
	    if (max(allgeno[!is.na(allgeno)]) <= 2) type <- "bc" else type <- "f2" 
	    cross <- list(geno = geno, pheno = pheno)
	    class(cross) <- c(type, "cross")
	    cross.type <- class(cross)[1] 
	    if (cross.type == "f2") max.gen <- 5 else max.gen <- 2 
	    u <- unique(allgeno) 
	    if (any(!is.na(u) & (u > max.gen | u < 1))) { 
	      err <- paste("There are stange values in the genotype data :", paste(sort(u), collapse = ":"), ".")
	      stop(err)
	    }
	    cross$pheno <- as.data.frame(cross$pheno)
	    if (!map.included) {
  	      for (i in 1:nchr(cross)) cross$geno[[i]]$map <- cross$geno[[i]]$map - min(cross$geno[[i]]$map)
            }
	    list(cross, (!map.included && estimate.map))
}


