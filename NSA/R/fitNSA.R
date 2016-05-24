###########################################################################/**
# @set "class=matrix"
# @RdocMethod fitNSA
# @alias fitNSA
# 
# @title "Find normal regions"
#
# \description{
#  @get "title".
# }
#
# @synopsis
#
# \arguments{
#  \item{data}{An Jx2 @numeric @array containing the allele specific copy number values of a number of J SNPs.}
#  \item{...}{Additional argument "chromosome" which indicates the chromosome to which the SNPs correspond to.
#   It is passed to internal NSA.}
# }
#
# \value{
#   Returns an Jx2 @numeric @array with the values of the level of hetezorigosity in one column and NA in the other.
# }                                   
#
#*/###########################################################################

setMethodS3("fitNSA", "matrix", function(data, ...) {
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Validate arguments
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Argument 'data':
  
  nbrOfDataPoints <- nrow(data);

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Fit segmentation model
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  fit <- fitOneNSA(data, ...);
  rm(data); # Not needed anymore

  verbose && enter(verbose, "Class of fitted object: ", class(fit)[1]);
  verbose && exit(verbose);

  fit;
}, private=TRUE)

fitOneNSA <- function(data, chromosome=NULL) 
{
  dataA <- data[,1]*(1-data[,2]);
  dataB <- data[,1]*data[,2]
  #LH Matrix
  LH <- LHMatrix(dataA,dataB);
  ans<-hist(LH,100,plot=FALSE);
  count<-ans$counts;
  mid<-ans$mids;
  a<-mid[order(-count)[1]];
  
  #Threshold1
  Thr1 <- 1/6*a+5/6;

  #LH Selection
  LH0<-LH;

  #Segmentation
  CNA.object<-CNA(genomdat=1*(LH0>=Thr1),array(chromosome,length(LH0)),maploc=1:length(LH0),data.type="binary",sampleid="sample");
  segment.smoothed.CNA.object <- segment(CNA.object, verbose =3, alpha = 1e-6);
  seg<-segment.smoothed.CNA.object$output;
  ini<-seg$loc.start;
  final<-seg$loc.end;
  value<-seg$seg.mean;

  LH0Seg<-array(0,length(LH0));
  for (k in 1:length(ini))
  {
    if (k==1)
      j=1;
    LH0Seg[ini[k]:final[k]]<-value[k];
    if (final[k]==length(LH0))
      j=j+1;
  }

  #Matrices Reconstruction
  LHSeg<-LH0Seg;

  Thr2=0.135;
  MAB<-(LHSeg>=Thr2);
  fit <- data;  
  fit[,1] <- LHSeg;
  fit[,2] <- NA;
  fit;
} # fitOne()

#Create LH Matrix
LHMatrix <- function(A,B)
{
  AB<-A+B;
  LH <- 2*pmin(A,B)/AB;
  return(LH);
}

############################################################################
# HISTORY:
# 2010-06-23 [MO]
# o Created.
############################################################################
