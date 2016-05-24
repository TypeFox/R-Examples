TSDistances <- function(x, y, tx=NULL, ty=NULL, distance, ...) {

if (is.character(tx)) {
  distance <- tx 
  tx <- NULL
}

#If x is given as a ts object, the values and the time index are extracted.
if (is.ts(x)) {
  tx<-as.numeric(time(x))
  x<-as.numeric(x)
}  
  
#If y is given as a ts object, the values and the time index are extracted.
if (is.ts(y)) {
  ty<-as.numeric(time(y))
  y<-as.numeric(y)
} 

#If x is given as a xts or a zoo object, the values and the time index are extracted.
if (is.zoo(x) | is.xts(x)) {
  tx<-as.numeric(index(x))
  x<-as.numeric(x)
}  
  
#If y is given as a xts or a zoo object, the values and the time index are extracted.
if (is.zoo(y) | is.xts(y)) {
  ty<-as.numeric(index(y))
  y<-as.numeric(y)
} 

#If x is given as a numerical vector but the time index is not provided
#a constant sampling rate is assumed.
if (is.numeric(x) & is.null(tx)) {
  tx <- c(1:length(x))
} 

#If y is given as a numerical vector but the time index is not provided
#a constant sampling rate is assumed.
if (is.numeric(x) & is.null(ty)) {
  ty <- c(1:length(y))
}

possible.distances <- c("euclidean", "manhattan", "minkowski", "infnorm", "ccor", "sts", "dtw", "lb.keogh", "edr", "erp", "lcss", "fourier", "tquest", "dissim", "acf", "pacf", "ar.lpc.ceps", "ar.mah", "ar.mah.statistic", "ar.mah.pvalue", "ar.pic", "cdm", "cid", "cor", "cort", "wav", "int.per", "per", "mindist.sax", "ncd", "pred", "spec.glk", "spec.isd", "spec.llr", "pdc", "frechet")

distance <- match.arg(distance, possible.distances)

#The distance is calculated.
d<-switch(distance, 
       "euclidean" = EuclideanDistance(x, y),
       "manhattan" = ManhattanDistance(x, y),
       "minkowski" = MinkowskiDistance(x, y, ...), 
       "infnorm" = InfNormDistance(x, y),
       "ccor" = CCorDistance(x, y, ...),
       "sts" = STSDistance(x, y, tx, ty),
       "dtw" = DTWDistance(x, y, ...),
       "lb.keogh" = LBKeoghDistance(x, y, ...),
       "edr" = EDRDistance(x, y, ...),
       "erp" = ERPDistance(x, y, ...),
       "lcss" = LCSSDistance(x, y, ...),
       "fourier" = FourierDistance(x, y, ...),
       "tquest" = TquestDistance(x, y, ...),
       "dissim" = DissimDistance(x, y, tx, ty),
       "acf" = ACFDistance(x, y, ...),
       "pacf" = PACFDistance(x, y, ...),
       "ar.lpc.ceps" = ARLPCCepsDistance(x, y, ...),
       "ar.mah" = ARMahDistance(x, y, ...),
       "ar.mah.statistic" = ARMahStatisticDistance(x, y, ...),
       "ar.mah.pvalue" = ARMahPvalueDistance(x, y, ...),
       "ar.pic" = ARPicDistance(x, y, ...),
       "cdm" = CDMDistance(x, y, ...),
       "cid" = CIDDistance(x, y, ...),
       "cor" = CorDistance(x, y, ...),
       "cort" = CortDistance(x, y, ...),
       "wav" = WavDistance(x, y, ...),
       "int.per" = IntPerDistance(x, y, ...),
       "per" = PerDistance(x, y, ...),
       "mindist.sax" = MindistSaxDistance(x, y, ...),
       "ncd" = NCDDistance(x, y, ...),
       "pred" = PredDistance(x, y, ...),
       "spec.glk" = SpecGLKDistance(x, y, ...),
       "spec.isd" = SpecISDDistance(x, y, ...),
       "spec.llr" = SpecLLRDistance(x, y, ...),
       "pdc" = PDCDistance(x, y, ...),
       "frechet" = FrechetDistance(x, y, tx, ty)
       )
  
return(d)

}


