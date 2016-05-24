inprod.fdata=function (fdata1,fdata2=NULL, w = 1, ...)   {
if (!inherits(fdata1,"fdata")) stop("No fdata class")
tt1<-fdata1[["argvals"]]
DATA1<-fdata1[["data"]]
nas1<-apply(fdata1$data,1,count.na)
rtt<-fdata1[["rangeval"]]
if (any(nas1)) {
   stop("fdata1 contain ",sum(nas1)," curves with some NA value \n")
   }
 if (is.null(fdata2)) {fdata2<-fdata1}
 else  if (!inherits(fdata2,"fdata")) stop("No fdata class")
 nas2<-apply(fdata2$data,1,count.na)
 if (any(nas2)) {
   stop("fdata2 contain ",sum(nas2)," curves with some NA value \n")
   }
 DATA2<-fdata2[["data"]]
 tt2<-fdata2[["argvals"]]
 if  (sum(tt1!=tt2)!=0) stop("Error: different discretization points in the input data.\n")
 numgr = nrow(DATA1); numgr2 = nrow(DATA2)
 dtt<-diff(tt1)
 eps<-as.double(.Machine[[1]]*10)
 inf<-dtt-eps;sup<-dtt+eps
 np<-length(tt1)
 if (all(dtt>inf) & all(dtt<sup)) {equi=TRUE}
 else equi=FALSE
 if ((length(w)!=np) & (length(w) != 1)) {
    stop("DATA ERROR: The weight vector hasn't the length of the functions\n")
 }
 mdist = array(0, dim = c(numgr, numgr2))
 for (i in 1:numgr) {
    for (ii in 1:numgr2) {
         f=w*(DATA1[i,]*DATA2[ii,])
         mdist[i,ii]=int.simpson2(tt1,f,equi,...)
}}
return(mdist)
}

