mest <- function(x, bend = 1.28, na.rm = FALSE){
#
#  Compute M-estimator of location using Huber's Psi.
#  The default bending constant is 1.28
#
if(na.rm)x<-x[!is.na(x)]
if(mad(x)==0)stop("MAD=0. The M-estimator cannot be computed.")
y<-(x-median(x))/mad(x)  #mad in splus is madn in the book.
A<-sum(hpsi(y,bend))
B<-length(x[abs(y)<=bend])
mest<-median(x)+mad(x)*A/B
repeat{
y<-(x-mest)/mad(x)
A<-sum(hpsi(y,bend))
B<-length(x[abs(y)<=bend])
newmest<-mest+mad(x)*A/B
if(abs(newmest-mest) <.0001)break
mest<-newmest
}
mest
}