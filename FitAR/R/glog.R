`glog` <-
function(x, a=1, InverseQ=FALSE) {
#see glog.nb for derivation of inverse
if (InverseQ) {
    out<-0.25*exp(-x)*(4*exp(2*x)-(a*a))
}
else
    out<-log((x + sqrt(x^2 + a^2))/2)
out
}

