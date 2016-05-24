tauvar<-function(x,cval=3){
#
# Compute the tau measure of scale as described in
# Yohai and Zamar (JASA, 1988, 83, 406-413).
# The computational method is described in Maronna and Zamar
# (Technometrics, 2002, 44, 307-317)
#  see p. 310
#
x<-elimna(x)
s<-qnorm(.75)*mad(x)
y<-(x-tauloc(x))/s
cvec<-rep(cval,length(x))
W<-apply(cbind(y^2,cvec^2),1,FUN="min")
val<-s^2*sum(W)/length(x)
val
}
