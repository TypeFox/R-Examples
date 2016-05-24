tauloc<-function(x,cval=4.5){
#
# Compute the tau measure of location as described in
# Yohai and Zamar (JASA, 83, 406-413).
#
x<-elimna(x)
s<-qnorm(.75)*mad(x)
y<-(x-median(x))/s
W<-(1-(y/cval)^2)^2
flag<-(abs(W)>cval)
W[flag]<-0
val<-sum(W*x)/sum(W)
val
}
