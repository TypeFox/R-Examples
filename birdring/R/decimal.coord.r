###################################################################################
# Computes decimal coordinates from degrees and minutes
# x: vector containing degrees and minutes
# Author: Fraenzi Korner, 9.9.2004, www.oikostat.ch and www.vogelwarte.ch
###################################################################################
decimal.coord<-function(x){
t.sign<-sign(x)
t.minuten<-abs(x)-abs(trunc(x))
t.dezimalstellen<-t.minuten*(5/3)
t.dec<-t.sign*(abs(trunc(x))+t.dezimalstellen)
t.dec
}
###################################################################################
