near<-function(x,pt,fr=1){
# determine which values in x are near pt
# based on fr * mad
m<-mad(x)
if(m==0){
temp<-idealf(x)
m<-(temp$qu-temp$ql)/(qnorm(.75)-qnorm(.25))
}
if(m==0)m<-sqrt(winvar(x)/.4129)
if(m==0)stop("All measures of dispersion are equal to 0")
dis<-abs(x-pt)
dflag<-dis <= fr*m
dflag
}
