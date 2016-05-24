optimal <-
function(C,x,e,c){
max<-0
for (i in 1:length(x))
if (x[i]*C[i]>max) max<-x[i]*C[i]
min<-max
for(i in 1:length(C))
if (C[i]<min) min<-C[i]
l<-length(c)
if (max-min*mean(x)<e/l) sol<-1 else sol<-0
return(sol)}
