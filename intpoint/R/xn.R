xn <-
function(x,C,p,a){
max<-0
for (i in 1:length(x))
if (abs(x[i]*C[i])>max) max<-abs(x[i]*C[i])
div<-a*p/max
s<-x+div
return(s)}
