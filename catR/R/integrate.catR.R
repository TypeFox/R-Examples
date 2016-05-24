integrate.catR<-function(x,y){
hauteur<-x[2:length(x)]-x[1:(length(x)-1)]
base<-apply(cbind(y[1:(length(y)-1)],y[2:length(y)]),1,mean)
res<-sum(base*hauteur)
return(res)}

