Robson1964.fun <-
function(dd, alpha, data.out){			
	if(data.out==T){
		alpha<-seq(0.02, 1, 0.01)
	}
	else{alpha<-alpha}
	names(dd)<-c("yrs", "sights")
	dd<-subset(dd, dd$sights>0)
	Tmin<-min(dd$yrs)
	Tmax<-max(dd$yrs)-Tmin
	Tn<-dd$yrs[length(dd$yrs)]-Tmin
	res<-Tmin+(Tn+((1-alpha)/alpha)*(Tn-(Tn-1)))
	if(data.out==T){
		res<-data.frame(yrs=rev(res), chance=rev(alpha))
		return(res)}
	else{return(res)}	}
