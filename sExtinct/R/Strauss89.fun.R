Strauss89.fun <-
function(dd, alpha){
		names(dd)<-c("yrs", "sights")
 		dd<-subset(dd, dd$sights>0)
		conf.level<-1-alpha		
		a<- min(dd$yrs)
		b<- max(dd$yrs)
		R<- b-a
		H<- length(dd$yrs)
		B<- ((1-conf.level)^(-1/(H-1)) - 1)
		rc<- B*R
		res<-data.frame(Estimate=(b + rc))
		class(res)<-"simpextmod"
		return(res)}
