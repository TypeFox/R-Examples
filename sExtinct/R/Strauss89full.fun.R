Strauss89full.fun <-
function (dd, alpha){
 		names(dd)<-c("yrs", "sights")
 		dd<-subset(dd, dd$sights>0)
 		alphs<-seq(0.01, 1, 0.01)
		conf.level<-1-alphs
		a<- min(dd$yrs)
		b<- max(dd$yrs)
		R<- b-a
		H<- length(dd$yrs)
		B<- ((1-conf.level)^(-1/(H-1)) - 1)
		rc<- B*R
		res<-data.frame(yrs=rev((b + rc)),chance=1-rev(conf.level))
		return(res)
 }
