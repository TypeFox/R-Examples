Solow1993.eq2.fun <-
function(dd){
		names(dd)<-c("yrs", "sights")
		Tmin <- min(dd$yrs) 
		Tn <- max(dd[(dd$sights>0),]$yrs)-Tmin 
		Tmax <- max(dd$yrs)-Tmin
		n<-sum(na.omit(dd$sights))
		res<-(Tn/Tmax)^n
		return(res)
		}
