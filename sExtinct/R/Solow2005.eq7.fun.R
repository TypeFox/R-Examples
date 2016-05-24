Solow2005.eq7.fun <-
function(dd){
		names(dd)<-c("yrs", "sights")
		dd$year<-dd$yrs-min(dd$yrs)
		sSolow<-sum(dd$year*dd$sights)
		tn<-max(dd[(dd$sights>0),]$year)
		TSolow<-max(dd$year)
		n<-sum(dd$sights)
		y<-tn/sSolow
		dummy<-0
		if(floor(1/y)>0){
			for(i in seq(floor(1/y))){
				dummy<-dummy+(-1)^(i-1)*choose(n,i)*(1-i*y)^(n-1)}}
		Fs1<-1-dummy
		y<-TSolow/sSolow
		dummy<-0
		if(floor(1/y)>0){
			for(i in seq(floor(1/y))){
				dummy<-dummy+(-1)^(i-1)*choose(n,i)*(1-i*y)^(n-1)}}
		Fs2<-1-dummy
		res<-Fs1/Fs2
		return(res)
		}
