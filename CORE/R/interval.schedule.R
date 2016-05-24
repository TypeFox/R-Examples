interval.schedule <-
function(v){
	if(!("start"%in%dimnames(v)[[2]])|!("end"%in%dimnames(v)[[2]])|
		length(dim(v))!=2)stop("interval.schedule: incorrect form of input\n")
	lv<-nrow(v)
	if(sum(sort(v[,"start"])==v[,"start"])!=lv)
			stop("interval.schedule: input incorrectly ordered\n")
	unsch<-lv
	schedule<-0
	vsched<-rep(0,lv)
	while(unsch>0){
		schedule<-schedule+1
		vsched[cummax(v[,"end"]*(vsched==0))==v[,"end"]]<-schedule
		unsch<-sum(vsched==0)
	}
	return(vsched)
}
