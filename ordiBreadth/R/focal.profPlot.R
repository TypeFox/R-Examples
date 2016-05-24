focal.profPlot <-
function(dat,focal.bug,dist.method="jaccard",col=c("black","red")){
	dat.m.f<-dat[-focal.bug,]
	temp<-which(apply(dat.m.f,2,sum)==0)
	if(length(temp)!=0){
		dat.m.f<-dat.m.f[,-temp]
		focal.bug.m.f<-dat[focal.bug,-temp]
			}else{focal.bug.m.f<-dat[focal.bug,]}
	specialization.m.f<-ordi.breadth(dat.m.f,dist.method=dist.method)
	focal.breadth<-hyp.group.dist.m.f(dat.m.f,grouping=as.numeric(focal.bug.m.f),distances=TRUE)$distances.all
		cols<-rep(col[1],length(focal.breadth))
		cols[which(focal.bug.m.f==1)]<-col[2]
			specdist<-focal.breadth#[which(focal.bug.m.f==1)]
				cols<-cols[order(specdist)]
				specdisto<-specdist[order(specdist)]
		plot(1:length(cols),specdisto,col=cols,pch=19,xlab="",ylab="distance",las=1,main=rownames(dat)[focal.bug])}
