null.breadth.focal <-
function(dat,dis.method="jaccard",reps=100){
	res<-vector("list",dim(dat)[1])
for(focal.bug in 1:dim(dat)[1]){
	cat("\n",focal.bug,"of",dim(dat)[1])
	dat.null<-dat[-focal.bug,]
	temp<-which(apply(dat.null,2,sum)==0)
		if(length(temp)!=0){
			dat.null<-dat.null[,-temp]
			focal.bug.m.f<-dat[focal.bug,-temp]
				}else{focal.bug.m.f<-dat[focal.bug,]}
		#dat.null is the matrix minus the focal species
		#focal.bug.m.f is the focal vector minus unique plants
		ug.scale.factor<-hyp.ordi.breadth(dat.null,grouping=rep(TRUE,dim(dat.null)[2]))
		observed.breadth<-hyp.ordi.breadth(dat.null,grouping=as.numeric(focal.bug.m.f))
		null.breadth.reps<-NA
			for(i in 1:reps){
				null.breadth.reps[i]<-hyp.ordi.breadth(dat.null,grouping=sample(as.numeric(focal.bug.m.f),size=length(focal.bug.m.f)))}
			res[[focal.bug]]$species<-rownames(dat)[focal.bug]
			res[[focal.bug]]$observed.breadth<-observed.breadth
			res[[focal.bug]]$scale.factor<-ug.scale.factor
			res[[focal.bug]]$observed.breadth.scaled<-observed.breadth/ug.scale.factor
			res[[focal.bug]]$totalplantrichness<-sum(dat[focal.bug,])
			res[[focal.bug]]$modplantrichness<-sum(focal.bug.m.f)
			res[[focal.bug]]$null<-null.breadth.reps
		}
		return(res)
		}
