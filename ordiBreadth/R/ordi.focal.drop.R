ordi.focal.drop <-
function(dat,dist.method="jaccard"){	
	res<-list()
	type.c="centroid"
	res<-vector("list",dim(dat)[1])
	spnames<-rownames(dat)
	ODB<-ordi.breadth(dat,dist.method=dist.method)
	for(i in 1:length(spnames)){
		res[[i]]$species<-spnames[i]
		res[[i]]$ODB<-ODB$tot.breadth[i]
		res[[i]]$ODB.scaled<-ODB$scaled.breadth[i]
	}
	
	for(focal.bug in 1:length(spnames)){
		cat("\n",focal.bug,"of",length(spnames))
	dat.m.f<-dat[-focal.bug,]
	temp<-which(apply(dat.m.f,2,sum)==0)
	if(length(temp)!=0){
		dat.m.f<-dat.m.f[,-temp]
		focal.bug.m.f<-dat[focal.bug,-temp]
			}else{focal.bug.m.f<-dat[focal.bug,]}
	specialization.m.f<-ordi.breadth(dat.m.f,dist.method=dist.method)
	
	temp<-hyp.group.dist.m.f(dat.m.f,grouping=as.numeric(focal.bug.m.f),distances=TRUE)
	focal.breadth<-temp$distances.all
	total.breadth<-temp$tot.breath
#		cols<-rep(col[1],length(focal.breadth))
#		cols[which(focal.bug.m.f==1)]<-col[2]
			specdist<-focal.breadth#[which(focal.bug.m.f==1)]
#				cols<-cols[order(specdist)]
				specdisto<-specdist[order(specdist)]
	
	scale.factor<-hyp.group.dist.m.f(dat.m.f,grouping=rep("YES",dim(dat.m.f)[2]))
	
	res[[focal.bug]]$focal.distances<-specdisto
#	res[[focal.bug]]$cols<-cols	
	res[[focal.bug]]$focal.breadth<-total.breadth
	res[[focal.bug]]$focal.scale.factor<-scale.factor
	res[[focal.bug]]$focal.scaled.breadth<-total.breadth/scale.factor
				}
	return(res)			
			}
