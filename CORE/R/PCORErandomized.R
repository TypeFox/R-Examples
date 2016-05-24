PCORErandomized <-
function(procno=1,COREobj,boundaries,nprocs=1,rngoffset=0){
minscore<-max(COREobj$minscore,min(COREobj$coreTable[,"score"]))
set.seed(COREobj$seedme)
COREobj$nshuffle<-COREobj$nshuffle-rngoffset
myshuffles<-COREobj$nshuffle%/%nprocs
shuffleres<-COREobj$nshuffle%%nprocs
shuffleskip<-myshuffles*(procno-1)+min(shuffleres,procno-1)
if(procno<=shuffleres)myshuffles<-myshuffles+1
weightList<-vector(mode="list",length=nrow(COREobj$coreTable))
weight<-COREobj$input[,"weight"]
for(i in 1:nrow(COREobj$coreTable)){
	za<-COREobj$input[COREobj$input[,"chrom"]==COREobj$coreTable[i,"chrom"],,drop=F]
	cigt<-za[,"start"]<=COREobj$coreTable[i,"start"]&
		za[,"end"]>=COREobj$coreTable[i,"end"]
	weight[COREobj$input[,"chrom"]==COREobj$coreTable[i,"chrom"]][cigt]<-0      
	weightList[[i]]<-matrix(ncol=2,
		data=c(which(COREobj$input[,"chrom"]==COREobj$coreTable[i,"chrom"])[cigt],
		weight[COREobj$input[,"chrom"]==COREobj$coreTable[i,"chrom"]][cigt]))
}
simscores<-matrix(ncol=myshuffles,nrow=nrow(COREobj$coreTable))
chrmax<-nrow(boundaries)
advanceRNG(randopt=COREobj$shufflemethod,nrand=shuffleskip+rngoffset,
	nevents=nrow(COREobj$input))
for(shuffle in 1:myshuffles){
if(COREobj$shufflemethod=="SIMPLE")z<-
	cbind(randomEventMoves(COREobj$input[,"end"]-COREobj$input[,"start"]+1,
	boundaries),COREobj$input[,"weight"])
if(COREobj$shufflemethod=="RESCALE")z<-
	cbind(randomRescaledEventMoves(COREobj$input[,c("start","end","chrom",
	"weight")],boundaries),COREobj$input[,"weight"])
dimnames(z)[[2]]<-c("start","end","chrom","weight")
chu<-unique(z[,"chrom"])
chc<-rep(0,length(chu))
for(i in 1:length(chu))chc[i]<-sum(z[z[, "chrom"]==chu[i],"weight"])
scoremat<-matrix(nrow=nrow(COREobj$coreTable),ncol=length(chu),data=0)
for(ich in 1:length(chu)){
	wza<-which(z[,"chrom"]==chu[rev(order(chc))][ich])
	za<-z[wza,c("start","end","weight"),drop=F]
	zaf<-za[za[,"weight"]>0,,drop=F]
	y<-cbind(c(zaf[,"end"]+1,zaf[,"start"]),c(-zaf[,"weight"],zaf[,"weight"]))
	y<-y[order(y[,1]),,drop=F]
	cy2<-cumsum(y[,2])
	mscore<-max(cy2)
	for(i in 1:nrow(COREobj$coreTable)){
		if(mscore<minscore){
			scoremat[i:nrow(COREobj$coreTable),ich]<-0
			break
		}	
		scoremat[i,ich]<-mscore
		if(i==nrow(COREobj$coreTable))break
		whereIwas<-which(wza%in%weightList[[i]][,1])
		if(length(whereIwas)>0){
			za[whereIwas,"weight"]<-0
			zaf<-za[za[,"weight"]>0,,drop=F]
			y<-cbind(c(zaf[,"end"]+1,zaf[,"start"]),c(-zaf[,"weight"],zaf[,"weight"]))
			y<-y[order(y[,1]),,drop=F]
			cy2<-cumsum(y[,2])
			mscore<-max(cy2)
		}
	}
}
simscores[,shuffle]<-apply(scoremat,1,max)
}
return(simscores)
}
