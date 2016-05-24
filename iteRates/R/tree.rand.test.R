tree.rand.test <- 
function(tree,reps=1000,mod.id=c(1,0,0,0),trace=TRUE){
	obs.cs<-comp.subs(tree,mod.id=mod.id,verbose=FALSE)#observed fit of comp.subs
	obsbt<-branching.times(tree)#observed branching times
	obs.p<-obs.cs$p.val
	obs.er<-obs.cs$Ev
	NullTrees<-treeSims(obsbt,reps)
	ncs<-Null.comp.subs.ef(NullTrees,mod.id=mod.id,trace=trace)
	obs.FPR<-sum(obs.p<0.05,na.rm=TRUE)/sum(!is.na(obs.p))
	FPR.q<-mean(ncs$FPR>obs.FPR)
	list(tree=tree,obs.p=obs.p,ncs=ncs,obs.detection=obs.FPR,p.detection=FPR.q)
}
