dsil<- function(ddist,o.hclr,o.relgr) {
# function introduced in dave 1.5
# silhouette plot by Otto rather than silhouette
	names<- as.character(o.hclr$labels)
	names<- strtrim(names,15)
	n<- length(o.relgr)
	ng<- length(table(o.relgr))
	dgr<- rep(1,ng)
	tgr<- seq(1,ng,1)
	meansilwidth<- rep(0,ng)
	swidth<- rep(0,n)
	second<- rep(0,n)
	large<- 10e10
	mddr<- as.matrix(ddist,nrow=n)
	for (i in 1:n) {
		mi<- mddr[i,]
#   length of current group i, then find ai
		lg<- length(o.relgr[o.relgr == o.relgr[i]])-1
		ai<- sum(mddr[i,o.relgr == o.relgr[i]])/lg
#   array dgr contains group distances to mi
		for (j in 1:ng) dgr[j]<- mean(mi[o.relgr == j])
		dgr[tgr == o.relgr[i]]<- large
		second[i]<- which.min(dgr)
		bi <- dgr[second[i]]
#   silhouette width
		swidth[i]<- (bi-ai)/max(ai,bi)
	}
# mean silhouette width
	meansilwidth<- rep(0,ng)
	for (i in 1:ng) meansilwidth[i]<- mean(swidth[o.relgr == i])
	o.sil<- list(o.relgr=o.relgr,second=second,swidth=swidth,meansilwidth=meansilwidth,names=names)
}
