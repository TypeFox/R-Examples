.area.between.curves <-
function(x, f1, f2, xrange=c(0,1))
{
	a<-0.0;
	for(i in 1:length(x)) {
		if(x[i]>=xrange[1] & x[i]<=xrange[2]) {
			if(i==1) {
				lhs<-0
			} else if(x[i-1]<xrange[1]) {
				lhs<-xrange[1]
			} else lhs<-x[i-1];
			if(i==length(x)) {
				rhs<-x[i]
			} else if(x[i+1]>xrange[2]) {
				rhs<-xrange[2];
			} else rhs<-x[i+1];
			a<-a+(f2[i]-f1[i])*(rhs-lhs)/2;
		} else if(i!=1) if(x[i-1]>=xrange[1] & x[i-1]<=xrange[2]) {
			y1<-f1[i-1]+(f1[i]-f1[i-1])*(xrange[2]-x[i-1])/(x[i]-x[i-1])
			y2<-f2[i-1]+(f2[i]-f2[i-1])*(xrange[2]-x[i-1])/(x[i]-x[i-1])
			a<-a+(y2-y1)*(xrange[2]-x[i-1])/2;
		} else if(i!=length(x)) if(x[i+1]>=xrange[1] & x[i+1]<=xrange[2]) {
			y1<-f1[i]+(f1[i+1]-f1[i])*(xrange[1]-x[i])/(x[i+1]-x[i])
			y2<-f2[i]+(f2[i+1]-f2[i])*(xrange[1]-x[i])/(x[i+1]-x[i])

			a<-a+(y2-y1)*(x[i+1]-xrange[1])/2;
		}
	}
	return(a)
}



disparity <- function(phy=NULL, data, index=c("avg.sq", "avg.manhattan", "num.states")){

	if(!is.null(phy)){
		if (!"phylo"%in%class(phy)) stop("Supply 'phy' as a 'phylo' object")
		td<-treedata(phy, data)
		phy=td$phy
		data=td$data
		desc=.cache.descendants(phy)$tips
		nb.tip <- length(td$phy$tip.label)
		nb.node <- td$phy$Nnode
		result<-numeric()

		for(i in 1:nb.node) {
			l<-desc[[nb.tip+i]]
			d<-td$data[phy$tip.label[l],]
			result[i]<-.disparity(d, index)
		}
		names(result)=nb.tip+1:nb.node
		return(result)

	} else {
		return(.disparity(data, index))
	}

}


.disparity <- function(data, index=c("avg.sq", "avg.manhattan", "num.states")){
	disp=match.arg(index, c("avg.sq", "avg.manhattan", "num.states"))

	if(disp=="avg.sq") {
		d<-dist(data, method="euclidean")^2
		r<-mean(d)
	}
	else if(disp=="avg.manhattan") {
		d<-dist(data, method="manhattan")
		r<-mean(d)
	}
	else if(disp=="num.states") {
		f<-function(x) length(unique(x))
		d<-apply(data, 2, f)
		r<-mean(d)
	}
	else r<-0;
	return(r)
}



ratematrix=function(phy, dat){
    .ic.sigma(phy=phy, data=dat)
}

.ic.sigma <-
function(phy, data)
{
	td<-treedata(phy, data, sort=TRUE)

	f<-function(x) pic(x, td$phy)
	ic<-apply(td$data, 2, function(x) {
			  names(x)=rownames(td$data)
			  f(x)
			  })
	r<-crossprod(ic, ic)/nrow(ic)
	return(r)
}


.dtt <-
function(phy, data, disp=c("avg.sq", "avg.manhattan", "num.states")){
	disp=match.arg(disp, c("avg.sq", "avg.manhattan", "num.states"))

	phy$node.label<-NULL
	td<-treedata(phy, data)
	phy2<-td$phy
	phy<-new2old.phylo(td$phy)

	result<-numeric()


	node.depth<-branching.times(phy2);
	stem.depth<-numeric();
	stem.depth[1]<-node.depth[1];
	for(i in 2:phy2$Nnode) {
		anc<-which(as.numeric(phy$edge[,2])==-i)
		stem.depth[i]<-node.depth[names(node.depth)==phy2$edge[anc,1]]
	}

	ltt<-sort(node.depth, decreasing=TRUE)
	node.depth<-node.depth/max(ltt);
	stem.depth<-stem.depth/max(ltt);
	ltt<-ltt/max(ltt);
	if(length(dim(td$data))==2) {
		d<-disparity(phy2, td$data, index=disp)
		result[1]<-d[1]
		for(i in 2:length(ltt)) {
			x<-d[stem.depth>=ltt[i-1]&node.depth<ltt[i-1]]
			if(length(x)==0) result[i]=0
			else result[i]<-mean(x);
		}
		result[length(ltt)+1]<-0;
		if(result[1]>0)
		result<-result/result[1];

	} else {
		if(length(dim(td$data))!=3)
		stop("Error in data");

		for(i in 1:dim(td$data)[3]) {
			pp<-as.matrix(td$data[,,i])
			d<-disparity(phy2, pp, index=disp)
			y<-numeric()

			y[1]<-d[1]
			for(j in 2:length(ltt)) {
				x<-d[stem.depth>=ltt[j-1]&node.depth<ltt[j-1]]
				if(length(x)==0) y[j]=0
				else y[j]<-mean(x);
			}
			y[length(ltt)+1]<-0;
			if(y[1]>0)
			y<-y/y[1];

			result<-cbind(result, y)
		}
	}

	return(result);
}

.dtt.polygon=function(mat, t, alpha=0.05){
	k=nrow(mat)
	dd=c(alpha/2, 1-alpha/2)

	tmp=sapply(1:k, function(x) quantile(mat[x,], probs=dd, na.rm=TRUE))
	yy=c(tmp[1,], rev(tmp[2,]))
	xx=c(t, rev(t))
	return(cbind(x=xx, y=yy))
}

.CI=function(CI=0.95){
    prob = (1 - CI)/2
    return(c(prob, 1-prob))
}

dtt<-function(phy, data, index=c("avg.sq", "avg.manhattan", "num.states"), mdi.range=c(0,1), nsim=0, CI=0.95, plot=TRUE, calculateMDIp=F){

	disp=match.arg(index, c("avg.sq", "avg.manhattan", "num.states"))

	td<-treedata(phy, data)

	dtt.data<-.dtt(td$phy, td$data, disp=disp)
	ltt<-sort(branching.times(td$phy), decreasing=TRUE)
	ltt<-c(0, (max(ltt)-ltt)/max(ltt));

	s<-ratematrix(td$phy, td$data)
	dtt.sims=NULL
	MDI=NULL

    ylim=c(range(pretty(dtt.data)))

	if(is.numeric(nsim)){
		if(nsim>0){
			sims<-sim.char(td$phy, s, nsim)

			dtt.sims<-.dtt(td$phy, sims)
			mean.sims<-apply(dtt.sims, 1, mean)
            median.sims<-apply(dtt.sims, 1, median)


			MDI<-unname(.area.between.curves(ltt, apply(dtt.sims, 1, median), dtt.data, sort(mdi.range)))
            names(MDI)=disp
			colnames(dtt.sims)=NULL

            yy=range(dtt.sims)
            ylim=range(c(ylim, yy))
		}
	}

	if(plot){
		plot(ltt, dtt.data, xlab="relative time", ylab="disparity", ylim=ylim, bty="n", type="n")

		if(!is.null(dtt.sims)){
			poly=.dtt.polygon(dtt.sims, ltt, alpha=1-CI)
			polygon(poly[,"x"], poly[,"y"], col=.transparency("lightgray", 0.5), border=NA)
			lines(ltt, median.sims, lty=2)
		}
        lines(ltt, dtt.data, type="l", lwd=2)

	}

	res=list(dtt=dtt.data, times=ltt, sim=dtt.sims, MDI=MDI)
	drp=sapply(res, function(x) is.null(x))
	if(any(drp)) res=res[-which(drp)]

	if(calculateMDIp) {
		pVal<-getMDIp(res)
		res<-c(res, MDIpVal=pVal)
	}
	return(res)
}

getMDIp<-function(dttRes) {
	foo<-function(x) {
		return(.area.between.curves(x= dttRes$times, f1=x, f2=dttRes$dtt))
	}
	mdis<-apply(dttRes$sim,1,foo)
	pVal<-length(which(mdis>=0))/length(mdis)
	return(pVal)
}
