bd.km=function(phy=NULL, time, n, missing = 0, crown=TRUE){
	.rate.estimate(phy=phy, time=time, n=n, epsilon=Inf, missing=missing, crown=crown, type="kendall-moran")
}

bd.ms=function(phy=NULL, time, n, missing = 0, crown=TRUE, epsilon=0){
	.rate.estimate(phy=phy, time=time, n=n, epsilon=epsilon, missing=missing, crown=crown, type="magallon-sanderson")
}

.rate.estimate <-
function(phy=NULL, time, n, epsilon = 0, missing = 0, crown=TRUE, type=c("magallon-sanderson", "kendall-moran"))
{
	type=match.arg(type, c("magallon-sanderson", "kendall-moran"))
	
	if(!is.null(phy)) {
		if (class(phy) != "phylo") 
		stop("'phy' must be a 'phylo' object")
       	
		if(is.null(phy$root.edge)) phy$root.edge=0
		
		if(type=="kendall-moran")
		{
			if(missing != 0)
			warning("Current implementation of Kendall-Moran estimate does not account for missing taxa")
			return(.kendallmoran.rate(phy))
		}
		if(!missing(time)) warning("'time' argument has been ignored")
		if(!missing(n)) warning("'n' argument has been ignored")
		
    	time<-max(branching.times(phy))+phy$root.edge
		n<-length(phy$tip.label)
		
		n<-n+missing;
		ctmp=ifelse(phy$root.edge==0, TRUE, FALSE)
		if(crown!=ctmp) {
			if(crown){
				warning(paste("Assuming 'crown' is ", ctmp, " due to non-zero 'root.edge' in 'phy'", sep=""))
			} else {
				warning(paste("Assuming 'crown' is ", ctmp, " due to missing 'root.edge' in 'phy'", sep=""))
			}
		}
		crown=ctmp
	}
	
    if(crown==TRUE) {
    	if(epsilon==0) {
			rate=(log(n)-log(2))/time
    	} else {
 			rate=1/time*(log(n/2*(1-epsilon^2)+
							 2*epsilon+1/2*(1-epsilon)*sqrt(n*(n*epsilon^2-
															   8*epsilon+2*n*epsilon+n)))-log(2))
    	}
		
    } else {
    	if(epsilon==0) {
			rate=log(n)/time
    	} else {
 			rate=1/time*log(n*(1-epsilon)+epsilon)
    	}
		
		
    }
	
	return(rate)
}

.kendallmoran.rate <-
function(phy)
{
	s<-sum(phy$edge.length)
	rate<-(length(phy$tip.label)-2)/s
	return(rate)
}

crown.p <- 
function(phy=NULL, time, n, r, epsilon)
{
	if(!is.null(phy)){
		phy$root.edge=0
		if(!missing(time)) warning("'time' argument has been ignored")
		time=max(heights.phylo(phy))
		if(!missing(n)) warning("'n' argument has been ignored")
		n=Ntip(phy)
	}
	
	b<-((exp((r*time)))-1)/((exp((r*time)))-epsilon)
	a<-epsilon*b
	p<-(((b^(n-2)))*((n*(1-a-b+(a*b)))+a+(2*b)-1))/(1+a)
	return(p)
}

stem.p <-
function(phy=NULL, time, n, r, epsilon)
{
	if(!is.null(phy)){
		if(is.null(phy$root.edge)) {
            warning("'phy' is assumed to have a 'root.edge' element of zero")
            phy$root.edge=0
        }
		if(!missing(time)) warning("'time' argument has been ignored")
		time=max(heights.phylo(phy))
		if(!missing(n)) warning("'n' argument has been ignored")
		n=Ntip(phy)
	}
	
	b<-((exp((r*time)))-1)/((exp((r*time)))-epsilon)
	p<-(b^(n-1))
	return(p)
}

stem.limits <-
function(time, r, epsilon, CI=0.95)
{
    q=1-CI
    prob=c(q/2, 1-(q/2))
	
	limits <- matrix(nrow=length(time), ncol=2)
	for (i in 1:length(time)) {
		beta <- (exp(r*time[i])-1)/(exp(r*time[i]) - epsilon) #From M&S 01 2b
		alpha <- epsilon * beta #From M&S 01 2a
		u <- (log(beta) + log(prob[1]))/log(beta) #From M&S 01 10a
		l <- (log(beta) + log(prob[2]))/log(beta)
		limits[i, 1] <- l
		limits[i, 2] <- u
	}
	rownames(limits)=time
	colnames(limits)=c("lb", "ub")
	return(limits)	
}



crown.limits <-  function(time, r, epsilon, CI=0.95)
{
    q=1-CI
    prob=c(q/2, 1-(q/2))
	
	limits <- matrix(nrow=length(time), ncol=2)
	for (i in 1:length(time))
	{
		foo<-function(x, prob)
		(crown.p(time=time[i], r=r, epsilon=epsilon, n=x)-prob)^2
		u<-optim(foo, par=exp(r*time), prob=prob[1], method="L-BFGS-B")
		l<-optim(foo, par=exp(r*time), prob=prob[2], method="L-BFGS-B")
		
		limits[i, 1] <- l$par
		limits[i, 2] <- u$par
	}
	rownames(limits)=time
	colnames(limits)=c("lb", "ub")
	return(limits)
}


rc <-
function(phy, plot=TRUE, ...)
{
	
	
	nb.tip <- length(phy$tip.label)
	nb.node <- phy$Nnode
    
	nd<-branching.times(phy);
	node.depth<-c(nd, rep(0, nb.tip))
	names(node.depth)<-c(names(nd), as.character(1:nb.tip))
	p<-numeric(nb.node)
	max.desc<-numeric(nb.node)
	num.anc<-numeric(nb.node)
	node.name<-character(nb.node)
	
	stem.depth<-numeric();
	stem.depth[1]<-node.depth[1];
	
	for(i in 2:length(node.depth)) {
		anc<-which(phy$edge[,2]==names(node.depth)[i])
		stem.depth[i]<-node.depth[names(node.depth)==phy$edge[anc,1]]
	}
	
	
	
	ltt<-sort(nd, decreasing=TRUE)
	
	for(i in 2:length(ltt)) {
		nn<-stem.depth>=ltt[i-1]&node.depth<ltt[i-1]
		anc<-sum(nn)
		desc<-numeric(anc)
		pp<-numeric(anc)
		num.anc[i]<-anc
		for(j in 1:anc) {
			desc[j]<-length(tips(phy, as.numeric(names(nn)[nn][j])))
		}
		max.desc[i]<-max(desc)
		p[i]<-.rcp(max.desc[i], nb.tip, anc)
		node.name[i]<-names(ltt[i])
	}
	num.anc[1]<-1
	max.desc[1]<-nb.tip
	p[1]<-1
	bonf.p<-pmin(p*length(ltt),1)
	res<-cbind(num.anc, max.desc, p, bonf.p)
	rownames(res)<-c("root", node.name[2:nb.node])
	
	if(plot) {
        plotter=function(bonf=FALSE, p.cutoff=0.05, ...){

            labels<-character(length(phy$edge.length))
            names(labels)<-as.character(1:length(labels)+nb.tip)
            if(bonf) {
                s<-which(res[,4]<p.cutoff)
            } else {
                s<-which(res[,3]<p.cutoff)
            }
            mark<-character(length(s))
            if(length(s)>0) {
                for(i in 1:length(s)) {
                    xx<-names(s)[i]
                    tt<-which(ltt==ltt[xx])
                    nn<-stem.depth>=ltt[tt-1]&node.depth<ltt[tt-1]
                    anc<-sum(nn)
                    desc<-numeric(anc)
                    for(j in 1:anc) {
                        desc[j]<-length(tips(phy, names(nn)[nn][j]))
                    }
                    bigone<-which(desc==max(desc))
                    mark[i]<-names(nn)[nn][bigone]
                }
            }
            labels[mark]<-"*"
            phy$node.label<-labels
            plot.phylo(phy, show.node.label=TRUE, no.margin=TRUE, ...)
        }
        plotter(...)
	}
	return(res)
}

.rcp <-
function(ni, n, k) 
{
	max<-floor((n-k)/(ni-1))
	sum=0
	for(v in 0:max) {
		term=(-1)^v*choose(k, v)*choose(n-v*(ni-1)-1, k-1)
		sum=sum+term
 	}
	return(1-(sum/choose(n-1, k-1)))
}


