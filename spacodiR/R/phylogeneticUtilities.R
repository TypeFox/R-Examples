#general statistical function for comparing two vectors on non-independent samples
#author: JM EASTMAN 2010; updated 03.06.2011 to remove NA

resamp.test <-
function(obs=obs, exp=exp,  mu=0, iter=10000, two.tailed=FALSE, na.rm=TRUE){
	if(na.rm) {
		oe=lapply(list(obs,exp), function(x) return(x[!is.na(x)]))
		obs=oe[[1]]
		exp=oe[[2]]
	}
	O=as.numeric(obs)[sample(1:length(obs), iter, replace=TRUE)]
	E=as.numeric(exp)[sample(1:length(exp), iter, replace=TRUE)]
	
	result=c(O-(E+mu))
	p=round(sum(result>=0)/iter, digits=nchar(iter))
	q=1-p
	
	if(two.tailed) {
		res=list(diffs=result, 2*(c(p,q)[which(c(p,q)==min(c(p,q)))]))
	} else {
		res=list(diffs=result, greater=p, lesser=q)
	}
	return(res)
}

#general phylogenetic utility for returning the first (usually, unless a polytomy exists) two descendants the supplied node 
#author: JM EASTMAN 2010

get.desc.of.node <-
function(node, phy) {
	return(phy$edge[which(phy$edge[,1]==node),2])
}


#general phylogenetic utility for returning all descendants (listed as given in phy$edge[,2]) of a node (excluding the supplied node) 
#author: JM EASTMAN 2010, based on GEIGER:::node.leaves()

get.descendants.of.node <-
function (node, phy, tips=FALSE) 
{
    node <- as.numeric(node)
    n <- Ntip(phy)
	if(node <= n) return(NULL)
    l <- c()
    d <- get.desc.of.node(node, phy)
    for (j in d) {
		l <- c(l, j)
        if (j > n) {
			l <- c(l, get.descendants.of.node(j, phy))
		}
    }
    if(tips) return(l[l<=n]) else return(l)
}


#determines if a value falls within a range [min,max]
#author: JM Eastman 2010

withinrange <-
function(x, min, max) {			 
	a=sign(x-min)
	b=sign(x-max)
	if(abs(a+b)==2) return(FALSE) else return(TRUE)
}




######
phy.deresolve=function(phy, time.range=c(0,0), relative=TRUE)
{
    if (is.null(phy$edge.length)) 
	stop("The tree does not appear to have branch lengths")
	if (class(phy)!="phylo") 
	stop("The tree does not appear to be a valid phylo object")
	if(length(time.range)>2)
	stop("Cannot interpret the range of time with more than two elements")
	if(length(time.range)==1)
	time.range=c(0,time.range)
	time.range=time.range[order(time.range)]
	bb <- branching.times(phy)
	if(relative) bb=bb/max(bb)
	inr=sapply(bb, function(x) withinrange(x, time.range[1], time.range[2]))
    ind <- as.numeric(names(bb[inr]))
	if(any(ind==Ntip(phy)+1)) ind=ind[-which(ind==Ntip(phy)+1)]
    n <- length(ind)
    if (!n) {
        return(phy)
	} else {
		ind.tmp = match(ind, phy$edge[,2])
		ind = ind.tmp[!is.na(ind.tmp)]
	}
	orig.edge=phy$edge
	orig.phy=phy
	ntips=Ntip(phy)
    reedge <- function(ancestor, des.to.drop) {
        wh <- which(phy$edge[, 1] == des.to.drop)
		dd <- which(orig.edge[, 2] == des.to.drop)
		dropped.branch <- phy$edge.length[dd]
		d.d <- c(get.desc.of.node(des.to.drop, orig.phy))
		if(length(d.d)) phy$edge.length[match(d.d, orig.edge[,2])]<<-phy$edge.length[match(d.d, orig.edge[,2])]+dropped.branch
		
        for (k in wh) {
            if (phy$edge[k, 2] %in% node.to.drop) {
                reedge(ancestor, phy$edge[k, 2])
            } else {
				phy$edge[k, 1] <<- ancestor
				
			}
        }
    }
	
    node.to.drop <- phy$edge[ind, 2]
    anc <- phy$edge[ind, 1]
    for (i in 1:n) {
        if (anc[i] %in% node.to.drop) 
		next
        reedge(anc[i], node.to.drop[i])
    }
    phy$edge <- phy$edge[-ind, ]
    phy$edge.length <- phy$edge.length[-ind]
    phy$Nnode <- phy$Nnode - n
    sel <- phy$edge > min(node.to.drop)
    for (i in which(sel)) phy$edge[i] <- phy$edge[i] - sum(node.to.drop < phy$edge[i])
    if (!is.null(phy$node.label)) 
	phy$node.label <- phy$node.label[-(node.to.drop - length(phy$tip.label))]
    phy
}

######
phy.nodetimes <-
function(phy, time.range=c(0,0), proportion=TRUE) {
	N=Ntip(phy)
	phy$node.label=NULL
	xx = c(rep(0,N), branching.times(phy))
	names(xx)[1:N]=1:N
	if(proportion) xx=xx/max(branching.times(phy))
	tt=sapply(xx,function(x) withinrange(x,min(time.range),max(time.range)))
	return(xx[tt])
}








