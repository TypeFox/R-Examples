## PRIOR DISTRIBUTIONS ##
# flexible generation of standardized probability masses for counts variable
dcount=function(x, FUN, ...){
	## current choices for FUN are dlunif (log-uniform), dtpois (truncated poisson), and other built-in functions for discrete distributions (dunif, dgeom, dbinom)
	y=FUN(x, log=TRUE, ...)
	attr(y, "count")=x
	attr(y, "cumsum")=cumsum(exp(y))
	y=.check.prior(y)
	yc=attr(y,"count")
	f=function(x) sapply(x, function(xi) if(xi%in%yc) y[which(yc==xi)] else NA)
	attr(f, "density")=y
	f
}

## PRIOR DISTRIBUTIONS ##
# random variate draw from normalized probability mass for counts variable
.rcount=function(n, FUN){
	FUN=.check.prior(FUN, count=TRUE)
	d=attr(FUN, "density")
	cs=attr(d, "cumsum")
	ct=attr(d, "count")
	r=runif(n)
	sapply(r, function(rv) ct[min(which(rv<cs))])
}

## PRIOR RATIO ##
# FUN: a function whose argument is 'x'; returns log-likelihood
# e.g., FUN=function(x) dexp(x, rate=1/10, log=TRUE)
# see 'controller()' for examples using counts variables
.dlnratio=function(cur, new, FUN){
	sum(FUN(new))-sum(FUN(cur))
}


.gbm.prior.lik=function(rates, jumpvar, njump, nshift, root, control){
	if(njump==0) jumpvar=0
	pR=sum(control$dlnRATE(rates))
	pJ=sum(control$dlnPULS(jumpvar))
	pNj=sum(control$dlnJUMP(njump))
	pNr=sum(control$dlnSHIFT(nshift))
	pT=sum(control$dlnROOT(root))
	prior=sum(c(pR,pJ,pNj,pNr,pT))
	prior
}


## PRIOR DISTRIBUTIONS ##
# standardize probability mass to sum to 1
.standardize.prior=function(x) {
	ct=attr(x,"count")
	xp=exp(x)
	sft=log(1/sum(xp))
	x=sft+x
	xp=xp*(1/sum(xp))
	attr(x,"count")=ct
	attr(x,"cumsum")=cumsum(xp)
	x
}

## PRIOR DISTRIBUTIONS ##
# ensure that prior distribution is proper for auteur
.check.prior=function(x, count=TRUE) {
	if(!count) {
		if(is.function(x)) x else stop("'x' must be a function returning Ln(probabilities).")
	} else {
		if(is.function(x)){
			y=x
			if(is.null(x<-attr(x,"density"))) stop("'x' must have a 'density' attribute.")
			x=.standardize.prior(x)
			attr(y,"density")=x
			x=y
		} else {
			if(is.null(attr(x,"count"))) stop("'x' must have a 'count' attribute.")
			if(is.null(attr(x,"cumsum"))) stop("'x' must have a 'cumsum' attribute.")
			if(!any(x<0)) stop("'x' does not appear to be a vector of Ln(probabilities).")
			if(abs(max(attr(x, "cumsum"))-1)>.Machine$double.eps) {
				x=.standardize.prior(x)
			}
		}
	}
    class(x)=c("gprior", class(x))
    x
}

## PROPOSAL UTILITY ##
.check.lim=function(x, lim=list(min=0,max=Inf), ...){
    
    ff=function(at.bound=TRUE){
        return(at.bound)
    }
    bb=ff(...)
    if(bb){
        if(all(x>=lim$min) & all(x<=lim$max)) return(TRUE) else return(FALSE)
    } else {
        if(all(x>lim$min) & all(x<lim$max)) return(TRUE) else return(FALSE)
    }
}


## PROPOSAL UTILITY ##
#scales phylogeny by relative evolutionary rates under Brownian motion
#author: JM EASTMAN 2011
.scale.brlen <-
function(phy, scalars){
	phy$edge.length=phy$edge.length*scalars
	return(phy)
}

.dunifn=function(n){
	FUN=function(x){
		if(x!=n) return(NA) else return(log(1))
	}
	
	g=0
	attr(g,"count")=n
	attr(g,"cumsum")=1
	attr(FUN,"density")=g
	return(FUN)
}

## PRIOR DISTRIBUTIONS ##
# truncated poisson probability mass vector
dtpois=function(x, min, max, log=TRUE, ...){
		
	if(any(x>max)){
		max=max(x)
		warning("'max' is inconsistent with 'x'")
	}
		
	if(any(x<min)){
		min=min(x)
		warning("'min' is inconsistent with 'x'")
	}
	
	y=dpois(min:max, log=TRUE, ...)
	yp=exp(y)
	if(log){
		sft=log(1/sum(yp))
		y=sft+y
		if(any(is.infinite(y))) warning("Probability mass for some 'x' is effectively zero.")

	} else {
		if(any(yp==0)) warning("Probability mass for some 'x' is effectively zero.")
		y=yp*(1/sum(yp))
	}
	
	nm=min:max
	z=y[match(x, nm)]

	return(z)
}

## PRIOR DISTRIBUTIONS ##
# zero-added log-uniform probability mass vector
dlunif=function(x, min=1, max, log=TRUE, dzero=NULL) {
#	dzero: density at zero
#	log: whether to use log density
#	max: upper limit on the range of integers
#	min: typically 1 unless 'x' does not involve 0
		
	if(any(x>max)){
		max=max(x)
		warning("'max' is inconsistent with 'x'")
	}
	
	if(min<1){
		tmp=max(c(1,min(x)))
		if(tmp!=min) {
			if(min!=0) warning("'min' is inconsistent with 'x'")
			min=tmp
		}
	}
	
	if(any(x==0) && !is.numeric(dzero)) stop("Probability mass at zero must be supplied.")
	
	const=log(max)-log(min)
	xx=min:max
	y=1/(xx*const)
	names(y)=xx
	
	if(is.numeric(dzero) && dzero>0){
		sumd=sum(y)
		dz=dzero*sumd
		names(dz)=0
		y=c(dz,y)
	}
	
	
	s=sum(y)
	y=y*(1/s)
	z=y[match(x, names(y))]
	if(!all(x%in%as.integer(names(z)))) stop("'x' must either be within range of 'min' and 'max' or be 0.")

	if(log) return(log(unname(z))) else return(unname(z))
}



## PROPOSAL MECHANISM ##
#rjmcmc proposal mechanism: proposal mechanism for traversing model complexity (by one parameter at a time)
#author: JM EASTMAN and J UYEDA 2013
.splitormerge <- function(x, delta, control, cache) {
	phy=cache$phy
	fdesc=cache$desc$fdesc
	adesc=cache$desc$adesc
	root=cache$root
	rootd=fdesc[[root]]
	nm=phy$edge[,2]
	
	bb=delta
	vv=x
	
	.shifts.simulation <- function(phy, exclude=NULL)
	{
		drp=phy$edge[,2]%in%exclude
		nm=phy$edge[!drp,2]
		nd=nm[sample(1:length(nm), 1)]
		nd
	}
	
	.splitrate <-
	function(value, n.desc, n.split, lim=list(min=0, max=Inf)){
		if(!.check.lim(value, lim)) stop("Rate appears out of bounds.")
		while(1) {
			u=runif(1, -n.desc*value, n.split*value)
			nr.desc=value+u/n.desc
			nr.split=value-u/n.split
			
			if(.check.lim(c(nr.desc, nr.split), lim)) break()
		}
		return(list(nr.desc=nr.desc, nr.split=nr.split))
	}
	
	.splitvalue <-
	function(cur.vv, n.desc, n.split, factor=log(2)){
		dev=cur.vv-.adjustvalue(cur.vv, factor)
		nr.desc=cur.vv + dev/n.desc
		nr.split=cur.vv - dev/n.split
		return(list(nr.desc=nr.desc, nr.split=nr.split))
	}
	
	
	s=.shifts.simulation(cache$phy, exclude=control$excludeSHIFT)
	marker=match(s, nm)
	if(s%in%rootd){
		tmp=bb
		tmp[marker]=1-tmp[marker]
		mr=match(rootd,nm)
		if(all(tmp[mr]==1)){
			sd=rootd[!rootd==s]
			s=sd[sample(1:length(sd),1)]
			marker=match(s,nm)
			new.bb=bb
			new.bb[marker]=1-new.bb[marker]
		} else {
			new.bb=tmp
		}
	} else {
		new.bb=bb
		new.bb[marker]=1-new.bb[marker]
	}
	new.vv=vv
	cur.vv=vv[marker]
    
	shifts=nm[bb>0]
	K=sum(bb)
	Nk=nrow(phy$edge)-length(control$excludeSHIFT)
	logspace=TRUE
	
	if(sum(new.bb)>sum(bb)) {			## add transition: SPLIT
        if(sum(new.bb)==Nk){
            return(list(x=x, delta=delta, lnpriorproposalRatio=0, decision="none")) ## CANNOT SPLIT
        }
		decision="split"
		n.desc=length(s.desc<-.opensubtree.phylo(s, phy, adesc, shifts))+1
		nd.desc=c(s.desc,s)
		n.split=sum(vv==cur.vv)-n.desc
		if(!logspace) {
			u=.splitvalue(cur.vv=cur.vv, n.desc=n.desc, n.split=n.split, factor=control$prop.width)
		} else {
			u=.splitrate(value=cur.vv, n.desc=n.desc, n.split=n.split, control$rate.lim)
		}
		nr.split=u$nr.split
		nr.desc=u$nr.desc
		new.vv[vv==cur.vv]=nr.split
		ms=match(nd.desc, nm)
		new.vv[ms]=nr.desc
		
        ###
        lnpriorproposalRatio=.lnpriorhastings_ratio.split(K=K, N=Nk, r=cur.vv, r_i=nr.split, r_j=nr.desc, n_i=n.split, n_j=n.desc, fun_k=control$dlnSHIFT, fun_v=control$dlnRATE, delta=1)
	} else {							## drop transition: MERGE
		decision="merge"
		ca.vv=length(which(vv==cur.vv))
		anc = .get.ancestor.of.node(s, phy)
		if(!is.root(anc, phy)) {			# base new rate on ancestral rate of selected branch
			anc.vv=as.numeric(vv[match(anc,nm)])
			na.vv=length(which(vv==anc.vv))
			nr=(anc.vv*na.vv+cur.vv*ca.vv)/(ca.vv+na.vv)
			new.vv[vv==cur.vv | vv==anc.vv]=nr
		} else {							# if ancestor of selected node is root, base new rate on sister node
			sister.tmp=.get.desc.of.node(anc,phy)
			sister=sister.tmp[sister.tmp!=s]
			sis.vv<-anc.vv<-as.numeric(vv[match(sister,nm)])
			ns.vv<-na.vv<-length(which(vv==sis.vv))
			nr=(sis.vv*ns.vv+cur.vv*ca.vv)/(ca.vv+ns.vv)
			new.vv[vv==cur.vv | vv==sis.vv]=nr
		}
		
        ###
        lnpriorproposalRatio=.lnpriorhastings_ratio.merge(K=K, N=Nk, r=nr, r_i=cur.vv, r_j=anc.vv, n_i=ca.vv, n_j=na.vv, fun_k=control$dlnSHIFT, fun_v=control$dlnRATE)
	}
	
	new.values=new.vv
	
	return(list(x=new.vv, delta=new.bb, lnpriorproposalRatio=lnpriorproposalRatio, decision=decision))
}

## PROPOSAL MECHANISM ##
#author: JM EASTMAN and J UYEDA 2013 
.lnpriorhastings_ratio.split=function(K, N, r, r_i, r_j, n_i, n_j, fun_k, fun_v, delta=1){
    #K: number of current shifts
    #N: number of tips in bifurcating tree
    #r: current rate
    #r_i: previous rate for class i
    #r_j: previous rate for class j
    #n_i: number of branches in class i
    #n_j: number of branches in class j
    #fun_k: log-prior function for shifts
    #fun_k: log-prior function for rates
    #delta: tuning parameter

    #  from J UYEDA
    #priors
    lpk1=fun_k(K+1)
    lpk=fun_k(K)
    
    lpri=fun_v(r_i)
    lprj=fun_v(r_j)
    lpr=fun_v(r)
    
    #proposals
    vi=-r_i*n_i
    vj=r_j*n_j
    
    lk1=log(K+1)
    ln2k=log(2*N-2-K)
    lvij=log(vj-vi)
    lnij=log(n_i+n_j)
    ldelta=log(delta)
    
    num=lpk1+lk1+lpri+lprj+lvij+lnij+ldelta
    den=lpk+ln2k+lpr
    
    num-den
}

## PROPOSAL MECHANISM ##
#author: JM EASTMAN and J UYEDA 2013 
.lnpriorhastings_ratio.merge=function(K, N, r, r_i, r_j, n_i, n_j, fun_k, fun_v){
    #K: number of current shifts
    #N: number of tips in bifurcating tree
    #r: merged rate
    #r_i: previous rate for class i
    #r_j: previous rate for class j
    #n_i: number of branches in class i
    #n_j: number of branches in class j
    #fun_k: log-prior function for shifts
    #fun_k: log-prior function for rates
    
    #  from J UYEDA
    #priors
    lpk1=fun_k(K-1)
    lpk=fun_k(K)
    
    lpri=fun_v(r_i)
    lprj=fun_v(r_j)
    lpr=fun_v(r)
    
    #proposals
    lninj=log(n_i*n_j)
    lninj2=2*(log(n_i+n_j))

    lk=log(K)
    ln2k=log(2*N-1-K)
    
    num=lpk1+ln2k+lpr+lninj
    den=lpk+lk+lpri+lprj+lninj2

    num-den
}



## PROPOSAL MECHANISM ##
#mcmc proposal mechanism: scale a value given a proposal width ('prop.width')
#author: JM EASTMAN 2010
.adjustvalue <-
function(value, prop.width) {
# mod 10.20.2010 JM Eastman
	
	vv=value
	if(runif(1)<0.95 | length(unique(vv))==1) {
		rch <- runif(1, min = -prop.width, max = prop.width)
		return(vv[]+rch)
	} else {
		return((vv-mean(vv))*runif(1,min=-prop.width, max=prop.width)+mean(vv))
	}
}



## PROPOSAL MECHANISM ##
#mcmc proposal mechanism: scale a rate (limited to only positive values) given a proposal width ('prop.width')
#author: JM EASTMAN 2010
.adjustrate <-
function(rate, prop.width) {
# mod 05.06.2011 JM Eastman
	vv=rate
	v=log(vv)
	rch <- runif(1, min = -prop.width, max = prop.width)
	v=exp(v[]+rch)
	v
}



## PROPOSAL MECHANISM ##
#proposal utility: move shift-point in tree for rjmcmc.bm()
#author: JM EASTMAN 2011
# fixed bug (10/2012) concerning shifts round the root
.adjustshift <- function(x, delta, rootv, control, cache){
	
	values=x
    rootv=rootv
	fdesc=cache$desc$fdesc
	phy=cache$phy
	root=cache$root
	
	
	root.des=fdesc[[root]]
	names(delta)<-names(values)<-phy$edge[,2]
	shifts=as.numeric(names(which(delta==1)))
	node=shifts[sample(1:length(shifts),1)]
	val=as.numeric(values[which(names(values)==node)])
	
	# select new node
	tmp=.treeslide(node=node, up=NA, use.edges=FALSE, cache=cache)
	newnode=tmp$node
	direction=tmp$direction
	
	# store all current shifts
	shifts=shifts[shifts!=node]
	
	# determine if new shift is non-permissible
	if(newnode%in%c(node,shifts) | newnode%in%control$excludeSHIFT) {
		return(list(new.delta=delta, new.values=values, lnHastingsRatio=0, direction="none"))
	} else {		
		
		# update delta
		new.delta=delta
		new.delta[match(c(newnode, node),names(delta))]=1-delta[match(c(newnode, node),names(delta))]
		
		# update values
		new.values=values
		if(direction=="up") {
			new.values=.assigndescendants(values, newnode, val, exclude=shifts, cache=cache)
			h = ((1/2)*(1/length(fdesc[[newnode]]))) / (1/2)
			lnh=log(h)
		} else if(direction=="down") {
			anc=.get.ancestor.of.node(node, phy)
			if(anc==root){
				anc.val=rootv
			} else {
				anc.val=as.numeric(values[which(names(values)==anc)])
			}
			new.values=.assigndescendants(values, anc, anc.val, exclude=shifts, cache=cache)
			new.values=.assigndescendants(new.values, newnode, val, exclude=shifts, cache=cache)
			h = (1/2) / ((1/2)*(1/length(fdesc[[node]])))
			lnh=log(h)
		} else if(direction=="root"){
            new.values=.assigndescendants(values, newnode, val, exclude=shifts, cache=cache)
            new.values=.assigndescendants(new.values, node, rootv, exclude=shifts, cache=cache)
			lnh=0
		}
		return(list(new.delta=new.delta, new.values=new.values, lnHastingsRatio=lnh, direction=direction))
	}
}



## PROPOSAL MECHANISM ##
#proposal utility: move jump-point in tree for mcmc.levy
#author: JM EASTMAN 2011
.adjustjump <-
function(jumps, add=FALSE, drop=FALSE, swap=FALSE, control, cache){
	# jump.table: orient as with phy$edge[,2]
	# jumps: list of nodes
	
	lim=control$jump.lim
	if(lim!=1) stop("Cannot yet accommodate more than a single jump per branch.")
	origjumps=jumps
	phy=cache$phy
	fdesc=cache$desc$fdesc
	nd=phy$edge[,2]
	xx=which(jumps>0)
	jnd=nd[xx]
	nnd=c(nd[!(nd%in%c(jnd, control$excludeJUMP))])
	lnh=0
		
	if(add==TRUE){
		
		ajnl=nnd[sample(1:length(nnd), 1)]
#		h=(1/(length(jnd)+1)) / (1/length(nnd))
#		lnh=log(h)
		mm=match(ajnl,nd)
		jumps[mm]=jumps[mm]+1
		
#		jj=.treeslide(node=NULL, use.edges=FALSE, cache=cache)$node
#		if(jj%in%control$excludeJUMP){
#			return(list(jumps=origjumps, lnHastingsRatio=0))
#		}
#		mm=match(jj,nd)
#		jumps[mm]=jumps[mm]+1
#		h=(1/(length(jnd)+1)) / (1/length(nnd))
#		lnh=log(h)
#		if(any(jumps>lim)) {
#			return(list(jumps=origjumps, lnHastingsRatio=0))
#		}
				
	} else if(drop==TRUE){
		
		if(length(jnd)) {
			ajnl=jnd[sample(1:length(jnd), 1)]
			mm=match(ajnl,nd)
			jumps[mm]=jumps[mm]-1
		}
		
		
#		h=(1/(length(nnd)+1)) / (1/length(jnd))
#		lnh=log(h)
				
#		if(!length(jnd)) {
#			return(list(jumps=origjumps, lnHastingsRatio=0))
#		}
		
		# select branch
		
#		if(drop==TRUE) {
#			
#			jumps[mm]=jumps[mm]-1
#			node.from=ajl
#			node.to=NULL
#			
#			h = (1/length(nnd)) / (1/length(jnd))
#			lnh=log(h)
#			todo="drop"
			
	} else if(swap==TRUE){
		if(length(nnd)){
			
			ajl = jnd[sample(1:length(jnd),1)]
			mm = match(ajl, nd)
			jumps[mm]=jumps[mm]-1
			
			ajnl=nnd[sample(1:length(nnd), 1)]
			nn=match(ajnl, nd)
			jumps[nn]=jumps[nn]+1
		}

#		h=(1/length(nnd)) / (1/length(jnd))
#		lnh=log(h)
		
		
		
		
#		ajnl<-tmp$node
#		if(ajnl%in%control$excludeJUMP){
#			return(list(jumps=origjumps, lnHastingsRatio=0))
#		} 
		
#		if(any(jumps>lim)) {
#			return(list(jumps=origjumps, lnHastingsRatio=0))
#		}
		
#		from=ajl
#		to=ajnl
#			print(c(from,to))
#		todo="swap"
		
#		direction=tmp$direction
#		if(from==to) {
#			lnh=0
#		} else {
#			if(direction=="up") {
#				h = ((1/2)*(1/length(fdesc[[to]]))) / (1/2)
#				lnh=log(h)
#			} else if(direction=="down") {
#				h = (1/2) / ((1/2)*(1/length(fdesc[[from]])))
#				lnh=log(h)
#			} else {
#				lnh=0
#			}
#		}

#			zz=sample(c("v","s","r"),1)
#			print(zz)
#			if(zz=="v"){
#				tmp<-.treeslide(node=ajl, up=NA, use.edges=FALSE, cache=cache) # up or down
#			} else if(zz=="s") {
#				anc=.get.ancestor.of.node(ajl, phy)								# sister
#				dd=fdesc[[anc]]
#				dd=dd[dd!=ajl]
#				ss=dd[sample(1:length(dd),1)]
#				tmp<-list(node=ss, direction="over")
#			} else {
#				print("here")
				#			}

			
		
	} else {
		stop("Must 'add', 'drop', or 'swap' jumps.")
	}
	
	return(list(jumps=jumps, lnHastingsRatio=lnh))
}


## PROPOSAL MECHANISM ##
#rjmcmc proposal mechanism: locally select an edge based on subtended or subtending branch lengths
#author: JM EASTMAN 2011
.treeslide <-
function(node=NULL, up=NA, use.edges=FALSE, cache){
	
	phy=cache$phy
	nd=phy$edge[,2]
	root=cache$root
	fdesc=cache$desc$fdesc
	edges=cache$edge.length.cs
	
	# choose random node if none given
	if(is.null(node)){
		if(!use.edges) node=nd[sample(1:length(nd), 1)] else node=nd[min(which(runif(1, max=max(edges))<edges))]
		return(list(node=node, direction="none"))
	} else {
		
		# choose direction if none given
		if(is.na(up)) up=as.logical(round(runif(1)))
		
		if(up){
			choice.tmp=.get.ancestor.of.node(node,phy)
			if(choice.tmp==root){	# choice traverses around root
				sisters.tmp=fdesc[[choice.tmp]]
				sister=sisters.tmp[which(sisters.tmp!=node)]
				if(length(sister)==1) {
					choice=sister
					direction="root"
				} else {
					if(!use.edges) sister.edges=rep(1,length(sister.edges)) else sister.edges=phy$edge.length[match(sister, nd)]
					see=cumsum(sister.edges)
					see=see/max(see)
					choice=sister[min(which(runif(1)<see))]
					direction="root"
				}
			} else {				# choice is rootward
				choice=choice.tmp
				direction="up"
			}
		} else {					# choice is tipward
			descendants=fdesc[[node]]
			if(!length(descendants)) {		# -- choice is a tip
				if(runif(1)<0.5) {
					choice=node
					direction="none"
				} else {
					tmp=.treeslide(node=node, up=TRUE, use.edges=use.edges, cache=cache)
                    choice=tmp$node
					direction=tmp$direction
				}
			} else {						# -- choice is an internal node
				if(!use.edges) desc.edges=rep(1,length(descendants)) else desc.edges=phy$edge.length[match(descendants, nd)]
				dee=cumsum(desc.edges)
				dee=dee/max(dee)
				choice=descendants[min(which(runif(1)<dee))]
				direction="down"
			}
		}
		
		return(list(node=choice, direction=direction))
	}	
}


## PROPOSAL UTILITY ##
.jumps.edgewise=function(phy){
	edge=phy$edge[,2]
	mm=match(1:(Ntip(phy)+Nnode(phy)), edge)
	jj=numeric(length(edge))
	jumps.edgewise=function(jump.table){
		jj.tmp=table(jump.table)
		jj[mm[as.numeric(names(jj.tmp))]]=jj.tmp
		jj
	}
	jumps.edgewise
}



## PROPOSAL UTILITY ##
#general phylogenetic utility: recurse down tree changing values of descendants to 'value' until an 'excluded' descendant subtree is reached
#author: JM EASTMAN 2011
.assigndescendants <-
function(vv, node, value, exclude=NULL, cache){
	phy=cache$phy
	adesc=cache$desc$adesc
	dd=c(node, .opensubtree.phylo(node, phy, adesc, exclude))
	vv[match(dd, phy$edge[,2])]=value
	return(vv)
}




## PROPOSAL UTILITY ##
.opensubtree.phylo=function (node, phy, adesc, exclude = NULL) 
{
	N=as.integer(Ntip(phy))
	n=as.integer(Nnode(phy))
	node=as.integer(node)
	exclude=as.integer(exclude)
	dat=list(N=N, n=n, node=node, exclude=exclude)
	res=.Call("open_subtree", dat=dat, desc=adesc, package="geiger")
	res
}




## PROPOSAL MECHANISM ##
#rjmcmc proposal mechanism for updating lineage-specific relative rates (i.e., a subvector swap)
#author: JM EASTMAN 2011
.tune.rate <-
function(rates, control) {
	prop.width=control$prop.width
	tuner=control$tune.scale
	lim=control$rate.lim
	
	ss=rates[sample(1:length(rates), 1)]
	ww=which(rates==ss)
	
	if(runif(1)<tuner){
		nn=.proposal.slidingwindow(ss, prop.width, lim)
	} else {
		nn=.proposal.multiplier(ss, prop.width, lim)
	}
	
	nv=nn$v
	new.rates=rates
	new.rates[ww]=nv
	lhr=nn$lnHastingsRatio
	
	return(list(values=new.rates, lnHastingsRatio=lhr))

}

## PROPOSAL MECHANISM ##
#rjmcmc proposal mechanism for updating a single numeric class from a vector
#author: JM EASTMAN 2011
.tune.value <-
function(values, control) {
	prop.width=control$prop.width
	tuner=control$tune.scale
	lim=control$root.lim
	
	ss=values[sample(1:length(values), 1)]
	ww=which(values==ss)	
	
	if(runif(1)<tuner){
		nn=.proposal.slidingwindow(ss, prop.width, lim)
	} else {
		nn=.proposal.multiplier(ss, prop.width, lim)
	}
	
	nv=nn$v
	lhr=nn$lnHastingsRatio
	lpr=0
	
	values[ww]=nv
	
	return(list(values=values, lnHastingsRatio=lhr, lnPriorRatio=lpr))
}




## PROPOSAL MECHANISM ##
#rjmcmc proposal mechanism for updating a single numeric class from a vector
#author: JM EASTMAN 2011
.tune.SE <-
function(values, control) {
	prop.width=control$prop.width
	tuner=control$tune.scale
	lim=control$se.lim
	
	ss=values[sample(1:length(values), 1)]
	ww=which(values==ss)	
	
	if(runif(1)<tuner){
		nn=.proposal.slidingwindow(ss, prop.width, lim)
	} else {
		nn=.proposal.multiplier(ss, prop.width, lim)
	}
	
	nv=nn$v
	lhr=nn$lnHastingsRatio
	lpr=0
	
	values[ww]=nv

	return(list(values=values, lnHastingsRatio=lhr, lnPriorRatio=lpr))
}



## PROPOSAL MECHANISM ##
#rjmcmc proposal mechanism: adjust a value (within bounds)
#author: JM EASTMAN 2011
.proposal.slidingwindow <- function(value, prop.width, lim=list(min=-Inf, max=Inf)){
	if(!.check.lim(value, lim)) stop("Values appear out of bounds.")
	min=lim$min
	max=lim$max
	
	while(1){
		u=runif(1)
		v=value+(u-0.5)*prop.width
		
		# reflect if out-of-bounds
		if(any(v>max)) {
			v[v>max]=max-(v[v>max]-max)
		}
		if(any(v<min)){
			v[v<min]=min-(v[v<min]-min)
		}
		if(.check.lim(v, lim)) break()
	}
		
	return(list(v=v, lnHastingsRatio=0))
}


## PROPOSAL MECHANISM ##
#rjmcmc proposal mechanism: scale a value with asymmetrically drawn multiplier
#author: JM EASTMAN 2011
# from Lakner et al. 2008 Syst Biol
.proposal.multiplier <- function(value, prop.width, lim=list(min=-Inf, max=Inf)){
	if(!.check.lim(value, lim)) stop("Values appear out of bounds.")
	
#	tmp=c(prop.width, 1/prop.width)
#	a=min(tmp)
#	b=max(tmp)
#	lambda=2*log(b)
	while(1){
		m=exp(prop.width*(runif(1)-0.5))
#		m=exp(lambda*(u-0.5))
		v=value*m
		if(.check.lim(v, lim)) break()
	}
	return(list(v=v, lnHastingsRatio=log(m)))
}



