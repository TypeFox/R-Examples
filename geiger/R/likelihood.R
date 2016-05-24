## replacing prepare.data.bm
make.bm.relaxed <- function(phy, dat, SE=NA, method=c("direct","vcv","reml")){
#   method=match.arg(method, c("direct","vcv","reml"))
    method=match.arg(method, c("direct","vcv"))

	lik=switch(method,
			   direct=.make.bm.relaxed.direct(phy, dat, SE),
			   vcv=.make.bm.relaxed.vcv(phy, dat, SE),
			   reml=.make.bm.relaxed.pic(phy, dat, SE))
	lik
}


#primary function for computing the likelihood of data, given a root state, VCV matrix, and Brownian motion model
#author: LJ HARMON 2009 and JM EASTMAN 2010
.bm.lik.fn.vcv <-
function(root, dat, vcv, SE) { 
# mod 12.02.2010 JM Eastman: using determinant()$modulus rather than det() to stay in log-space
	y=dat
	b <- vcv
	if(any(SE!=0)) diag(b)=diag(b)+SE^2
	w <- rep(root, nrow(b))
	num <- -t(y-w)%*%solve(b)%*%(y-w)/2
	den <- 0.5*(length(y)*log(2*pi) + as.numeric(determinant(b)$modulus))
	return((num-den)[1,1])
}


#primary function for computing the likelihood of data, using REML, under Brownian motion model 
#author: LJ HARMON, LJ REVELL, and JM EASTMAN 2011
#'rphy' is rate-scaled tree in pruning-wise order (as 'ic')
.bm.lik.fn.reml <- 
function(rphy, ic) {
    new.var=.pic_variance.phylo(rphy)
	reml=dnorm(ic, mean=0, sd=sqrt(new.var), log=TRUE)
	return(sum(reml))
}


#compute expected PIC variance given tree: used for .bm.lik.fn.reml()  
#author: JM EASTMAN 2011
.pic_variance.phylo <- function(phy)
{
	phy$node.label=NULL
	nb.tip <- Ntip(phy)
    nb.node <- phy$Nnode
	if(!is.binary.tree(phy)) stop("'phy' is not fully dichotomous.")
	
    phy <- reorder(phy, "postorder")
	
	ans <- .C("pic_variance", as.integer(nb.tip), as.integer(nb.node),
              as.integer(phy$edge[, 1]), as.integer(phy$edge[, 2]),
              as.double(phy$edge.length), double(nb.node),
              PACKAGE = "geiger")
	
    var <- ans[[6]]
    var
}

.make.bm.relaxed.vcv <- function(phy, dat, SE=NULL){
    ## NOT currently handling given node states
    
	cache=.prepare.bm.univariate(phy, dat, nodes=NULL, SE=SE, control=list(binary=FALSE))
	
	check.argn=function(rates, root){
		if(length(rates)!=(cache$n.node+cache$n.tip-1) && length(root)==1) stop("Supply 'rates' as a vector of rate scalars for each branch, and supply 'root' as a single value.")
	}
	
    adjvar = as.integer(attributes(cache$y)$adjse)
 
    if(any(adjvar==1)){ # adjustable SE
		likSE <- function(rates, root, SE) {
			check.argn(rates, root)
            tt=.scale.brlen(cache$phy, rates)
            vcv=.paths.phylo(tt, vcv=TRUE)
            xSE=cache$SE[mm<-match(colnames(vcv), names(cache$SE))]
            adjvar=adjvar[mm]
            xSE[which(adjvar==1)]=SE
            dat=cache$dat[match(colnames(vcv), names(cache$dat))]
            .bm.lik.fn.vcv(root, dat, vcv, xSE)
        }
		attr(likSE, "cache") <- cache
		attr(likSE,"argn")=list(rates=cache$phy$edge[,2], root="root", SE="SE")
		class(likSE)=c("rbm","bm","function")
		return(likSE)
	} else { # unadjustable SE
		lik <- function(rates, root){
            check.argn(rates, root)
            tt=.scale.brlen(cache$phy, rates)
            vcv=.paths.phylo(tt, vcv=TRUE)
            SE=cache$SE[match(colnames(vcv),names(cache$SE))]
            dat=cache$dat[match(colnames(vcv),names(cache$dat))]
            .bm.lik.fn.vcv(root, dat, vcv, SE)
        }
        attr(lik, "cache") <- cache
        attr(lik,"argn")=list(rates=cache$phy$edge[,2], root="root")
        class(lik)=c("rbm","bm","function")
        return(lik)
    }
}

.make.bm.relaxed.pic <- function(phy, dat, SE=NULL){
    
    ## NOT currently handling given node states
	cache=.prepare.bm.univariate(phy, dat, nodes=NULL, SE=SE)
	
	ic=pic(cache$dat, cache$phy, scaled=FALSE)
	cache$ic=ic
	
	check.argn=function(rates){
		if(length(rates)!=(cache$n.node+cache$n.tip-1)) stop("Supply 'rates' as a vector of rate scalars for each branch.")
	}
	
    adjvar = as.integer(attributes(cache$y)$adjse)

    xSE=cache$SE
    N=Ntip(cache$phy)
    subSE=match(1:N, cache$phy$edge[,2])

    
    if(any(adjvar==1)){ # adjustable SE
		likSE <- function(rates, SE) {
			check.argn(rates)
            xSE[adjvar==1]=SE
            rphy=.scale.brlen(cache$phy, rates)
            rphy$edge.length[subSE]=rphy$edge.length[subSE]+xSE^2
            .bm.lik.fn.reml(rphy, cache$ic)
        }
		attr(likSE, "cache") <- cache
		attr(likSE,"argn")=list(rates=cache$phy$edge[,2], SE="SE")
		class(likSE)=c("rbm","bm","function")
		return(likSE)
	} else {
        lik <- function(rates){
            check.argn(rates)
            rphy=.scale.brlen(cache$phy, rates)
            rphy$edge.length[subSE]=rphy$edge.length[subSE]+xSE^2
            .bm.lik.fn.reml(rphy, cache$ic)
        }
        attr(lik,"cache")=cache
        attr(lik,"argn")=list(rates=cache$phy$edge[,2])
        class(lik)=c("rbm","bm","function")
    }

	lik
}


.make.bm.relaxed.direct <- function (phy, dat, SE=NULL) 
{

    ## NOT currently handling given node states
    
	cache=.prepare.bm.univariate(phy, dat, nodes=NULL, SE=SE)
    N = cache$n.tip
    n = cache$n.node
	z = length(cache$len)
    rootidx = as.integer(cache$root)
    adjvar = as.integer(attributes(cache$y)$adjse)
	given = as.integer(attributes(cache$y)$given)
    given[rootidx]=1
    
    datc_init = list(
                len = as.numeric(cache$len),
				intorder = as.integer(cache$order[-length(cache$order)]), 
				tiporder = as.integer(1:N), 
				root = rootidx, 
                y = as.numeric(cache$y[1, ]),
				n = as.integer(z),
                given = as.integer(given),
				descRight = as.integer(cache$children[, 1]),
				descLeft = as.integer(cache$children[, 2]),
                drift=0
    )
	
	ll.bm.direct <- function(pars, datc) {
        out = .Call("bm_direct", dat = datc, pars = pars, package = "geiger")
#       vals = c(out$initM[rootidx], out$initV[rootidx], out$lq[rootidx])
        loglik <- sum(out$lq)
#       intermediates=FALSE
#       if (intermediates) {
#          attr(loglik, "intermediates") <- intermediates
#          attr(loglik, "vals") <- vals
#       }
        return(loglik)
    }
    class(ll.bm.direct) <- c("bm.direct", "bm", "function")
	
	vv = numeric(N+n)
    mm = match(cache$phy$edge[, 2], 1:(N+n))
	check.argn=function(rates, root){
		if(length(rates)!=(N+n-1) && length(root)==1) stop("Supply 'rates' as a vector of rate scalars for each branch, and supply 'root' as a single value.")
	}
	   
    var = as.numeric(cache$y[2, ]^2)
    
	## LIKELIHOOD FUNCTION
	if(any(adjvar==1)){ # adjustable SE
		likSE <- function(rates, root, SE) {
			check.argn(rates, root)
			vv[mm] = rates
			datc_se=datc_init
			var[which(adjvar==1)]=SE^2
            datc_se$var=var
            
            datc_se$y[rootidx]=root
            
			ll = ll.bm.direct(pars = vv,  datc_se)
			return(ll)
		}
		attr(likSE, "cache") <- cache
		attr(likSE,"argn")=list(rates=cache$phy$edge[,2], root="root", SE="SE")
		class(likSE)=c("rbm","bm","function")
		return(likSE)
	} else { # unadjustable SE
		lik <- function(rates, root) {
			check.argn(rates, root)
			vv[mm] = rates
            datc_init$var=var
            datc_init$y[rootidx]=root

			ll = ll.bm.direct(pars = vv,  datc_init)
			return(ll)
		}
		attr(lik, "cache") <- cache
		attr(lik,"argn")=list(rates=cache$phy$edge[,2], root="root")
		class(lik)=c("rbm","bm","function")
		return(lik)
	}	
}

.cache.tree <- function (phy) 
{
	ordxx=function (children, is.tip, root) 
# from diversitree:::get.ordering
	{
		todo <- list(root)
		i <- root
		repeat {
			kids <- children[i, ]
			i <- kids[!is.tip[kids]]
			if (length(i) > 0) 
            todo <- c(todo, list(i))
			else break
		}
		as.vector(unlist(rev(todo)))
	}
	
    edge <- phy$edge
    edge.length <- phy$edge.length
    idx <- seq_len(max(edge))
    n.tip <- Ntip(phy)
    tips <- seq_len(n.tip)
    root <- n.tip + 1
    is.tip <- idx <= n.tip
	
	desc=.cache.descendants(phy)
    children <- desc$fdesc
	if(!max(sapply(children, length) == 2)){
		children=NULL
		order=NULL
		binary=FALSE
	} else {
		children <- rbind(matrix(NA, n.tip, 2), t(matrix(unlist(children), nrow=2)))
		order <- ordxx(children, is.tip, root)	
		binary=TRUE
	}
	
    len <- edge.length[mm<-match(idx, edge[, 2])]
	
	ans <- list(tip.label = phy$tip.label, node.label = phy$node.label,
				len = len, children = children, order = order, 
				root = root, n.tip = n.tip, n.node = phy$Nnode, tips = tips, 
				edge = edge, edge.length = edge.length, nodes = phy$edge[,2], binary = binary, desc = desc)
    ans
}


.prepare.bm.univariate=function(phy, dat, nodes=NULL, SE=NA, control=list(binary=TRUE, ultrametric=FALSE)){

    ## CONTROL OBJECT
    ct=list(binary=TRUE, ultrametric=FALSE)
    ct[names(control)]=control
    
    ## MATCHING and major problems
    td=treedata(phy, dat, sort=TRUE, warnings=FALSE)
    phy=reorder(td$phy, "postorder")
    if(ct$binary) if(!is.binary.tree(phy)) stop("'phy' should be a binary tree")
    if(ct$ultrametric) if(!is.ultrametric(phy)) stop("'phy' should be an ultrametric tree")
    if(is.null(phy$edge.length)) stop("'phy' must include branch lengths in units of time")

    if(ncol(td$data)>1) stop("'dat' should be univariate")
    dat=td$data[,1]

	## RESOLVE SE
    seTMP=structure(rep(NA, length(dat)), names=names(dat))
    
	if(is.null(SE)) SE=NA
    
    if(length(SE)>1){
        if(is.null(names(SE))) stop("'SE' should be a named vector")
        if(!all(names(dat)%in%names(SE))) stop("names in 'SE' must all occur in names of 'dat'")
        seTMP[names(SE[names(dat)])]=SE[names(dat)]
        SE=seTMP
    } else {
        if(is.numeric(SE)){
            seTMP[]=SE
            SE=seTMP
        } else {
            SE=seTMP
        }
    }
    
    if(!all(is.na(SE) | SE >= 0)) stop("'SE' values should be positive (including 0) or NA")
  
    ## CACHE tree
    cache=.cache.tree(phy)
    N=cache$n.tip
    n=cache$n.node
    m<-s<-g<-numeric(N+n)
    
    ## RESOLVE data: given trait values (m and g) and SE (s) for every node (tips and internals)
    g[1:N]=1
    m[]=NA; m[1:N]=dat
    s[1:N]=SE

    ## RESOLVE nodes
    if(!is.null(nodes)){
        nn=(N+1):(N+n)
        vec=.cache.y.nodes(m, s, g, nn, phy, nodes=nodes)
    } else {
        vec=rbind(m=m, s=s)
        attr(vec, "given")=g
        attr(vec, "adjse")=as.numeric(is.na(s))[1:N]
    }
     	
	cache$SE=SE
	cache$dat=dat[match(phy$tip.label, names(dat))]
	cache$phy=phy
    
    cache$y=vec

    return(cache)
}

.cache.y.nodes=function(m, s, g, nn, phy, nodes){
    
    if(is.numeric(nodes) & is.vector(nodes)){
        if(!all(names(nodes)%in%nn)) stop("'nodes' must have (integer) names corresponding to the internal nodes of 'phy'")
        nodes=data.frame(cbind(node=as.integer(names(nodes)), mean=nodes, SE=0), stringsAsFactors=FALSE)
    } else {
        if(!all(c("taxon1", "taxon2", "mean", "SE")%in%colnames(nodes))){
            flag=FALSE
            if(!all(c("mean", "SE")%in%colnames(nodes)) | is.null(rownames(nodes))){
                flag=TRUE
            } else if(!all(rr<-as.integer(rownames(nodes))%in%nn)){
                flag=TRUE
            }
            if(flag) stop("'nodes' must minimally have column names: 'taxon1', 'taxon2', 'mean', and 'SE'")
            nodes=as.data.frame(nodes)
            nodes$node=as.integer(rownames(nodes))
        } else {
            nodes=as.data.frame(nodes)
            if(!is.numeric(nodes$mean) | !is.numeric(nodes$SE)){
                stop("'nodes' must have numeric vectors for 'mean' and 'SE'")
            }
            
            if(!all(zz<-unique(c(as.character(nodes$taxon1), as.character(nodes$taxon2)))%in%phy$tip.label)){
                stop(paste("Some taxa appear missing from 'phy':\n\t", paste(zz[!zz%in%phy$tip.label], collapse="\n\t", sep=""), sep=""))
            }
            
            nodes$node=apply(nodes[,c("taxon1", "taxon2")], 1, .mrca, phy=phy)
        }
        
        if(!length(unique(nodes$node))==nrow(nodes)) {
            stop("Some nodes multiply constrained:\n\t", paste(nodes$node[duplicated(nodes$node)], collapse="\n\t", sep=""), sep="")
        }
    }
    
    nidx=nodes$node
    if(any(g[nidx]==1->zz)) stop("Some nodes already constrained:\n\t", paste(nidx[which(zz)], collapse="\n\t", sep=""), sep="")
    
    m[nidx]=as.numeric(nodes$mean)
    s[nidx]=as.numeric(nodes$SE)
    g[nidx]=1
    
 	vec=rbind(m=m, s=s)
   	attr(vec, "given")=g
   	attr(vec, "adjse")=as.numeric(is.na(s))
    
	vec
}

.proc.lnR=function(gen, subprop, cur.lnL, new.lnL, lnp, lnh, heat=1, control){
	if(is.infinite(cur.lnL)) stop("Starting point has exceptionally poor likelihood.")
	if(is.finite(new.lnL)) {
		lnLikelihoodRatio = new.lnL - cur.lnL
	} else {
		new.lnL=-Inf
		lnLikelihoodRatio = -Inf
	}

	if(control$sample.priors) lnLikelihoodRatio=0

	lnR=(heat * lnLikelihoodRatio) + (heat * lnp) + lnh

	r=.assess.lnR(lnR)	
	
	if(r$error) .error.rjmcmc(gen, subprop, cur.lnL, new.lnL, lnLikelihoodRatio, lnp, lnh, control$errorlog)
	return(r)
}


#general phylogenetic utility for determining whether to accept ('r'=1) or reject ('r'=0) a proposed model in Markov chain Monte Carlo sampling
#author: JM EASTMAN 2010
.assess.lnR <- function(lnR) {
	if(is.na(lnR) | is.infinite(abs(lnR))) {
		error=TRUE
		r=0
	} else {
		if(lnR < -20) {
			r=0 
		} else if(lnR >= 0) {
			r=1 
		} else {
			r=exp(lnR)
		}
		error=FALSE
	}
	return(list(r=r, error=error))
}



