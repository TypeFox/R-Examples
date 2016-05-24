.get.simulation.matrix=function(phy){
	N=Ntip(phy)
	n=nrow(phy$edge)
	dd=.cache.descendants(phy)$tips
	m=matrix(0, N, n)
	edg=phy$edge.length
	idx=phy$edge[,2]
	for(x in 1:length(edg)){
		edge=edg[x]
		m[dd[[idx[x]]],x]=sqrt(edge)
	}
	m
}

.check.Qmatrix=function(Q){
    m=unique(dim(Q))
    if(length(m)>1) stop("'Q' must be a square matrix")
    didx=1 + 0L:(m - 1L) * (m + 1)
    if(!all(abs(rowSums(Q))<0.000001)) stop("rows of 'Q' must sum to zero")
    if(!all(Q[didx]<=0)) stop("diagonal elements of 'Q' should be negative")
    if(!all(Q[-didx]>=0)) stop("off-diagonal elements of 'Q' should be positive")
}

.make.modelmatrix=function(m, model=c("BM", "speciational", "discrete")){
	model=match.arg(model, c("BM", "speciational", "discrete"))
	if(model=="discrete"){
		if(is.matrix(m)){
			m=list(m)
            for(j in 1:length(m)){
                .check.Qmatrix(m[[j]])
            }
		}
	} else {
		if(is.numeric(m)) m=as.matrix(m) else stop("Supply 'm' as a matrix of rates")
        if(any(diag(m)<0)) stop("'m' appears to have negative variance component(s)")
	}
	return(m)
}

sim.char <- function(phy, par, nsim=1, model=c("BM", "speciational", "discrete"), root=1)
{
	model=match.arg(model, c("BM", "speciational", "discrete"))
	
	model.matrix=.make.modelmatrix(par, model)
	
	nbranches<-nrow(phy$edge)
	nspecies<-Ntip(phy)
    
    if(length(root)>1) stop("'root' should be a single value")
	
	if(model%in%c("BM", "speciational"))
	{
		m<-.get.simulation.matrix(phy)
		
		if(model=="speciational") {
			m[m>0]<-1.0;
		}
		
		nchar<-nrow(model.matrix)
		rnd<-t(mvrnorm(nsim*nbranches, mu=rep(0, nchar), Sigma=model.matrix))
        rnd<-array(rnd, dim=c(nchar, nbranches, nsim));
		
		simulate<-function(v, root) (m %*% as.matrix(v))+root;
		
		result<-apply(rnd, 1, simulate, root)
		result<-aperm(array(result, dim=c(nspecies, nsim, nchar)), c(1, 3, 2))
		
		rownames(result)<-phy$tip.label;
	} else {
        rt=nspecies+1
        zphy=reorder.phylo(phy, "postorder")
        el=zphy$edge.length
		nchar<-length(model.matrix);
		result<-array(0, dim=c(nspecies, nchar, nsim))
        .get.state=function(s, p){
            pp=cumsum(p[s,])
            min(which(runif(1)<pp))
        }
		for(j in 1:nchar) {
            m=model.matrix[[j]]
            if(!root%in%c(1:nrow(m))) stop(paste("'root' must be a character state from 1 to ", nrow(m), sep=""))
            p=lapply(el, function(l) matexpo(m*l))
            
			for(k in 1:nsim) {
                node.value<-numeric(nspecies+Nnode(zphy))
                node.value[rt]<-root
	   			for(i in nbranches:1) {
                    cur=zphy$edge[i,2]
                    anc=zphy$edge[i,1]
                    curp=p[[i]]
                    s=node.value[anc]
					node.value[cur]=.get.state(s, curp)
	   			}
	   			result[,j,k]<-node.value[1:nspecies]
	   		}
		}
		rownames(result)<-zphy$tip.label;
		
	}
	return(result);
}


.check.stoppingcrit=function(time, taxa){
	flag=FALSE
	if(is.na(time)) time=0
	if(is.na(taxa)) taxa=0
	if(!( (tm<-is.numeric(time)) | (tx<-is.numeric(taxa)))) flag=TRUE 
	if(!flag) if(time==0 & taxa==0) flag=TRUE
	if(!flag) if(time!=0 & taxa!=0) flag=TRUE
	if(flag) stop("Either 'time.stop' or 'taxa.stop' must be specified as a stopping criterion")
	
	ww=which(c(tm, tx))
	xx=which(c(time, taxa)[ww]!=0)
	return(c("time", "taxa")[ww][xx])
}



.set.seed.clock <- function(print=FALSE){
	date = date()
 	seed1 = as.numeric(strsplit(substring(date,12,19),":")[[1]])%*%c(1,100,10000)
 	seed <- runif(1, min=0, max=50) * seed1
 	set.seed(seed)
 	if(print) cat("seed = ", seed, "\n");
 	seed[1,1]
}


sim.bdtree <- function (b=1, d=0, stop=c("taxa", "time"), n=100, t=4, seed=0, extinct=TRUE) {
# December 6 2005 Jason T. Weir
# Modified by Luke J. Harmon
# Modified by JM Eastman
# The following simulates Yule trees to a given time T
	
	stop = match.arg(stop, c("taxa", "time"));
	time.stop <- taxa.stop <- 0;
	if (stop == "taxa") taxa.stop = n + 1;
	if (stop == "time") time.stop = t;
	if (time.stop == 0 & taxa.stop == 0) stop("Stopping criterion ('n' or 't') must be provided");
	if (seed == 0) {
		seed = .set.seed.clock(print = FALSE);
	} else {
		set.seed(seed); # this condition was missing JWB
	}
	return.all.extinct = extinct;
	
	while (1) {
		edge <- rbind(c(1, 2), c(1, 3)); # this is a starting edge matrix
		edge.length <- rep(NA, 2);
		stem.depth <- numeric(2);
		alive <- rep(TRUE, 2); # marker for live lineages
		t <- 0; # time at any point in the tree
		next.node <- 4;
		
############
		repeat {
			if (taxa.stop) {
				if (sum(alive) >= taxa.stop) break;
			}
			if (sum(alive) == 0) break;
			dt <- rexp(1, sum(alive) * (b + d));
			t <- t + dt;
			if (time.stop) {
				if (t >= time.stop) {
					t <- time.stop;
					break;
				}
			}
			r <- runif(1);
			if (r <= b/(b + d)) { ###4 #this creates a bifucation in the tree
				random_lineage <- round(runif(1, min = 1, max = sum(alive)));
				e <- matrix(edge[alive,], ncol = 2);
				parent <- e[random_lineage,2];
				alive[alive][random_lineage] <- FALSE;
				edge <- rbind(edge, c(parent, next.node), c(parent, next.node + 1));
				next.node <- next.node + 2;
				alive <- c(alive, TRUE, TRUE);
				stem.depth <- c(stem.depth, t, t);
				x <- which(edge[,2] == parent);
				edge.length[x] <- t - stem.depth[x];
				edge.length<-c(edge.length, NA, NA)
			}###4
			
			else {###4 This terminates one of the current lineages on the tree
				random_lineage <- round(runif(1, min = 1, max = sum(alive)));
				edge.length[alive][random_lineage] <- t - stem.depth[alive][random_lineage];
				alive[alive][random_lineage] <- FALSE;
			}###4
		}#1A
		
		if (return.all.extinct == TRUE | sum(alive) > 1) break;
	}
	edge.length[alive] <- t - stem.depth[alive];
	n <- -1;
	for (i in 1:max(edge)) {
		if (any(edge[,1] == i)) {
			edge[which(edge[,1] == i), 1] <- n;
			edge[which(edge[,2] == i), 2] <- n;
			n <- n - 1;
		}
	}
	
	edge[edge > 0] <- 1:sum(edge > 0);
	tip.label <- 1:sum(edge > 0);
	mode(edge) <- "character";
	mode(tip.label) <- "character";
	obj <- list(edge = edge, edge.length = edge.length, tip.label=tip.label);
	class(obj) <- "phylo";
	obj <- old2new.phylo(obj);
	obj <- read.tree(text = write.tree(obj));
	attr(obj, "seed") = seed;
	if (stop == "taxa") {
		drp = obj$edge[min(which(obj$edge.length == 0)),2];
		obj = .drop.tip(obj, obj$tip.label[drp]);
	}
	obj$tip.label = paste("s", 1:Ntip(obj), sep = "");
	return (obj);
}






######################################################################################
##STOCHASTIC SIMULATION OF A TIME-HOMOGENOUS BIRTH-DEATH PROCESS
##PARAMETERS: N0=starting number of lineages, ks=speciation rate, ep=relative rate of extinction, finaltime=endpoint for simulation
######################################################################################
sim.bd <- function (b=1, d=0, n0=1, times=0:4, seed=0) {
	
	if (seed == 0) {
		seed = .set.seed.clock(print = FALSE);
	} else {
		set.seed(seed);
	}
	
	n <- n0;
	t <- 0;

	times <- unique(c(0, sort(times)));
	pop <- numeric(length(times));
	pop[] <- 0;
	pop[1] <- n;

	i <- 2;
	while (1) {
		if (t > times[i]) {
			m = max(which(t > times));
			pop[i:m] = n;
			i = m + 1;
		}

		waittime <- rexp(1, rate = n * (b + d));
		t <- t + waittime;
		if (t > max(times)) {
			break;
		} else {
			ran <- runif(1);
			if (ran <= b/(b + d)) {
				n <- n + 1;
			} else {
				n <- n - 1;
			}
		}
		if (n == 0) break;
	}
	pop[i:length(times)] = n;
	res = cbind(times, pop);
	colnames(res) = c("time", "n");
	return (res);
}



