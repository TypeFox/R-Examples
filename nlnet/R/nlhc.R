# source this code to conduct the Noninear Hierarchical Clustering based on general dependency hiararchical cLustering
# array is the data matrix with no missing values
# hamil.method is the method to find the hamiltonian path. It is passed onto the function tsp of library TSP
#				to use linkern method, the user needs to install concord as instructed in TSP
# use.normal.approx: whether to use the normal approximation to for the null hypothesis. If TRUE, normal approximation is used for every feature, AND all covariances are assumed to be zero. If FALSE, generates permutation based null distribution - mean vector and a variance-covariance matrix.
# normalization: the normalization method for the array. There are three choices - "standardize" means removing the mean of each row and make the standard deviation one; "normal_score" means normal score transformation; "none" means do nothing. In that case we still assume some normalization has been done by the user such that each row has approximately mean 0 and sd 1. 
# combine.linear: whether linear association should be found by correlation to combine with nonlinear association found by DCOL. The two pieces of information is combined at the start to initiate the distance matrix. 
# use.traditional.hclust: whether traditional agglomerative clustering should be used. 
# method.traditional.hclust: the method to pass on to hclust() if traditional method is chosen.
#
#
# Returns a hclust object same as the output of hclust().

################################################################################

nlhc<-function(array, hamil.method="nn", concorde.path=NA, use.normal.approx=FALSE, normalization="standardize", combine.linear=TRUE, use.traditional.hclust=FALSE, method.traditional.hclust="average")
{
    #library(TSP)
	
	if(hamil.method=="linkern")
	{
		if(is.na(concorde.path))
		{
			message("linkern is chosen but concorde path wasn't given.")
			return(-1)
		}else{
			concorde_path(concorde.path)
		}
	}
	
	gene.specific.null<-function(array, B=500)
	{
		null.mat<-matrix(0, nrow=nrow(array), ncol=B)
		l<-ncol(array)
		d.array<-array[,1:(l-1)]
		for(i in 1:B)
		{
			this.order<-sample(l, l, replace=FALSE)
			for(j in 1:(l-1)) d.array[,j]<-abs(array[,this.order[j+1]]-array[,this.order[j]])
			null.mat[,i]<-apply(d.array, 1, sum)
		}
		r<-cbind(apply(null.mat, 1, mean), apply(null.mat, 1, sd))
		return(r)
	}
	
	gene.specific.p<-function(null.distr, new.d)
	{
		for(i in 1:length(new.d))
		{
			new.d[i]<-pnorm(new.d[i], mean=null.distr[i,1], sd=null.distr[i,2], lower.tail=TRUE)
		}
		return(new.d)
	}
	
	sim.emp.distri<-function(m, normalization="standardize", n.perm=1e5)
	{
		r<-rep(0, n.perm)
		if(normalization == "normal_score")
		{
			a<-qnorm((1:m)/(1+m))
			for(i in 1:n.perm)
			{
				a2<-sample(a, m, replace=F)
				r[i]<-sum(abs(diff(a2)))
			}
		}else{
			for(i in 1:n.perm)
			{
				a2<-rnorm(m)
				r[i]<-sum(abs(diff(a2)))
			}
			
		}
		return(r)
	}
	
	scol.matrix.order<-function(a,x) # x is the vector, a is the matrix, find ordered distance of rows.of.a|x
	{
		if(is.null(nrow(a)) | nrow(a) == 1)
		{
			a<-as.vector(a)
			a<-a[order(x)]
			d<-a[2:length(a)]-a[1:(length(a)-1)]
			dd<-sum(abs(d),na.rm=T)
		}else{
			a<-a[,order(x)]
			d<-a[,2:ncol(a)]-a[,1:(ncol(a)-1)]
			dd<-apply(abs(d),1,sum,na.rm=T)
		}
		return(dd)
	}
	
	scol.matrix<-function(a, direction=2)	# when direction is 1, scol.matrix[i,j] = SCOL(a[i,], a[j,]), j|i
	{		
		
		rdmat<-matrix(0, ncol=nrow(a), nrow=nrow(a))
		for(j in 1:nrow(a))
		{
			rdmat[j,]<-scol.matrix.order(a, a[j,])
		}
		
		if(direction == 2)
		{
			rdmat.diff<-rdmat-t(rdmat)
			sel<-which(rdmat.diff > 0)
			rdmat[sel]<-t(rdmat)[sel]
		}
		return(rdmat)
	}
	
	hamil<-function(z, method.1="linkern")  # z is a matrix, with column being genes
	{
		#library(TSP)
		d<-dist(z, method="man")
		tsp<-TSP(d)
		tsp<-insert_dummy(tsp, label = "cut")
		tour<-solve_TSP(tsp, method=method.1)
		path <- cut_tour(tour, "cut")
		
		return(list(path, attributes(tour)$tour_length))
	}
	
	normscore.row<-function(a)
	{
		#library(coin)
		b<-t(apply(a, 1, normal_trafo))
		return(b)
	}
	
	normrow<-function(a)
	{
		m<-apply(a,1,mean,na.rm=T)
		s<-apply(a,1,sd,na.rm=T)
		a<-(a-m)/s
		return(a)
	}
	
	find.min.pos<-function(d)
	{
		pos<-which(d==min(d))[1]
		pos.x<-pos %% nrow(d)
		if(pos.x == 0) pos.x<-nrow(d)
		pos.y<-floor((pos-1)/nrow(d)+1)
		pos<-c(pos.x, pos.y)
		return(pos)
	}
	
	reorg.tree<-function(h)
	{
		h2<-h
		used<-rep(F, nrow(h$merge))
		used[1]<-T
		corresponder<-h$height*0
		corresponder[1]<-1
		
		for(i in 2:nrow(h$merge))
		{
			next.height<-min(h$height[!used])
			sel<-which(h$height==next.height & !used)[1]
			corresponder[i]<-sel
			this<-h$merge[sel,]
			if(this[1] > 0) this[1]<-which(corresponder == this[1]) 
			if(this[2] > 0) this[2]<-which(corresponder == this[2]) 
			
			used[sel]<-T
			h2$merge[i,]<-this
			sel.2<-which(!used)
			h2$height[i]<-h$height[sel]
		}
		return(h2)
	}
	
# initiation
	n<-nrow(array)
	m<-ncol(array)
	
	if(use.normal.approx & normalization=="none") normalization<-"standardize"
	
	if(normalization=="normal_score")
	{ 
		a<-normscore.row(array)
	}else{
		if(normalization=="standardize")
		{ 
			a<-normrow(array) 
		}else{ 
			if(normalization=="none") 
			{
				a<-array
			}else{ 
				return("wrong normalization method.") 
			}
		}
	}
	
	if(!use.normal.approx)
	{
		null.distr<-gene.specific.null(a)
	}else{
		m.emp<-2/sqrt(pi)*(m-1)
		s.emp<-sqrt(m-1)*4*(2*pi+3*sqrt(3)-9)/3/pi
		null.distr<-cbind(rep(m.emp, nrow(a)), rep(s.emp, nrow(a)))
	}
	
	orig.a<-a		# a records cluster profiles; orig.a records gene profiles
	
	sim.mat<-scol.matrix(a,direction=1)  ## similarity matrix by SCOL, asymmetric, column given row
	gc()
	d.mat<-sim.mat
	if(!use.normal.approx)
	{
		for(i in 1:nrow(d.mat)) d.mat[i,]<-log10(gene.specific.p(null.distr, sim.mat[i,]))
	}else{
		d.mat<-log10(pnorm(d.mat, mean=null.distr[1,1], sd=null.distr[1,2], lower.tail=TRUE))
	}
	diag(d.mat)<-0		# d.mat (column given row) records log10 p-value between clusters
	
	if(combine.linear)
	{
		d.linear<-cor(t(a))
		d.linear<- -abs(d.linear)*sqrt((m-2)/(1-d.linear^2))
		d.linear<-2*pt(d.linear, df=m-2, lower.tail=T)
		d.linear<-log10(d.linear)
		diag(d.linear)<-0
		min.non.inf<-min(d.linear[d.linear != -Inf])
		d.linear[d.linear == -Inf]<-min.non.inf
		
		sel<-which(d.linear < d.mat)
		d.mat[sel]<-d.linear[sel]
		rm(d.linear)
	}
	gc()	
	
	if(!use.traditional.hclust)
	{
		grps<-lab<--(1:nrow(a)) 
# grps records which cluster the gene belongs to
# lab points to which row of a has a cluster profile
		
		ptr<-1
		
		h.merge<-matrix(0,ncol=2,nrow=n-1)
		h.height<-rep(0,n-1)
		h.order<-NULL
		
# iteration
# every round, select the lowest distance and merge
		h<-hclust(as.dist(d.mat[1:2,1:2]))
		
		while(ptr<nrow(a))
		{
			pos<-find.min.pos(d.mat)
			h.height[ptr]<-min(d.mat)
			h.merge[ptr,]<-lab[pos]
			this.labs<-sort(lab[pos])
			
			if(sum(this.labs < 0) == 2)
			{
				h.order<-c(h.order, this.labs)
				mem<--this.labs
			}else{
				if(sum(this.labs < 0) == 1)
				{
					this.members<--which(grps==this.labs[2])
					mem<-c(-this.members, -this.labs[1])
					to.insert<-max(which(h.order %in% this.members))
					if(to.insert == length(h.order))
					{
						h.order<-c(h.order, this.labs[1])
					}else{
						h.order<-c(h.order[1:to.insert], this.labs[1], h.order[(to.insert+1):length(h.order)])
					}
				}else{
					this.members.1<--which(grps==this.labs[1])
					this.members.2<--which(grps==this.labs[2])
					mem<--c(this.members.1, this.members.2)
					this.range.1<-range(which(h.order %in% this.members.1))
					this.range.2<-range(which(h.order %in% this.members.2))
					if(this.range.1[2] > this.range.2[2])
					{
						this.range.0<-this.range.2
						this.range.2<-this.range.1
						this.range.1<-this.range.0
					}
					if(this.range.2[1] == this.range.1[2]+1)
					{}else{
						range.between<-c(this.range.1[2]+1, this.range.2[1]-1)
						if(this.range.2[2] == length(h.order))
						{
							h.order<-c(h.order[1:this.range.1[2]], h.order[this.range.2[1]:this.range.2[2]], h.order[range.between[1]:range.between[2]])
						}else{
							h.order<-c(h.order[1:this.range.1[2]], h.order[this.range.2[1]:this.range.2[2]], h.order[range.between[1]:range.between[2]], h.order[(this.range.2[2]+1):length(h.order)])
						}
					}
				}
			}
			
			lab[pos[1]]<-ptr
			lab[pos[2]]<-NA
			a[pos[2],]<-rep(NA, m)
			grps[mem]<-ptr
			
			trav<-hamil(t(orig.a[mem,]), method.1=hamil.method)
			prof<-qnorm((1:m)/(m+1))[order(trav[[1]])]
			a[pos[1],]<-prof
			
			prof<-matrix(prof, nrow=1)
			
			all.given.curr<-scol.matrix.order(orig.a, prof)
			all.given.curr<-gene.specific.p(null.distr,all.given.curr)
			curr.given.all<-d.mat[,mem]
			
			d.mat[pos[1],]<-log10(all.given.curr)
			d.mat[pos[2],]<-NA
			
			uniq.grps<-unique(grps)
			all.given.curr<-qnorm(all.given.curr)
			grp.given.curr<-all.given.curr
			grp.given.curr[is.na(lab)]<-NA
			
			grp.given.curr[which(grps<0)]<-log10(pnorm(grp.given.curr[which(grps<0)], lower.tail=T))
			
			for(i in uniq.grps[uniq.grps>0])
			{
				this<-mean(all.given.curr[which(grps == i)])
				grp.given.curr[which(lab == i)]<-log10(pnorm(this, lower.tail=T))
			}
			
			if(length(mem) > 1)
			{
				curr.given.all<-qnorm(10^curr.given.all)
				curr.given.all<-apply(curr.given.all,1,mean)
				curr.given.all<-log10(pnorm(curr.given.all, lower.tail=T))
			}
			
			new.sim<-apply(cbind(grp.given.curr, curr.given.all), 1, min, na.rm=F)
			new.sim[is.na(new.sim)]<-0
			d.mat[pos[1], ] <- d.mat[, pos[1]] <- new.sim
			diag(d.mat)<-0
			d.mat[pos[2],]<-d.mat[,pos[2]]<-rep(0, n)
			
			ptr<-ptr+1
		}
		
		h.height.0<-h.height
		for(i in (n-1):2)
		{
			for(k in 1:2)
			{
				if(h.merge[i,k]>0)
				{
					if(h.height[h.merge[i,k]] >= h.height[i]) h.height[h.merge[i,k]]<-h.height[i]-0.01/nrow(array)
				}
			}
		}
		
		h$merge<-h.merge
		h$order<--h.order
		h$height<-h.height
		h$height.0<-h.height.0
		h$call[1]<-"NLHC"
		h$call[2]<-"SCOL matrix"
		h$method<-NULL
		
		return(reorg.tree(h))
	}else{
		rdmat.diff<-d.mat-t(d.mat)
		sel<-which(rdmat.diff > 0)
		d.mat[sel]<-t(d.mat)[sel]
		
		h<-hclust(as.dist(d.mat), method=method.traditional.hclust)
		return(h)			
	}
}

