spacodi.calc <-function(sp.plot, phy = NULL, sp.traits = NULL, all.together=TRUE, prune=TRUE, pairwise=FALSE,...){
	if(!missing(phy) && !missing(sp.traits)) stop("Please supply either a phylogeny or a trait matrix, but not both.")

	rematrix=function(arr, names){
		n=sqrt(length(arr))
		counter=1
		M=matrix(0,n,n)
		for(i in 1:n){
			j=i+1
			while(j<=n) {
				M[i,j]=arr[counter]
				counter=counter+1
				j=j+1
			}
		}
		M[lower.tri(M)]=t(M)[lower.tri(t(M))]
		dimnames(M)[[1]]<-dimnames(M)[[2]]<-names
		return(M)
	}
	
	# determine type of abundance
	stripped=unname(unlist(c(sp.plot)))
	if(all(!is.na(match(stripped, c(0,1))))) {
		abundt = 0		# presence|absence
	} else if(all(!is.na(match(range(stripped), c(0,1)))) && length(unique(stripped))>2) {
		abundt = 1		# relative abundance
	} else {
		abundt = 2		# abundance (n.individuals)
	}
	
	# iter is 1 unless more than a single trait is being used for separate analyses
	iter=1
	traits.tmp<-distmat<-NULL
	
	# INTERPRET whether to compute trait diversity or phylogenetic diversity & prepare data
	if(!missing(phy)) {
		# distmat is phylogenetic: Pst
		sp.data=match.spacodi.data(sp.plot=sp.plot, phy=phy, prune=prune, ...)
		sp.plot=sp.data$sp.plot
		phy=sp.data$sp.tree
		if(check.distmat(phy)) distmat <- phy else distmat <- cophen(phy)
		
	} else if(missing(sp.traits) & missing(phy)) {
		# distmat is null: Ist
		distmat <- matrix(1, ncol = nrow(sp.plot), nrow = nrow(sp.plot))
		
	} else if(!missing(sp.traits)){
		# distmat is trait-based: Tst
		if(class(sp.traits)=="dist") sp.traits=as.matrix(sp.traits)
		if(ncol(sp.traits)==1) all.together=TRUE
		if(all(is.null(names(sp.traits)))) names(sp.traits)=paste("trt",seq(1:ncol(sp.traits)),sep="")	
		if(all.together==TRUE | check.distmat(sp.traits)) {
			sp.data=match.spacodi.data(sp.plot=sp.plot, sp.traits=sp.traits, prune=prune, ...)
			sp.plot=sp.data$sp.plot
			sp.traits=sp.data$sp.traits
			if(check.distmat(sp.traits)) distmat=sp.traits else distmat=as.matrix(dist(sp.traits))
		} else {	
			iter=ncol(sp.traits)
			traits.tmp=lapply(1:ncol(sp.traits), function(x) {
				   trait.tt<-data.frame(sp.traits[,x])
				   row.names(trait.tt)<-row.names(sp.traits)
				   sp.data<-match.spacodi.data(sp.plot=sp.plot, sp.traits=trait.tt, prune=prune, ...)
				   distmat <- as.matrix(dist(sp.data$sp.traits))
				   return(list(distmat=distmat, sp.plot=sp.data$sp.plot))
				   }
				   )
		}
	} else if(is.null(distmat)){ 
		stop("Cannot decipher input object(s).")
	}
	
	if(is.null(names(sp.plot))) pnames=paste("plt",seq(1,ncol(sp.plot)),sep=".") else pnames=names(sp.plot)

	
	# PREPARE output
	gen.out=list()
	prw.out=list()
	
	for(tt in 1:iter) {	
		if(is.null(traits.tmp)) {
			sp.plot <-	as.matrix(sp.plot)
		} else {
			sp.plot <-	as.matrix(traits.tmp[[tt]]$sp.plot)
			distmat <-	as.matrix(traits.tmp[[tt]]$distmat)
		}
		
		diag(distmat) <- 0
		np <- ncol(sp.plot)
		out <- NA
		out <- .C("spacodi", 
				  np = as.integer(np),
				  ns = as.integer(nrow(sp.plot)),
				  sp.plot = as.double(as.vector(sp.plot)),
				  distmat = as.double(as.vector(as.matrix(distmat))),
				  abundtype = as.integer(abundt),
				  Ndclass = as.integer(0),
				  dclass = as.double(c(0,0)),
				  Ist = 0,
				  Pst = 0,
				  Bst = 0,
				  PIst = 0,
				  pairwiseIst = as.double(as.vector(matrix(0,np,np))),
				  pairwisePst = as.double(as.vector(matrix(0,np,np))),
				  pairwiseBst = as.double(as.vector(matrix(0,np,np))),
				  pairwisePIst = as.double(as.vector(matrix(0,np,np))),
				  PACKAGE="spacodiR"
				  )
		
		# compile results for phylogenetic turnover
		if(missing(phy) & missing(sp.traits)) {
			r.out=as.numeric(out[8])
			names(r.out)="Ist"
			gen.out[[tt]]=as.data.frame(t(r.out))
			if(pairwise) {
				prw.out[[tt]]=list(rematrix(unlist(out[12]),pnames))
				names(prw.out[[tt]])<-paste("pairwise","Ist",sep=".")
			}
		} else if(!missing(phy)){
			r.out=as.numeric(c(out[8:11]))
			names(r.out)<-c("Ist","Pst","Bst","PIst")
			gen.out[[tt]]=as.data.frame(t(r.out))
			if(pairwise) {
				prw.out[[tt]]=lapply(out[12:15], function(x) rematrix(x,pnames))
				names(prw.out[[tt]])<-paste("pairwise",c("Ist","Pst","Bst","PIst"),sep=".")
			}
		} else if(!missing(sp.traits)) {
			r.out=as.numeric(c(out[8:11]))
			names(r.out)=c("Ist","Tst","Ust","TAUst")
			gen.out[[tt]]=as.data.frame(t(r.out))
			if(pairwise) {
				prw.out[[tt]]=lapply(out[12:15], function(x) rematrix(x,pnames))
				names(prw.out[[tt]])<-paste("pairwise",c("Ist","Tst","Ust","TAUst"),sep=".")
			}
		} 
	} 
	if(iter>1 & !all.together) {
		RES=lapply(1:iter, function(x) {
				   if(pairwise) return(c(gen.out[[x]], prw.out[[x]])) else return(c(gen.out[[x]]))
				   }
				   )
		names(RES)<-names(sp.traits)
	} else {
		prw.out=unlist(prw.out,recursive=FALSE)
		gen.out=unlist(gen.out,recursive=FALSE)
		if(pairwise) RES=c(gen.out,prw.out) else RES=gen.out
	}
	
	return(unlist(list(RES,sp.plot=list(sp.plot),sp.tree=list(phy),sp.traits=list(sp.traits)),recursive=FALSE)) 

} 

		
