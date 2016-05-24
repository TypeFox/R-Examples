#Expected value generator for the TrendyM and TrendMS models

#written by Jeremy M. Beaulieu

expected.trendy<-function(phy, edges, Rate.mat, root.state, simmap.tree=FALSE, scaleHeight=FALSE, root.value){
	
	n=max(phy$edge[,1])
	ntips=length(phy$tip.label)
	if(is.null(root.state)) {
		root.state<-which(edges[dim(edges)[1],]==1)-5
		edges<-edges[-1*dim(edges)[1],]
	}
	if(simmap.tree==TRUE){
		k=length(colnames(phy$mapped.edge))
	}
	if(simmap.tree==FALSE){
		mm<-dim(edges)
		k<-length(6:mm[2])
	}
	pp <- prop.part(phy)
	edges = edges
	oldregime=root.state
	nodevar=rep(0,max(edges[,3]))
	
	#Parameters:
	alpha=Rate.mat[1,]
	sigma=Rate.mat[2,]
	trend=Rate.mat[3,]
	
	E <- matrix(0,ntips,k)
	for(j in 1:k){
		n.cov=matrix(0, n, 1)			
		#Weights for each species per regime
		for(i in 1:length(edges[,1])){
			anc = edges[i, 2]
			oldtime=edges[i,4]
			newtime=edges[i,5]
			if(simmap.tree==TRUE){
				if(scaleHeight==TRUE){
					currentmap<-phy$maps[[i]]/max(nodeHeights(phy))
				}
				else{
					currentmap<-phy$maps[[i]]
				}					
			}
			if(simmap.tree==TRUE){
				nodevar[i]=0
				for (regimeindex in 1:length(currentmap)){
					regimeduration <- currentmap[regimeindex]
					newtime <- oldtime + regimeduration
					regimenumber <- which(colnames(phy$mapped.edge)==names(currentmap)[regimeindex])
					if (regimenumber == j){
						nodevar[i] <- trend[regimenumber]*(newtime - oldtime)
					}
					else{
						nodevar[i]=nodevar[i]
					}
					oldtime <- newtime
				}
			}
			if(simmap.tree==FALSE){
				if(anc%in%edges[,3]){
					start=which(edges[,3]==anc)
					oldregime=which(edges[start,6:(k+5)]==1)
				}
				else{
					#For the root:
					oldregime=root.state
				}	
				newregime=which(edges[i,6:(k+5)]==1)
				if(oldregime==newregime){
					if(newregime==j){
						nodevar[i] = trend[oldregime]*(newtime - oldtime)
					}
					else{
						nodevar[i]=0
					}
				}
				else{
					halftime=newtime-((newtime-oldtime)/2)
					epoch1 = trend[oldregime]*(newtime - oldtime)
					oldtime=halftime
					newtime=newtime
					epoch2=trend[newregime]*(newtime - oldtime)
					if(oldregime==j){
						nodevar[i]=epoch1
					}
					if(newregime==j){
						nodevar[i]=epoch2
					}
					if(!newregime==j && !oldregime==j){
						nodevar[i] = 0
					}
				}
			}
			n.cov[edges[i,3],]=nodevar[i]
		}
		e.piece<-mat.gen(phy,n.cov,pp)
		E[1:(ntips),j]<-diag(e.piece)
	}
	#The expected values are E[a] = root _ time_1 *trend_1 + time_2 * trend_2 + ...
	E <- rowSums(E) + root.value
	E
}

##Matrix generating function taken from vcv.phylo in ape:
mat.gen<-function(phy,piece.wise,pp){
	phy <- reorder(phy, "pruningwise")
    n <- length(phy$tip.label)
    anc <- phy$edge[,1]
    des <- phy$edge[,2]
    ep <- piece.wise[,1]
	comp <- numeric(n + phy$Nnode)
    mat <- matrix(0, n, n)
	
	for (i in length(anc):1) {
        focal <- comp[anc[i]]
        comp[des[i]] <- focal + ep[des[i]]
        j <- i - 1L
        while (anc[j] == anc[i] && j > 0) {
            left <- if (des[j] > n) pp[[des[j] - n]] else des[j]
            right <- if (des[i] > n) pp[[des[i] - n]] else des[i]
            mat[left, right] <- mat[right, left] <- focal
            j <- j - 1L
        }
    }
    diag.elts <- 1 + 0:(n - 1)*(n + 1)
    mat[diag.elts] <- comp[1:n]
	
	mat
}


