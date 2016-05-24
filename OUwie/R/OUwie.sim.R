##OUwie Simulator##

#written by Jeremy M. Beaulieu

#Simulates the Hansen model of continuous characters evolving under discrete selective 
#regimes. The input is a tree of class "phylo" that has the regimes as internal node labels 
#and a data file the contains the regime states for each species. The trait file must be in 
#the following order: Species names then Regime. The user must specify the parameters values
#for each simulation (i.e. alpha, sigma.sq, theta0, theta). 

##The following examples assume 2 selective regimes and different models can be specified:
#single rate Brownian motion BM1: alpha=c(0,0); sigma.sq=c(0.9); theta0=0; theta=c(0,0)
#two rate Brownian motion BMS: alpha=c(0,0); sigma.sq=c(0.45,.9); theta0=0; theta=c(0,0)
#global OU (OU1): alpha=c(0.1,0.1); sigma.sq=c(0.9,0.9); theta0=1; theta=c(1,1)
#normal OU (OUSM): alpha=c(0.1,0.1); sigma.sq=c(0.9,0.9); theta0=0; theta=c(1,2)
#multiple sigmas (OUSMV): alpha=c(0.1,0.1); sigma.sq=c(0.45,0.9); 
#multiple alphas (OUSMA): alpha=c(0.5,0.1); sigma.sq=c(0.9,0.9); theta0=0; theta=c(1,2)
#multiple alphas and sigmas (OUSMVA): alpha=c(0.5,0.1); sigma.sq=c(0.45,0.9); theta0=0; theta=c(1,2)

OUwie.sim <- function(phy, data=NULL, simmap.tree=FALSE, scaleHeight=FALSE, alpha, sigma.sq, theta0, theta, mserr="none"){

	if(simmap.tree==FALSE){
		#This is annoying, but the second column has to be in there twice otherwise, error.
        if(mserr=="none"){
            data <- data.frame(data[,2], data[,2], row.names=data[,1])
        }
        if(mserr=="known"){
            data <- data.frame(data[,2], data[,3], row.names=data[,1])
        }

		data <- data[phy$tip.label,]
		
		n=max(phy$edge[,1])
		ntips=length(phy$tip.label)
		
		int.states<-factor(phy$node.label)
		phy$node.label=as.numeric(int.states)
		tip.states<-factor(data[,1])
		data[,1]<-as.numeric(tip.states)
		tot.states<-factor(c(phy$node.label,as.character(data[,1])))
		k<-length(levels(tot.states))
		
		regime=matrix(rep(0,(n-1)*k), n-1, k)

		#Obtain root state and internal node labels
		root.state<-phy$node.label[1]
		int.state<-phy$node.label[-1]
		
		#New tree matrix to be used for subsetting regimes
		edges=cbind(c(1:(n-1)),phy$edge,nodeHeights(phy))
		if(scaleHeight==TRUE){
			edges[,4:5]<-edges[,4:5]/max(nodeHeights(phy))
		}
		edges=edges[sort.list(edges[,3]),]

		mm<-c(data[,1],int.state)

		regime <- matrix(0,nrow=length(mm),ncol=length(unique(mm)))
		#Generates an indicator matrix from the regime vector
		for (i in 1:length(mm)) {
			regime[i,mm[i]] <- 1 
		}
		#Finishes the edges matrix
		edges=cbind(edges,regime)
		
		#Resort the edge matrix so that it looks like the original matrix order
		edges=edges[sort.list(edges[,1]),]

		oldregime=root.state
		
		alpha=alpha
		alpha[alpha==0] = 1e-10
		sigma=sqrt(sigma.sq)
		theta=theta

		x <- matrix(0, n, 1)
		TIPS <- 1:ntips
		ROOT <- ntips + 1L
		x[ROOT,] <- theta0
		
		for(i in 1:length(edges[,1])){
			anc = edges[i,2]
			desc = edges[i,3]
			oldtime=edges[i,4]
			newtime=edges[i,5]
			if(anc%in%edges[,3]){
				start=which(edges[,3]==anc)
				oldregime=which(edges[start,6:(k+5)]==1)
			}
			else{
				#For the root:
				oldregime=oldregime
			}	
			newregime=which(edges[i,6:(k+5)]==1)

			if(oldregime==newregime){
				x[edges[i,3],]=x[edges[i,2],]*exp(-alpha[oldregime]*(newtime-oldtime))+(theta[oldregime])*(1-exp(-alpha[oldregime]*(newtime-oldtime)))+sigma[oldregime]*rnorm(1,0,1)*sqrt((1-exp(-2*alpha[oldregime]*(newtime-oldtime)))/(2*alpha[oldregime]))
			}
			else{
				halftime=newtime+((oldtime-newtime)/2)
				epoch1=x[edges[i,2],]*exp(-alpha[oldregime]*(halftime-oldtime))+(theta[oldregime])*(1-exp(-alpha[oldregime]*(halftime-oldtime)))+sigma[oldregime]*rnorm(1,0,1)*sqrt((1-exp(-2*alpha[oldregime]*(halftime-oldtime)))/(2*alpha[oldregime]))
				oldtime=halftime
				newtime=newtime
				x[edges[i,3],]=epoch1*exp(-alpha[newregime]*(newtime-oldtime))+(theta[newregime])*(1-exp(-alpha[newregime]*(newtime-oldtime)))+sigma[newregime]*rnorm(1,0,1)*sqrt((1-exp(-2*alpha[newregime]*(newtime-oldtime)))/(2*alpha[newregime]))
			}
		}
		
		sim.dat<-matrix(,ntips,3)
		sim.dat<-data.frame(sim.dat)
		
		sim.dat[,1]<-phy$tip.label
		sim.dat[,2]<-data[,1]
		sim.dat[,3]<-x[TIPS]
        
        if(mserr == "known"){
            for(i in TIPS){
                sim.dat[i,3] <- rnorm(1,sim.dat[i,3],data[i,2])
            }
        }
		colnames(sim.dat)<-c("Genus_species","Reg","X")
        
	}
	if(simmap.tree==TRUE){
		n=max(phy$edge[,1])
		ntips=length(phy$tip.label)
		
		k=length(colnames(phy$mapped.edge))

		regimeindex<-colnames(phy$mapped.edge)
		##Begins the construction of the edges matrix -- similar to the ouch format##
		#Makes a vector of absolute times in proportion of the total length of the tree
		branch.lengths=rep(0,(n-1))
		branch.lengths[(ntips+1):(n-1)]=branching.times(phy)[-1]/max(branching.times(phy))
		
		#Obtain root state and internal node labels
		root.state<-which(colnames(phy$mapped.edge)==names(phy$maps[[1]][1]))
		
		#New tree matrix to be used for subsetting regimes
		edges=cbind(c(1:(n-1)),phy$edge,nodeHeights(phy))
		if(scaleHeight==TRUE){
			edges[,4:5]<-edges[,4:5]/max(nodeHeights(phy))
		}
		edges=edges[sort.list(edges[,3]),]
		
		#Resort the edge matrix so that it looks like the original matrix order
		edges=edges[sort.list(edges[,1]),]
		
		oldregime=root.state
		oldtime=0
		
		alpha=alpha
		sigma=sqrt(sigma.sq)
		theta=theta
		
		n.cov=matrix(rep(0,n*n), n, n)
		nodecode=matrix(c(ntips+1,1),1,2)
		
		x <- matrix(0, n, 1)
		TIPS <- 1:ntips
		ROOT <- ntips + 1L
		x[ROOT,] <- theta0
		
		for(i in 1:length(edges[,1])){
			
			if(scaleHeight==TRUE){
				currentmap<-phy$maps[[i]]/max(nodeHeights(phy))
			}
			else{
				currentmap<-phy$maps[[i]]
			}
			oldtime=edges[i,4]
			
			if(length(phy$maps[[i]])==1){
				regimeduration<-currentmap[1]
				newtime<-oldtime+regimeduration
				regimenumber<-which(colnames(phy$mapped.edge)==names(currentmap)[1])
				x[edges[i,3],]=x[edges[i,2],]*exp(-alpha[regimenumber]*(newtime-oldtime))+(theta[regimenumber])*(1-exp(-alpha[regimenumber]*(newtime-oldtime)))+sigma[regimenumber]*rnorm(1,0,1)*sqrt((1-exp(-2*alpha[regimenumber]*(newtime-oldtime)))/(2*alpha[regimenumber]))
			}
			if(length(phy$maps[[i]])>1){
				regimeduration<-currentmap[1]
				newtime<-oldtime+regimeduration
				regimenumber<-which(colnames(phy$mapped.edge)==names(currentmap)[1])
				x[edges[i,3],]=x[edges[i,2],]*exp(-alpha[regimenumber]*(newtime-oldtime))+(theta[regimenumber])*(1-exp(-alpha[regimenumber]*(newtime-oldtime)))+sigma[regimenumber]*rnorm(1,0,1)*sqrt((1-exp(-2*alpha[regimenumber]*(newtime-oldtime)))/(2*alpha[regimenumber]))				
				oldtime<-newtime
				for (regimeindex in 2:length(currentmap)){
					regimeduration<-currentmap[regimeindex]
					newtime<-oldtime+regimeduration
					regimenumber<-which(colnames(phy$mapped.edge)==names(currentmap)[regimeindex])
					x[edges[i,3],]=x[edges[i,3],]*exp(-alpha[regimenumber]*(newtime-oldtime))+(theta[regimenumber])*(1-exp(-alpha[regimenumber]*(newtime-oldtime)))+sigma[regimenumber]*rnorm(1,0,1)*sqrt((1-exp(-2*alpha[regimenumber]*(newtime-oldtime)))/(2*alpha[regimenumber]))
					oldtime<-newtime
					newregime<-regimenumber
				}
			}
		}
		
		sim.dat<-matrix(,ntips,2)
		sim.dat<-data.frame(sim.dat)
		
		sim.dat[,1]<-phy$tip.label
		sim.dat[,2]<-x[TIPS,]
        if(mserr == "known"){
            for(i in 1:TIPS){
                sim.dat[i,2] <- rnorm(1,sim.dat[i,2],data[i,2])
            }
        }
		colnames(sim.dat)<-c("Genus_species","X")
	}
	sim.dat
}



