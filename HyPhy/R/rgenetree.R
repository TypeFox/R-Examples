rgenetree <-
function(n,spec.phy,lams,mus,root=NULL,genetips=NULL,alltips=FALSE)
{
	if(length(lams)==1) lams<-rep(lams,length(mus))
	BsEs<-matrix(NA,dim(spec.phy$edge)[1],5)
	As<-mus/lams
	Exps<-exp((lams-mus)*spec.phy$edge.length)
	OnemUs<-(1-As)/(Exps-As)
	r0<-lams==mus
	if(any(r0)) OnemUs[r0]<-1/(lams[r0]*spec.phy$edge.length[r0]+1)
	BsEs[,5]<-OnemAUs<-Exps*OnemUs
	BsEs[,4]<-Us<-1-OnemUs

	for (i in dim(spec.phy$edge)[1]:1)
	{
		if (spec.phy$edge[i,2]<spec.phy$edge[1])  BsEs[i,3]<-0
		else BsEs[i,3]<-prod(BsEs[which(spec.phy$edge[,1]==spec.phy$edge[i,2]),2])

		Denom<-1-Us[i]*BsEs[i,3]
		BsEs[i,2]<-1-OnemAUs[i]*(1-BsEs[i,3])/Denom
		BsEs[i,1]<-1-OnemUs[i]/Denom
	}
	colnames(BsEs)<-c("B","E1","E2","U","OmAU")


	if (!is.null(genetips) && length(genetips)>1)
	{
		PGNgM<-PGgM<-PGM2M1gN<-PGgN<-vector("list",dim(spec.phy$edge)[1])
		gs<-rep(0,dim(spec.phy$edge)[1])

		for (i in dim(spec.phy$edge)[1]:1)
		{

			
			PGNgM[[i]]<-array(0,dim=rep(gs[i]+1,2))

			if (spec.phy$edge[i,2]<spec.phy$edge[1]){
				gs[i]<-genetips[spec.phy$edge[i,2]]
				PGgM[[i]]<-PGgN[[i]]<-rep(0,gs[i]+1)
				PGgN[[i]][gs[i]+1]<-1
			}
			else
			{
				nexts<-which(spec.phy$edge[,1]==spec.phy$edge[i,2])
				OnemEs.next<-1-BsEs[nexts,2]
				gs[i]<-sum(gs[nexts])

				PGgM[[i]]<-PGgN[[i]]<-rep(0,gs[i]+1)
				PGM2M1gN[[i]]<-array(0,dim=c(gs[i],gs[nexts])+1)
					## dimensions are N, M1, M2

				for (Ni in 0:gs[i]) #I believe that his will work for a 0, but not sure
				{
	
					for (M1 in max(0,Ni-gs[nexts[2] ]):min(Ni,gs[nexts[1] ]))
						for (M2 in (Ni-M1):min(Ni,gs[nexts[2] ]))
							PGM2M1gN[[i]][Ni+1,M1+1,M2+1]<-dbinom(M1,Ni,OnemEs.next[1]/(1-prod(BsEs[nexts,2])))*
								dbinom(M2+M1-Ni,M1,OnemEs.next[2])*
								PGgM[[nexts[1] ]][M1+1]*PGgM[[nexts[2] ]][M2+1]
					PGM2M1gN[[i]][Ni+1,,]<-cumsum(PGM2M1gN[[i]][Ni+1,,])
					PGgN[[i]][Ni+1]<-PGM2M1gN[[i]][Ni+1,,][length(PGM2M1gN[[i]][Ni+1,,])]	
				}
			}
	
			PGNgM[[i]]<-array(0,dim=c(gs[i],gs[i])+1)
			if (gs[i]==0) PGNgM[[i]][1,1]<-PGgM[[i]][1]<-1
			else for (Mi in 1:gs[i])
			{
				for (Ni in Mi:gs[i])
					PGNgM[[i]][Mi+1,Ni+1]<-dnbinom(Ni-Mi,Mi,1-BsEs[i,1])*PGgN[[i]][Ni+1]
				PGNgM[[i]][Mi+1,]<-cumsum(PGNgM[[i]][Mi+1,])
				PGgM[[i]][Mi+1]<-PGNgM[[i]][Mi+1,gs[i]+1]
			}
		
		}

		nexts<-which(spec.phy$edge[,1]==spec.phy$edge[1])
		OnemEs.next<-1-BsEs[nexts,2]

		gR<-sum(gs[nexts])
		PGgNR<-rep(0,gR+1)
		PGM2M1gNR<-array(0,dim=c(gR,gs[nexts])+1)
				## dimensions are N, M1, M2

		for (Ni in 1:gR) #I believe that his will work for a 0, but not sure
		{

			for (M1 in max(0,Ni-gs[nexts[2] ]):min(Ni,gs[nexts[1] ]))
				for (M2 in (Ni-M1):min(Ni,gs[nexts[2] ]))
					PGM2M1gNR[Ni+1,M1+1,M2+1]<-dbinom(M1,Ni,OnemEs.next[1]/(1-prod(BsEs[nexts,2])))*
						dbinom(M2+M1-Ni,M1,OnemEs.next[2])*
						PGgM[[nexts[1] ]][M1+1]*PGgM[[nexts[2] ]][M2+1]
			PGM2M1gNR[Ni+1,,]<-cumsum(PGM2M1gNR[Ni+1,,])
			PGgNR[Ni+1]<-PGM2M1gNR[Ni+1,,][length(PGM2M1gNR[Ni+1,,])]	
		}

		if (is.null(root)) root<-rep(1,gR)
		if (length(root)>1) PNgGR<-cumsum(root*PGgNR[-1])
	}
	else if (!is.null(genetips))
	{
		Ts<-alpha<-omega<-rep(NA,dim(spec.phy$edge)[1])
		PGNgM<-PGgM<-PGG1M2M1gN<-PGgN<-vector("list",dim(spec.phy$edge)[1])
		beta<-alltips*1
		T0<-length(spec.phy$tip.label)

		for (i in dim(spec.phy$edge)[1]:1)
		{
			if (spec.phy$edge[i,2]<spec.phy$edge[1]) Ts[i]<-1
			else{
				nexts<-which(spec.phy$edge[,1]==spec.phy$edge[i,2])
				OnemEs.next<-1-BsEs[nexts,2]

				Ts[i]<-sum(Ts[nexts])
			}

			alpha[i]<-Ts[i]*alltips
			omega[i]<-genetips+(Ts[i]-T0)*alltips
			PGgM[[i]]<-PGgN[[i]]<-matrix(0,omega[i]+1,omega[i]+1)
			PGNgM[[i]]<-array(0,dim=rep(omega[i]+1,3))
	
			if (spec.phy$edge[i,2]<spec.phy$edge[1])
				PGgN[[i]][cbind(alpha[i]:omega[i],alpha[i]:omega[i])+1]<-1
			else
			{
				PGG1M2M1gN[[i]]<-array(0,dim=c(omega[i],omega[nexts],omega[nexts[1]],omega[i])+1)
					## dimensions are N, M1, M2, G1, G for this node
				for (Ni in 1:omega[i])
				{

					for (M1 in max(beta,Ni-omega[nexts[2] ]):min(Ni,omega[nexts[1] ]))
						for (M2 in max(beta,Ni-M1):min(Ni,omega[nexts[2] ],omega[i]-M1))
						{
							PM2M1gN<-dbinom(M1,Ni,OnemEs.next[1]/(1-prod(BsEs[nexts,2])))*
								dbinom(M2+M1-Ni,M1,OnemEs.next[2])
							for(G1 in max(M1,alpha[nexts[1] ]):min(omega[nexts[1] ],omega[i]-M2))
								for (Gi in max(G1+M2,alpha[i]):omega[i])
									PGG1M2M1gN[[i]][Ni+1,M1+1,M2+1,G1+1,Gi+1]<-PM2M1gN*PGgM[[nexts[1] ]][M1+1,G1+1]*PGgM[[nexts[2] ]][M2+1,Gi+1-G1]
						}
					for (Gi in Ni:omega[i])
					{
						PGG1M2M1gN[[i]][Ni+1,,,,Gi+1]<-cumsum(PGG1M2M1gN[[i]][Ni+1,,,,Gi+1])
						PGgN[[i]][Ni+1,Gi+1]<-PGG1M2M1gN[[i]][Ni+1,,,,Gi+1][length(PGG1M2M1gN[[i]][Ni+1,,,,Gi+1])]
					}	
				}
				PGG1M2M1gN[[i]][1]<-PGgN[[i]][1]<-1-alltips
			}

			PGNgM[[i]]<-array(0,dim=c(omega[i],omega[i],omega[i])+1)
			for (Mi in 1:omega[i])
				for (Ni in Mi:omega[i])
					for (Gi in Ni:omega[i])
						PGNgM[[i]][Mi+1,Ni+1,Gi+1]<-dnbinom(Ni-Mi,Mi,1-BsEs[i,1])*PGgN[[i]][Ni+1,Gi+1]
			for (Mi in 1:omega[i])
				for (Gi in Mi:omega[i])
				{
					PGNgM[[i]][Mi+1,,Gi+1]<-cumsum(PGNgM[[i]][Mi+1,,Gi+1])
					PGgM[[i]][Mi+1,Gi+1]<-PGNgM[[i]][Mi+1,omega[i]+1,Gi+1]
				}
			PGNgM[[i]][1]<-PGgM[[i]][1]<-1-alltips
		
		}


		nexts<-which(spec.phy$edge[,1]==spec.phy$edge[1])
		OnemEs.next<-1-BsEs[nexts,2]

		alphaR<-T0
		PGgMR<-PGgNR<-matrix(0,genetips+1)

		PGG1M2M1gNR<-array(0,dim=c(genetips,omega[nexts],omega[nexts[1]])+1)
				## dimensions are N, M1, M2, G1, No G, couse know it for this node
		for (Ni in 1:genetips)
		{
			for (M1 in max(beta,Ni-omega[nexts[2] ]):min(Ni,omega[nexts[1] ]))
				for (M2 in max(beta,Ni-M1):min(Ni,omega[nexts[2] ]))
				{
					PM2M1gN<-dbinom(M1,Ni,OnemEs.next[1]/(1-prod(BsEs[nexts,2])))*
						dbinom(M2+M1-Ni,M1,OnemEs.next[2])
					for(G1 in max(M1,alpha[nexts[1] ]):omega[nexts[1] ])
							PGG1M2M1gNR[Ni+1,M1+1,M2+1,G1+1]<-PM2M1gN*PGgM[[nexts[1] ]][M1+1,G1+1]*PGgM[[nexts[2] ]][M2+1,genetips+1-G1]
				}
			PGG1M2M1gNR[Ni+1,,,]<-cumsum(PGG1M2M1gNR[Ni+1,,,])
			PGgNR[Ni+1]<-PGG1M2M1gNR[Ni+1,,,][length(PGG1M2M1gNR[Ni+1,,,])]
	
		}

		if (is.null(root)) root<-rep(1/genetips,genetips)
		if(length(root)>1 & length(root)<genetips) root<-c(root,rep(0,genetips-length(root)))
		if (length(root)>1) PNgGR<-cumsum(root*PGgNR[-1])
	}
	else if (alltips)
	{
		FS<-FE<-F2<-counts<-OmEE<-OmES<-vector("list",dim(spec.phy$edge)[1])

		for(i in dim(spec.phy$edge)[1]:1)
		{
			if (spec.phy$edge[i,2]<spec.phy$edge[1])
			{
				counts[[i]]<-c(0,1)
				OmEE[[i]]<-c(0,1)
				OmES[[i]]<-c(0,BsEs[i,5])
			}
			else 
			{
				nexts<-which(spec.phy$edge[,1]==spec.phy$edge[i,2])

				for (j in 1:length(OmES[[nexts[1]]]))
				{
					counts[[i]]<-c(counts[[i]],counts[[nexts[1]]][j]+counts[[nexts[2]]])
					OmEE[[i]]<-c(OmEE[[i]],1-(1-OmES[[nexts[1]]][j])*(1-OmES[[nexts[2]]]))
				}
				OmES[[i]]<-BsEs[i,5]*OmEE[[i]]/(1-BsEs[i,4]*(1-OmEE[[i]]))

			}
			FE[[i]]<-1-OmEE[[i]]/OmEE[[i]][length(OmEE[[i]])]
			FS[[i]]<-1-OmES[[i]]/OmES[[i]][length(OmES[[i]])]
			F2[[i]]<- (-1)^counts[[i]]
		}

		nexts<-which(spec.phy$edge[,1]==spec.phy$edge[1])
		countR<-OmER<-NULL	
		for (j in 1:length(OmES[[1]]))
		{
			countR<-c(countR,counts[[1]][j]+counts[[nexts[2]]])
			OmER<-c(OmER,1-(1-OmES[[1]][j])*(1-OmES[[nexts[2]]]))
		}
		FR<-1-OmER/OmER[length(OmER)]
		F2R<- (-1)^countR


		if (is.null(root)) stop("Must set either root or genetips to value other than Null")
		else if (length(root)>1) 
		{
			choose.R<-sum(F2R*FR)*root[1]
			for (i in 2:length(root))
				choose.R[i]<-choose.R[i-1]+root[i]*sum(F2R*FR^(i))
		}

	}



	out<-list()
	for (i in 1:n)
	{

		MsNs<-array(NA,dim=dim(spec.phy$edge))
		if (!is.null(genetips) && length(genetips)>1)
		{
			if (length(root)==1) RootN<-root ##Modify when complicate root
			else RootN<-which(PNgGR>=runif(1,0,PNgGR[gR]))[1]

			chooses<-which(array(PGM2M1gNR[RootN+1,,],dim=dim(PGM2M1gNR)[-1])>runif(1,0,PGgNR[RootN+1]),arr.ind=T)
			if (!is.null(dim(chooses))) chooses<-chooses[1,]

			nexts<-which(spec.phy$edge[,1]==spec.phy$edge[1])
			MsNs[nexts[1],1]<-chooses[1]-1
			MsNs[nexts[2],1]<-chooses[2]-1

			for (j in 1:dim(spec.phy$edge)[1])
			{
				if (spec.phy$edge[j,2]<spec.phy$edge[1]) MsNs[j,2]<-genetips[spec.phy$edge[j,2]]
				else
				{
					nexts<-which(spec.phy$edge[,1]==spec.phy$edge[j,2])
					if (MsNs[j,1]==0)
						MsNs[j,2]<-MsNs[nexts,1]<-0

					else{
						MsNs[j,2]<-which(PGNgM[[j]][MsNs[j,1]+1,]>runif(1,0,PGgM[[j]][MsNs[j,1]+1]))[1]-1
	
						chooses<-which(matrix(PGM2M1gN[[j]][MsNs[j,2]+1,,],gs[nexts[1]]+1,gs[nexts[2]]+1)>runif(1,0,PGgN[[j]][MsNs[j,2]+1]),arr.ind=T)
						if (!is.null(dim(chooses))) chooses<-chooses[1,]
	
						MsNs[nexts[1],1]<-chooses[1]-1
						MsNs[nexts[2],1]<-chooses[2]-1
					}
				}
			}
		}
		else if (!is.null(genetips))
		{
			gammas<-rep(NA,dim(spec.phy$edge)[1])

			if (length(root)==1) RootN<-root ##Modify when complicate root
			else RootN<-which(PNgGR>=runif(1,0,PNgGR[genetips]))[1]

			chooses<-which(PGG1M2M1gNR[RootN+1,,,]>runif(1,0,PGgNR[RootN+1]),arr.ind=T)
			if (!is.null(dim(chooses))) chooses<-chooses[1,]

			nexts<-which(spec.phy$edge[,1]==spec.phy$edge[1])
			gammas[nexts[1]]<-chooses[3]-1
			gammas[nexts[2]]<-genetips+1-chooses[3]
			MsNs[nexts[1],1]<-chooses[1]-1
			MsNs[nexts[2],1]<-chooses[2]-1

			for (j in 1:dim(spec.phy$edge)[1])
			{
				if (spec.phy$edge[j,2]<spec.phy$edge[1]) MsNs[j,2]<-gammas[j]
				else
				{
					nexts<-which(spec.phy$edge[,1]==spec.phy$edge[j,2])
					if (MsNs[j,1]==0)
						MsNs[j,2]<-gammas[nexts]<-MsNs[nexts,1]<-0

					else{
						MsNs[j,2]<-which(PGNgM[[j]][MsNs[j,1]+1,,gammas[j]+1]>runif(1,0,PGgM[[j]][MsNs[j,1]+1,gammas[j]+1]))[1]-1

						chooses<-which(PGG1M2M1gN[[j]][MsNs[j,2]+1,,,,gammas[j]+1]>runif(1,0,PGgN[[j]][MsNs[j,2]+1,gammas[j]+1]),arr.ind=T)
						if (!is.null(dim(chooses))) chooses<-chooses[1,]

						gammas[nexts[1]]<-chooses[3]-1
						gammas[nexts[2]]<-gammas[j]+1-chooses[3]
						MsNs[nexts[1],1]<-chooses[1]-1
						MsNs[nexts[2],1]<-chooses[2]-1
					}
				}

			}
		}
		else if (alltips) 
		{
			if (length(root)==1) RootN<-root
			else RootN<-which(choose.R>=runif(1,0,choose.R[length(root)]))[1]

			pp<-runif(1,0,sum(F2R*FR^RootN))
			nexts<-which(spec.phy$edge[,1]==spec.phy$edge[1])
			OnemEs.next<-1-BsEs[nexts,2]
		
			MsNs[1,1]<-1
			MsNs[nexts[2],1]<-max(0,RootN-2)
			while(pp>0)
			{	
				MsNs[nexts[2],1]<-MsNs[nexts[2],1]+1
				if(MsNs[nexts[2],1]>RootN ){
					MsNs[1,1]<-MsNs[1,1]+1
					MsNs[nexts[2],1]<-max(1,RootN-MsNs[1,1])	
				}
				pp<-pp-dbinom(MsNs[nexts[1],1],RootN,OnemEs.next[1]/(1-prod(BsEs[nexts,2])))*
					dbinom(MsNs[nexts[2],1]+MsNs[nexts[1],1]-RootN,MsNs[nexts[1],1],OnemEs.next[2])*
					sum(F2[[nexts[1] ]]*FS[[nexts[1] ]]^MsNs[nexts[1],1])*
					sum(F2[[nexts[2] ]]*FS[[nexts[2] ]]^MsNs[nexts[2],1])
			}

			for (j in 1:dim(spec.phy$edge)[1])	
			{
	
				pp<-runif(1,0,sum(F2[[j]]*FS[[j]]^MsNs[j,1]))

				MsNs[j,2]<-MsNs[j,1]-1
				while(pp>0)
				{

					MsNs[j,2]<-MsNs[j,2]+1
					pp<-pp-sum(F2[[j]]*FE[[j]]^MsNs[j,2])*dnbinom(MsNs[j,2]-MsNs[j,1],MsNs[j,1],1-BsEs[j,1])
				}

				if (spec.phy$edge[j,2]>spec.phy$edge[1])
				{
					pp<-runif(1,0,sum(F2[[j]]*FE[[j]]^MsNs[j,2]))
					nexts<-which(spec.phy$edge[,1]==spec.phy$edge[j,2])
					OnemEs.next<-1-BsEs[nexts,2]
		
					MsNs[nexts[1],1]<-1
					MsNs[nexts[2],1]<-max(0,MsNs[j,2]-2)
					while(pp>0)
					{
						MsNs[nexts[2],1]<-MsNs[nexts[2],1]+1
						if(MsNs[nexts[2],1]>MsNs[j,2] ){
							MsNs[nexts[1],1]<-MsNs[nexts[1],1]+1
							MsNs[nexts[2],1]<-max(1,MsNs[j,2]-MsNs[nexts[1],1])	
						}
						pp<-pp-dbinom(MsNs[nexts[1],1],MsNs[j,2],OnemEs.next[1]/(1-BsEs[j,3]))*
							dbinom(sum(MsNs[nexts,1])-MsNs[j,2],MsNs[nexts[1],1],OnemEs.next[2])*
							sum(F2[[nexts[1] ]]*FS[[nexts[1] ]]^MsNs[nexts[1],1])*
							sum(F2[[nexts[2] ]]*FS[[nexts[2] ]]^MsNs[nexts[2],1])
					}
				}
			}
		}
		else
		{
			if (is.null(root)) stop("Must set root or genetips to value other than NULL")
			if (length(root)==1) RootN<-root
			else RootN<-which(root>=runif(1,0,root[length(root)]))[1]

			nexts<-which(spec.phy$edge[,1]==spec.phy$edge[1])
			OnemEs.next<-1-BsEs[nexts,2]
				
			MsNs[nexts[1],1]<-rbinom(1,RootN,(OnemEs.next[1])/(1-prod(BsEs[nexts,2])))
			MsNs[nexts[2],1]<-rbinom(1,MsNs[nexts[1],1],OnemEs.next[2])+RootN-MsNs[nexts[1],1]

			for(j in 1:dim(spec.phy$edge)[1])
			{
				if (MsNs[j,1]==0) MsNs[j,2]<-0
				else MsNs[j,2]<-rnbinom(1,MsNs[j,1],1-BsEs[j,1])+MsNs[j,1]
				
				nexts<-which(spec.phy$edge[,1]==spec.phy$edge[j,2])
				if (length(nexts)>1){

					if (MsNs[j,2]==0) MsNs[nexts,1]<-0
					else{
						OnemEs.next<-1-BsEs[nexts,2]
					
						MsNs[nexts[1],1]<-rbinom(1,MsNs[j,2],(OnemEs.next[1])/(1-BsEs[j,3]))
						MsNs[nexts[2],1]<-rbinom(1,MsNs[nexts[1],1],OnemEs.next[2])+MsNs[j,2]-MsNs[nexts[1],1]
					}
				}
			}
		}

		edge<-rbind(c(0,spec.phy$edge[1]),spec.phy$edge)
		edge.length<-c(1/0,spec.phy$edge.length)
		MsNs<-rbind(c(1,RootN),MsNs)

		edge1<-matrix(NA,sum(MsNs[edge[,2]<spec.phy$edge[1],2])*2-1,2)
		loose<-matrix(FALSE,dim(edge)[1],dim(edge1)[1])
		unused<-dim(edge1)[1]
		tips<-NULL
		nextnode<- -1
		for(j in dim(edge)[1]:1)
		{
			if (MsNs[j,2]==0) loose[j,]<-F
			else if (edge[j,2]<spec.phy$edge[1])
			{
				using<-unused+1-(MsNs[j,2]:1)
				edge1[using,2]<-length(tips)+(1:MsNs[j,2])

				loose[j,using]<-T
				unused<-unused-MsNs[j,2]

				tips<-c(tips,rep(spec.phy$tip.label[edge[j,2]],MsNs[j,2]))
			}
			else
			{
				nexts<-which(edge[,1]==edge[j,2])
				combN<-sum(MsNs[nexts,1])-MsNs[j,2]
				loose[j,]<-loose[nexts[1],]|loose[nexts[2],]
				if (combN>0)
				{
					for (k in 1:2)
					{
						if (sum(loose[nexts[k],])==1) using<-which(loose[nexts[k],])
						else using<-sample(which(loose[nexts[k],]),combN)

						edge1[using,1]<-nextnode +1 -(1:combN)
						loose[j,using]<-F
					} 
					using<-unused+1-(combN:1)
					edge1[using,2]<-nextnode +1 -(1:combN)

					loose[j,using]<-T
					unused<-unused-combN
					nextnode<- nextnode -combN
			
				}
			}
			combN<-MsNs[j,2]-MsNs[j,1]
			if (combN>0) for (k in 1:combN)
			{
				using<-sample(which(loose[j,]),2)
				edge1[using,1]<-nextnode
				loose[j,using]<-F
				edge1[unused,2]<-nextnode
				loose[j,unused]<-T
				unused<-unused-1
				nextnode<- nextnode -1
			}

		}

		fixedge<-function(rown,firstn,nextn)
		{
			if (edge1[rown,2]>0) return(matrix(c(firstn,edge1[rown,2]),1,2))
			edget<-c(firstn,nextn)
			nexts<-which(edge1[,1]==edge1[rown,2])
			edget<-rbind(edget,fixedge(nexts[1],nextn,nextn+1))
			matrix(rbind(edget,
				fixedge(nexts[2],nextn,max(c(nextn,edget))+1)),ncol=2)
		}
		edge2<-fixedge(1,0,length(tips)+1)[-1,]
		out[[i]]<-list(edge=edge2,tip.label=tips,Nnode=edge2[1]-2)
		class(out[[i]])<-"phylo"

	}
	if (n==1) return(out[[1]])
	class(out)<-"multiPhylo"
	out
}
