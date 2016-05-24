TrackCluster <- function(List,Selection,nrclusters=NULL,followMaxComps=FALSE,followClust=TRUE,fusionsLog=TRUE,WeightClust=TRUE,names=NULL,SelectionPlot=TRUE,Table=TRUE,CompleteSelectionPlot=FALSE,ClusterPlot=FALSE,cols=NULL,legendposx=0.5,legendposy=2.4,plottype="new",location=NULL){
	
	ClusterDistribution.2<-function(List,Selection,nrclusters,followMaxComps,followClust,fusionsLog,WeightClust,names){
		FoundClusters=list()	
		FoundCl=list()
		Matrix=ReorderToReference(List,nrclusters,fusionsLog,WeightClust,names)
		
		ClusterNumber=NULL
		if(class(Selection)=="numeric"){
			ClusterNumber=Selection
			Selection=colnames(Matrix)[which(Matrix[1,]==Selection)]
		}
		
		
		if(followClust==TRUE){			
			cluster.interest=NULL
			m=1
			while(m<=dim(Matrix)[1] & is.null(cluster.interest)){
				clustersfound=unique(Matrix[m,which(colnames(Matrix)%in%Selection)])
				if(length(clustersfound)==1){
					cluster.interest=clustersfound
				}
				else{
					m=m+1
				}
			}
			if(is.null(cluster.interest)){
				message("This selection is not found to be part of a cluster. FollowMaxComps will be put to true and the function will proceed")
				followMaxComps=TRUE
				followClust=FALSE
			}
			
		}	

		for(i in 1:dim(Matrix)[1]){
			
			
			temp=list()
			temp[[1]]=Selection
			names(temp)[1]="Selection"
			
			clusternumbers=unique(Matrix[i,which(colnames(Matrix)%in%Selection)])
			
			nr.clusters=length(clusternumbers)
			temp[[2]]=nr.clusters
			names(temp)[2]="nr.clusters"
			
			min.together=min(table(Matrix[i,which(colnames(Matrix)%in%Selection)]))
			max.together=max(table(Matrix[i,which(colnames(Matrix)%in%Selection)]))
			
			nr.min.max.together=c(min.together,max.together)
			temp[[3]]=nr.min.max.together
			names(temp)[3]="nr.min.max.together"
			
			min.perc.together <- min.together/length(Selection) *100
			max.perc.together <- max.together/length(Selection) *100
			
			perc.min.max.together =c(min.perc.together,max.perc.together) 
			temp[[4]]=perc.min.max.together	
			names(temp)[4]="perc.min.max.together"
			
			temp[[5]]=list()
			names(temp)[5]="AllClusters"
			
			for(a in 1:length(clusternumbers)){
				temp[[5]][[a]]=list()
				names(temp[[5]])[a]=paste("Cluster",clusternumbers[a],sep=" ")
				
				temp[[5]][[a]][[1]]=clusternumbers[a]
				temp[[5]][[a]][[2]]=names(which(Matrix[i,]==clusternumbers[a])) #complete cluster
				temp[[5]][[a]][[3]]=intersect(Selection,temp[[5]][[a]][[2]]) #Objects from original selection in this cluster
				temp[[5]][[a]][[4]]=temp[[5]][[a]][[2]][which(!(temp[[5]][[a]][[2]] %in% Selection))] #Objects extra to this cluster
				names(temp[[5]][[a]])=c("clusternumber","Complete cluster","Objects from original selection in this cluster","Objects extra to this cluster")					
			}
			
			if(followMaxComps==TRUE){
				
				maxcluster=names(which(table(Matrix[i,which(colnames(Matrix)%in%Selection)])==max(table(Matrix[i,which(colnames(Matrix)%in%Selection)]))))
				temp[[6]]=maxcluster
				names(temp)[6]="Cluster with max Objects"
				complabels=rownames(Matrix[i,which(colnames(Matrix)%in%Selection & Matrix[i,]==as.numeric(maxcluster)),drop=FALSE])
				temp[[7]]=complabels
				names(temp)[7]="Complabels"
				complete.new.cluster=names(Matrix[i,which(Matrix[i,]==as.numeric(maxcluster))])
				temp[[8]]=complete.new.cluster
				names(temp)[8]="Complete.new.cluster"
				extra.new.cluster=complete.new.cluster[which(!(complete.new.cluster %in% Selection))]
				temp[[9]]=extra.new.cluster
				names(temp)[9]="Extra.new.cluster"
				
			}
			
			if(followClust==TRUE){
								
				temp[[6]]=cluster.interest
				names(temp)[6]="Cluster"
				complabels=rownames(Matrix[i,which(colnames(Matrix)%in%Selection & Matrix[i,]==as.numeric(cluster.interest)),drop=FALSE])
				temp[[7]]=complabels
				names(temp)[7]="Complabels"
				complete.new.cluster=names(Matrix[i,which(Matrix[i,]==as.numeric(cluster.interest))])
				temp[[8]]=complete.new.cluster
				names(temp)[8]="Complete.new.cluster"
				extra.new.cluster=complete.new.cluster[which(!(complete.new.cluster %in% Selection))]
				temp[[9]]=extra.new.cluster
				names(temp)[9]="Extra.new.cluster"
					
			}
						
			FoundClusters[[i]]=temp	
		}
		
		if(ClusterPlot==TRUE & !(is.null(ClusterNumber))){
			
			for(i in 1:dim(Matrix)[1]){
				temp=list()
				if(i==1){
					temp[[1]]=names(Matrix[i,which(Matrix[i,]==ClusterNumber)])
					names(temp)[1]=paste("Cluster ", ClusterNumber,sep="")
					PrevCluster=temp[[1]]
				}
				else{
					temp[[1]]=names(Matrix[i,which(Matrix[i,]==ClusterNumber)])
					names(temp)[1]=paste("Cluster ",  ClusterNumber,sep="")
					
					Diss=list()
					DissComps=NULL
					if(length((which(!(PrevCluster%in%temp[[1]]))!=0))){
						DissComps=PrevCluster[(which(!(PrevCluster%in%temp[[1]])))]
					}				
					if(!(is.null(DissComps))){ #Dissapeared comps: what is missing form PrevClust, where to did they move?
						cl=i
						for(t in 1:nrclusters){
							TempComps=which(DissComps%in%names(which(Matrix[cl,]==t)))
							disscl=NULL
							if(length(TempComps)!=0){
								disscl=list()
								disscl[[1]]=DissComps[TempComps]
								disscl[[2]]=t
								names(disscl)=c("Cpds","Cluster")
							}
							Diss[[t]]=disscl
						}
					}
					if(length(Diss)!=0){
						r=c()
						for(l in 1:length(Diss)){
							if(is.null(Diss[[l]])){
								r=c(r,l)
							}
						}
						if(!is.null(r)){
							Diss=Diss[-r]
						}
						for(l in 1:length(Diss)){
							names(Diss)[l]=paste("Comps Diss To Cluster " ,Diss[[l]][[2]],sep="")
						
						}
					}
					temp[[2]]=Diss
					names(temp)[2]="Dissapeared"
									
					Joined=list()
					JoinComps=NULL
					if(length((which(!(temp[[1]]%in%PrevCluster))!=0))){
						JoinComps=temp[[1]][(which(!(temp[[1]]%in%PrevCluster)))]
					}				
					if(!(is.null(JoinComps))){ #Dissapeared comps: what is missing form PrevClust, where to did they move?
						prevcl=i-1
						for(t in 1:nrclusters){
							TempComps=which(JoinComps%in%names(which(Matrix[prevcl,]==t)))
							joincl=NULL
							if(length(TempComps)!=0){
								joincl=list()
								joincl[[1]]=JoinComps[TempComps]
								joincl[[2]]=t
								names(joincl)=c("Cpds","Cluster")
								
							}
							Joined[[t]]=joincl
						}
					}
					if(length(Joined)!=0){
						r=c()
						for(l in 1:length(Joined)){
							if(is.null(Joined[[l]])){
								r=c(r,l)
							}
						}
						if(!is.null(r)){
							Joined=Joined[-r]
						}
						
						for(l in 1:length(Joined)){
							names(Joined)[l]=paste("Comps Joined From Cluster " ,Joined[[l]][[2]],sep="")
						
						}
					}
					temp[[3]]=Joined
					names(temp)[3]="Joined"	
					PrevCluster=temp[[1]]
				}
				FoundCl[[i]]=temp
			}
			
		}
		
		FoundClusters[[length(FoundClusters)+1]]=FoundCl

		FoundClusters[[length(FoundClusters)+1]]=Matrix
		
		return(FoundClusters)
	}	
	
	Found=ClusterDistribution.2(List,Selection,nrclusters,followMaxComps,followClust,fusionsLog,WeightClust,names=names)
	
	Matrix=Found[[length(Found)]]
	Found=Found[-length(Found)]
	
	FoundCl=Found[[length(Found)]]
	Found=Found[-length(Found)]
	
	if(is.null(names)){
		for(j in 1:dim(Matrix)[1]){
			names[j]=paste("Method",j,1)
		}
	}
	
	ClusterNumber=NULL
	if(class(Selection)=="numeric"){
		ClusterNumber=Selection
		Selection=colnames(Matrix)[which(Matrix[1,]==Selection)]

	}
	
	plottypein<-function(plottype,location){
		if(plottype=="pdf" & !(is.null(location))){
			grDevices::pdf(paste(location,".pdf",sep=""))
		}
		if(plottype=="new"){
			grDevices::dev.new()
		}
		if(plottype=="sweave"){
			
		}
	}
	plottypeout<-function(plottype){
		if(plottype=="pdf"){
			grDevices::dev.off()
		}
	}
	
	
	if(SelectionPlot==TRUE){
	
		
		if(followMaxComps==TRUE){
			lab1=c("Maximum of compounds of original cluster together")
			labelcluster=c()
			for(z in 1:length(Found)){
				labelcluster=c(labelcluster,Found[[z]]$Cluster)				
			}
		}
		else{lab1=c("Number of compounds still in original cluster")}
		
		nrcluster=c()
		nrcomps=c()
		for(j in 1:length(Found)){
			nrcluster=c(nrcluster,Found[[j]]$nr.clusters)
			nrcomps=c(nrcomps,length(Found[[j]]$Complabels))
		}
		
		if(is.null(ClusterNumber)){
			xl=c(0,length(Found)+0.5)
		}
		else{
			xl=c(1,length(Found)+0.5)
		}
		
		if(!(is.null(location))){
			location=paste(location,"_SelectionPlot.pdf",sep="")
			
		}
		plottypein(plottype,location)
		
		graphics::plot(type="n",x=0,y=0,xlim=xl,ylim=c(0,max(nrcluster,nrcomps)+2.0),xlab="",ylab="",xaxt="n",yaxt="n",cex.lab=1.25)
		graphics::lines(x=seq(1,length(Found)),y=nrcluster,lty=1,col="red",lwd=1.5)
		graphics::points(x=seq(1,length(Found)),y=nrcluster,pch=19,col="red",cex=1.5)
		
		if(is.null(ClusterNumber)){
			graphics::lines(x=c(0,1),y=c(length(Selection),nrcomps[1]),lty=1,col="black",lwd=1.5)
			graphics::points(x=0,y=length(Selection),pch=19,col="black",cex=1.5)
		}
		
		graphics::lines(x=seq(1,length(Found)),y=nrcomps,lty=1,col="blue",lwd=1.5)
		graphics::points(x=seq(1,length(Found)),y=nrcomps,pch=19,col="blue",cex=1.5)
		
		if(is.null(ClusterNumber)){
			graphics::text(0,length(Selection), "S",cex=1.5,pos=1,col="black",font=2)	
		}
		if(SelectionPlot==TRUE & followMaxComps==TRUE)
			graphics::text(seq(1,length(Found)),nrcomps, labelcluster,cex=1.5,pos=1,col="black",font=2)
		
		if(is.null(ClusterNumber)){
			graphics::axis(1,at=seq(0,length(Found)),labels=c("Selection",names),las=2,cex=1)
		}
		else{
			graphics::axis(1,at=seq(1,length(Found)),labels=c(names),las=2,cex=1)
		}
		graphics::axis(2,at=seq(0,max(nrcluster,nrcomps)+2.0),labels=seq(0,max(nrcluster,nrcomps)+2.0),cex=1)
		

		graphics::legend(legendposx,max(nrcluster,nrcomps)+legendposy,lty=c(1,1,0),pch=c(19,19,0),col=c("blue","red","black"),legend=c(lab1,"Number of clusters original cluster divided amongst","Cluster number"),bty="n",cex=1.2)
		
		plottypeout(plottype)
	}
	if(CompleteSelectionPlot==TRUE){
		if(!(is.null(location))){
			location=paste(location,"_CompleteSelectionPlot.pdf",sep="")
			
		}
		plottypein(plottype,location)
		nrcluster=c(1)
		nrcomps=c(length(Selection))
		for(j in 1:length(Found)){
			nrcluster=c(nrcluster,Found[[j]]$nr.clusters)
			nrcomps=c(nrcomps,length(Found[[j]]$Complabels))
		}
		
	
	
		graphics::plot(type="n",x=0,y=0,xlim=c(0,length(Found)+0.5),ylim=c(0,max(nrcluster,nrcomps)+1.5),xlab="",ylab="Number of Compounds",xaxt="n",yaxt="n")
		
		xnext=c()
		ynext=c()
		colorsp=c()
		
		if(is.null(ClusterNumber)){
			p=seq(0,length(Found))
		}
		else{
			p=seq(0,length(Found))
		}
		
		for(m in p){
			if(m==0){
				if(!(is.null(ClusterNumber))){
					howmany=length(Found[[1]]$AllClusters)
				}
				else{
					howmany=1
				}
			}
			else{
				howmany=length(Found[[m]]$AllClusters)
			}

			if(m==0){
				xnext=0.1
				ynext=c(length(Selection))
				if(!(is.null(ClusterNumber))){
					colorsp=cols[ClusterNumber]
				}
				else{
					colorsp=c("black")
				}
				
					
				graphics::points(x=xnext,y=ynext,col=colorsp,pch=19,cex=1.25)
				if(!(is.null(ClusterNumber))){
					labelcluster=ClusterNumber
				}
				else{
					labelcluster="S"
				}
				position=3
				if(!(is.integer(ynext[length(ynext)]))){
						position=1
				}
				graphics::text(xnext,ynext,labelcluster,cex=1.5,pos=position,col="black",font=2)					
				
				L1=list()
				for(n in 1:howmany){
					if(is.null(ClusterNumber)){				
						L1[[n]]=Selection
					}
					else{
						L1[[n]]=Found[[1]]$AllClusters[[n]][[3]]
					}
				}
				
			}
			else{
				xprev=xnext
				yprev=ynext
				ynext=c()
				colorsp=c()
				xnext=rep(seq(1,length(Found))[m],howmany)
				for(n in 1:howmany){										
					ynext=c(ynext,length(Found[[m]]$AllClusters[[n]][[3]]))
					
					colorsp=c(colorsp,cols[Found[[m]]$AllClusters[[n]]$clusternumber])
					
					if(length(ynext)>1){
						for(t in 1:(length(ynext)-1)){
							if(ynext[t]==ynext[length(ynext)]){
								ynext[length(ynext)]=ynext[length(ynext)]-0.3
							}
						}	
					}
					
					graphics::points(x=xnext[n],y=ynext[length(ynext)],col=colorsp[length(colorsp)],pch=19,cex=1.25)
					labelcluster=Found[[m]]$AllClusters[[n]]$clusternumber
					position=3
					if(!(is.integer(ynext[length(ynext)]))){
						position=1
					}
					graphics::text(xnext[n],ynext[length(ynext)],labelcluster,cex=1.5,pos=position,col="black",font=2)					
				}		
				
				
				L2=L1				
				L1=list()
				for(n in 1:howmany){
					L1[[n]]=Found[[m]]$AllClusters[[n]][[3]]					
				}
				
				for(q in 1:length(L1)){
					for(p in 1:length(L2)){
						if(length(which(L2[[p]] %in% L1[[q]])) != 0){
							graphics::segments(x0=xprev[p],y0=yprev[p],x1=xnext[q],y1=ynext[q],col=colorsp[q],lwd=2)
						}							
					}
				}
			}
		}
		if(is.null(ClusterNumber)){
			graphics::axis(1,at=seq(0,length(Found)),labels=c("Selection",names),las=2,cex=1.5)
		}
		else{
			graphics::axis(1,at=seq(0,length(Found)-1),labels=c(names),las=2,cex=1.5)
		}
		graphics::axis(2,at=seq(0,max(nrcluster,nrcomps)),labels=seq(0,max(nrcluster,nrcomps)),cex=1.5)
		graphics::legend(legendposx,max(nrcluster,nrcomps)+legendposy,pch=c(0),col=c("black"),legend=c("Cluster number"),bty="n",cex=1.2)
		
		
		plottypeout(plottype)
	}
	
	if(ClusterPlot==TRUE & !(is.null(ClusterNumber))){
		
		#FoundCl=Found[[length(Found)]]
		
		if(!(is.null(location))){
			location=paste(location,"_SelectionPlot.pdf",sep="")
		
		}
		plottypein(plottype,location)
		nrcluster=c(1)
		nrcomps=c(length(Selection))
		for(j in 1:length(FoundCl)){
			nrcomps=c(nrcomps,length(FoundCl[[j]][[1]]))
		}
		
		
		graphics::plot(type="n",x=0,y=0,xlim=c(0,length(FoundCl)-0.5),ylim=c(0,max(nrcluster,nrcomps)+1.5),xlab="",ylab="Number of Compounds",xaxt="n",yaxt="n")
		
		xnext=c()
		ynext=c()
		colorsp=c()
		
	
		p=seq(1,length(FoundCl))
		
		for(m in p){

			if(m==1){
				xnext=0.1
				ynext=length(FoundCl[[1]][[1]])
				colorsp=cols[ClusterNumber]
				graphics::points(x=xnext,y=ynext,col=colorsp,pch=19,cex=1.25)
				
				labelcluster=ClusterNumber
				
				position=3
				if(!(is.integer(ynext[length(ynext)]))){
					position=1
				}
				graphics::text(xnext,ynext,labelcluster,cex=1.5,pos=position,col="black",font=2)		
				
				
			}
			
			else{
				xprev=xnext
				yprev=ynext
				
				ynext=c()
				colorsp=c()

				xnext=m-1
				ynext=length(FoundCl[[m]][[1]])
				colorsp=cols[ClusterNumber]
				
				graphics::points(x=xnext,y=ynext,col=colorsp,pch=19,cex=1.25)
				graphics::segments(x0=xprev,y0=yprev,x1=xnext,y1=ynext,col=colorsp)
				labelcluster=ClusterNumber
				
				position=3
				if(!(is.integer(ynext[length(ynext)]))){
					position=1
				}
				graphics::text(xnext,ynext,labelcluster,cex=1.5,pos=position,col="black",font=2)
				
				#Dissapeared Compounds
				if(length(FoundCl[[m]]$Dissapeared)!=0){
					
					xdiss=rep(xnext,length(FoundCl[[m]]$Dissapeared))
					ydiss=c()
					colorsd=c()
					labs=c()
					for(d in 1:length(FoundCl[[m]]$Dissapeared)){
						ydiss=c(ydiss,length(FoundCl[[m]]$Dissapeared[[d]][[1]]))
						colorsd=c(colorsd,cols[FoundCl[[m]]$Dissapeared[[d]][[2]]])
						labs=c(labs,FoundCl[[m]]$Dissapeared[[d]][[2]])
					}
					
					if(length(ydiss)>1){
						for(t in 1:(length(ydiss)-1)){
							if(ydiss[t]==ydiss[length(ydiss)]){
								ydiss[length(ydiss)]=ydiss[length(ydiss)]-0.3
							}
						}	
					}
					
					graphics::points(x=xdiss,y=ydiss,col=colorsd,pch=19,cex=1.25)
					labelcluster=labs
					position=3
					if(!(is.integer(ydiss[length(ydiss)]))){
						position=1
					}
					graphics::text(xdiss,ydiss,labelcluster,cex=1.5,pos=position,col="black",font=2)	
					
					for(p in 1:length(ydiss)){
						graphics::segments(x0=xprev,y0=yprev,x1=xdiss[p],y1=ydiss[p],col=colorsd[p])
					}	
				
				}	
				
				
				if(length(FoundCl[[m]]$Joined)!=0){
					
					xjoin=rep(xprev,length(FoundCl[[m]]$Joined))
					yjoin=c()
					colorsj=c()
					labs=c()
					for(d in 1:length(FoundCl[[m]]$Joined)){
						yjoin=c(yjoin,length(FoundCl[[m]]$Joined[[d]][[1]]))
						colorsj=c(colorsj,cols[FoundCl[[m]]$Joined[[d]][[2]]])
						labs=c(labs,FoundCl[[m]]$Joined[[d]][[2]])
					}
					
					if(length(yjoin)>1){
						for(t in 1:(length(yjoin-1))){
							if(yjoin[t]==yjoin[length(yjoin)]){
								yjoin[length(yjoin)]=yjoin[length(yjoin)]-0.3
							}
						}	
					}
					
					graphics::points(x=xjoin,y=yjoin,col=colorsj,pch=19,cex=1.25)
					labelcluster=labs
					position=3
					if(!(is.integer(yjoin[length(yjoin)]))){
						position=1
					}
					graphics::text(xjoin,yjoin,labelcluster,cex=1.5,pos=position,col="black",font=2)	
					
					for(p in 1:length(yjoin)){
						graphics::segments(x0=xjoin[p],y0=yjoin[p],x1=xnext,y1=ynext,col=colorsj[p])
					}	
					
				}	
							
			}
					
		}
		graphics::axis(1,at=c(0.1,seq(1,length(FoundCl)-1)),labels=c(names),las=2,cex=1.5)
		graphics::axis(2,at=seq(0,max(nrcomps)),labels=seq(0,max(nrcomps)),cex=1.5)
		graphics::legend(legendposx,max(nrcluster,nrcomps)+legendposy,pch=c(0),col=c("black"),legend=c("Cluster number"),bty="n",cex=1.2)
		
		plottypeout(plottype)
		
	}	
	
	if(Table==TRUE & SelectionPlot==TRUE){
		SharedComps=Selection
		Extra=list()
		temp=c()
		for(a in 1:length(Found)){
			SharedComps=intersect(SharedComps,Found[[a]]$Complabels)
		}
		for(a in 1:length(Found)){
			Extra[[a]]=Found[[a]]$Complete.new.cluster[which(!(Found[[a]]$Complete.new.cluster%in%SharedComps))]
			names(Extra)[a]=names[a]
			if(followMaxComps==TRUE){
				names(Extra)[a]=paste(names(Extra)[a],"_",labelcluster[a],sep="")
			}
			temp=c(temp,length(Extra[[a]]))
		}
		
		ExtraOr=Selection[which(!(Selection%in%SharedComps))]
		
		collength=max(length(SharedComps),length(ExtraOr),temp)
		
		if(length(SharedComps)<collength){
			spaces=collength-length(SharedComps)
			SharedComps=c(SharedComps,rep(" ",spaces))
		}
		
		if(length(ExtraOr)<collength){
			spaces=collength-length(ExtraOr)
			ExtraOr=c(ExtraOr,rep(" ",spaces))
		}
		
		for(b in 1:length(Extra)){
			if(length(Extra[[b]])<collength){
				spaces=collength-length(Extra[[b]])
				Extra[[b]]=c(Extra[[b]],rep(" ",spaces))
			}
			
		}		
		table=data.frame(ExtraOr=ExtraOr,SharedComps=SharedComps)
		for(b in 2:length(Extra)){
			table[1+b]=Extra[[b]]
			colnames(table)[1+b]=names(Extra)[b]
		}
		
	}
	else{table=list()}
	
	names(Found)=names
	Found[[length(Found)+1]]=table
	names(Found)[length(Found)]="Table"
	return(Found)	
}
