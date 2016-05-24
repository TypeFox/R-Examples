#--------------------------

#' profile plot for bic vs outside bic either by gene or condtion.
#' plots either for a single or for all in one data frame.
#' profile plot of the biclust genes and the remain genes 
#' @ dset a parameter for the dataset
#' @ bres a parameter of the bicluster object
#' @ plot.bic a vector of different plot names.
#' @teta; @ph; values -360<teta/ph< 360, rotates the 3D plot
#' @mname; method name;'biclust','fabia' or 'isa2'
#' @bnum; biclust number; must be an existed biclust. 

#-------------------------------------#
profileBic <- function(dset,bres,mname=c("fabia","isa2","biclust","bicare"),bplot="all",gby="genes",bnum=1,teta=120,ph=30,fabia.thresZ=0.5,fabia.thresL=NULL,BClabel=TRUE,gene.lines=NULL,condition.lines=NULL){
	
	# Small extra for the GUI
	if(bplot=="threeD"){bplot <- "3D"}	
	
	highlight.g <- gene.lines
	highlight.c <- condition.lines
	


	#check if the bic number to plotted is specified and
	if(any(!bplot %in% c("all","boxplot","lines","3D","histogram"))) {
		stop("`bplot' must be one of `boxplot', `lines', `3D','histogram' or `all'")
	}
	if(any(!mname %in% c("fabia","isa2","biclust","bicare"))){
		stop("`mname' must be one of `fabia',`isa2', 'biclust' or 'bicare'")
	}
	if(any(!gby %in% c("genes","conditions"))) {
		stop("`gby' must be one of `genes', `conditions'")
	}
	if(any(!mname %in% c("fabia","isa2","biclust","bicare"))){
		stop("`mname' must be one of `fabia',`isa2', 'biclust' or 'bicare'")
	} 
	#par(mfrow=c(1,1))
	ind.gc<-indexedBic(dset,bres,mname,bnum,fabia.thresZ=fabia.thresZ,fabia.thresL=fabia.thresL)# returns the required indecies based on thier method names
	indg<-ind.gc[[1]]
	indc<-ind.gc[[2]]
	
		
	# gby conditions	
	if(gby=="conditions"){
		#group the genes in to two.
		cnams <- colnames(dset)
		gnams <- rownames(dset)
		grp <- rep(1, length(cnams))
		grp[indc] <- 2
		d<-dset[indg, order(grp, decreasing=TRUE)]
		dbc<-dset[indg,(grp==2)]
		
		if(bplot=="lines"){
			
			
			if(BClabel){
				mar.temp <- par()$mar
				
				max.g.nchar <- max(sapply(gnams,FUN=nchar))
				if(max.g.nchar>11){
					par(mar=par()$mar+c(0,0.2*(max.g.nchar-11),0,0))
				}
				max.c.nchar <- max(sapply(cnams,FUN=nchar))
				if(max.c.nchar>14){
					par(mar=par()$mar+c(0.15*(max.c.nchar-11),0,0,0))
					
				}
			}
			
			col.line <- rep("green3",nrow(d))
			col.bcline <- rep("red",nrow(d))
			lwd.line <- rep(1,nrow(d))
			
			if(!is.null(highlight.g)){
				
				if(class(highlight.g)=="character"){
					sel.g <- sapply(highlight.g,FUN=function(x){
								temp <- which(x==gnams[indg])
								if(length(temp)==0){stop(paste0("gene.lines contains genes which are not in bicluster ",bnum),call.=FALSE)}
								return(temp)
							})
					
				}
				if(class(highlight.g)=="numeric"){
					sel.g <- sapply(highlight.g,FUN=function(x){
								temp <- which(x==indg)
								if(length(temp)==0){stop(paste0("gene.lines contains genes which are not in bicluster ",bnum),call.=FALSE)}
								return(temp)
							})
				}
				
				
				col.line[-sel.g] <- "grey"
				col.bcline[-sel.g] <- "grey"
				lwd.line[sel.g] <- 2
			}
			if(!is.null(highlight.c)){
				
				if(class(highlight.c)=="character"){
					
					sel.c <- sapply(highlight.c,FUN=function(x){
								temp <- which(x==cnams[indc])
								if(length(temp)==0){stop(paste0("condition.lines contains conditions which are not in bicluster ",bnum),call.=FALSE)}
								return(temp)
							})
				}
				if(class(highlight.c)=="numeric"){
					
					sel.c <- sapply(highlight.c,FUN=function(x){
								temp <- which(x==indc) 
								if(length(temp)==0){stop(paste0("condition.lines contains conditions which are not in bicluster ",bnum),call.=FALSE)}
								return(temp)
							})
				}
				

			}
			
			
			matplot(y =t(d),type ="n",col="green3", xlab="Condtions",ylab="Expression", axes=T, pch = rep(1, ncol(d)))
			
			if(is.null(highlight.g)){
				matlines(y = t(d), type = "l",lty = rep(1, nrow(d)) ,col="green3", lwd = 1, pch = 1)
				matlines(y = t(dbc), type = "l",lty = rep(1, nrow(d)) ,col="red", lwd = 1, pch = 1)
			}
			else{
				matlines(y = t(d)[,c((1:length(indg))[-sel.g],sel.g)], type = "l",lty = rep(1, nrow(d)) ,col=col.line[c((1:length(indg))[-sel.g],sel.g)], lwd = lwd.line[c((1:length(indg))[-sel.g],sel.g)], pch = 1)
				matlines(y = t(dbc)[,c((1:length(indg))[-sel.g],sel.g)], type = "l",lty = rep(1, nrow(d)) ,col=col.bcline[c((1:length(indg))[-sel.g],sel.g)], lwd = lwd.line[c((1:length(indg))[-sel.g],sel.g)], pch = 1)
				
			}
			
			
			legend("topright",c( "Biclust conditions","Outside Condtions"), col=c("red","green3"),lty=c(1,1))
			
			if(BClabel){
				if(is.null(highlight.c)){
					axis(side=1,at=c(1:length(indc)),labels=colnames(dbc),las=2,cex.axis=0.7)
				}
				else{
					axis(side=1,at=c(1:length(indc))[-sel.c],labels=colnames(dbc)[-sel.c],las=2,cex.axis=0.7,col.axis="grey")
					axis(side=1,at=c(1:length(indc))[sel.c],labels=colnames(dbc)[sel.c],las=2,cex.axis=0.7,col.axis="black")
				}
				if(is.null(highlight.g)){
					axis(side=2,at=d[,1],labels=gnams[indg],las=1,cex.axis=0.7)
				}
				else{
					axis(side=2,at=d[,1][-sel.g],labels=gnams[indg][-sel.g],las=1,cex.axis=0.7,col.axis="grey")
					axis(side=2,at=d[,1][sel.g],labels=gnams[indg][sel.g],las=1,cex.axis=0.7,col.axis="black")
				}
				
#				print(par()$mar)
				par(mar=mar.temp)
#				print(par()$mar)
			}
#			box()
		}
		
		if(bplot=="boxplot"){
			boxplot.matrix(d,col="green3",main="",axes=T,lty=1)
			boxplot.matrix( dbc,add=T,col="red",axes=F,lty=1)
			legend("topright",c( "Bic conditions","Outbic Condtions"), col=c("red","green3"),lty=c(1,1))
			box()
			
		}
			
		if(bplot=="histogram"){
			#hist 
			hist(d,col="green3",xlab="",main="",lty=1)
			hist( dbc,col="red",add=T,xlab="",main="",lty=1)
			legend("topright",c( "Bic conditions","Outbic Condtions"), col=c("red","green3"),lty=c(1,1))
			box()
			
		}
			
		if(bplot=="3D"){
			d1<-c(1:nrow(d))
			d2<-c(1:ncol(d))
			
			fill <- matrix("green3", nrow = nrow(d), ncol = ncol(d))
			fill[,(grp==2)] <- "red";fill<-sort(fill, decreasing=T)
			persp(d1,d2,d,theta = teta, phi = ph, expand = 0.5, col = fill,
				 ltheta = 120, shade = 0.75, ticktype = "detailed",
				,xlab="Genes",ylab="Condtions",zlab="Gene expression")
			legend("topright",c( "Bic conditions","Outbic Condtions"), col=c("red","green3"),lty=c(1,1))
			box()

		}
		if(bplot=="all"){
			#par(mfrow=c(2,2))
			.checkcurrentgrid(4,2,2)
			profileAll(dset,indg,indc,grp,gby)
		}
	}
	
	#gby genes
	if(gby=="genes"){
		#group the genes in to two.
		gnams <- rownames(dset)
		cnams <- colnames(dset)
		grp <- rep(1, length(gnams))
		grp[indg] <- 2
		d<-dset[order(grp, decreasing=T),indc]
		dbc<-dset[(grp==2),indc]
		
		if(bplot=="lines"){
			
			if(BClabel){
				mar.temp <- par()$mar
				
				max.g.nchar <- max(sapply(gnams,FUN=nchar))
				if(max.g.nchar>14){
					par(mar=par()$mar+c(0.15*(max.g.nchar-11),0,0,0))
				}
				max.c.nchar <- max(sapply(cnams,FUN=nchar))
				if(max.c.nchar>11){
					par(mar=par()$mar+c(0,0.2*(max.c.nchar-11),0,0))
					
				}
			}	
			
			col.line <- rep("black",nrow(dset))
			col.bcline <- rep("red",nrow(dset))
			lwd.line <- rep(1,nrow(dset))
			
			if(!is.null(highlight.c)){
				if(class(highlight.c)=="character"){
					
					sel.c <- sapply(highlight.c,FUN=function(x){
								temp <- which(x==cnams[indc])
								if(length(temp)==0){stop(paste0("condition.lines contains conditions which are not in bicluster ",bnum),call.=FALSE)}
								return(temp)
							})
				}
				if(class(highlight.c)=="numeric"){
					
					sel.c <- sapply(highlight.c,FUN=function(x){
								temp <- which(x==indc) 
								if(length(temp)==0){stop(paste0("condition.lines contains conditions which are not in bicluster ",bnum),call.=FALSE)}
								return(temp)
							})
				}
				

				col.line[-sel.c] <- "grey"
				col.bcline[-sel.c] <- "grey"
				lwd.line[sel.c] <- 2
			}
			if(!is.null(highlight.g)){
				if(class(highlight.g)=="character"){
					sel.g <- sapply(highlight.g,FUN=function(x){
								temp <- which(x==gnams[indg])
								if(length(temp)==0){stop(paste0("gene.lines contains genes which are not in bicluster ",bnum),call.=FALSE)}
								return(temp)
							})
					
				}
				if(class(highlight.g)=="numeric"){
					sel.g <- sapply(highlight.g,FUN=function(x){
								temp <- which(x==indg)
								if(length(temp)==0){stop(paste0("gene.lines contains genes which are not in bicluster ",bnum),call.=FALSE)}
								return(temp)
							})
				}		
			}
			
			
			matplot(y = d, type = "n",xlab="Genes",ylab="Expression", axes=T, pch = rep(1, ncol(dset)))
			
			if(is.null(highlight.c)){
				matlines(y = d, type = "l", col="black",lty = rep(1, nrow(dset)),lwd = 1, pch = 1) 
				matlines(y = dbc, type = "l",col="red",lty = rep(1, nrow(dset)) , lwd = 1, pch = 1)
			}
			else{
				matlines(y = d[,c((1:length(indc))[-sel.c],sel.c)], type = "l", col=col.line[c((1:length(indc))[-sel.c],sel.c)],lty = rep(1, nrow(dset)),lwd = lwd.line[c((1:length(indc))[-sel.c],sel.c)], pch = 1) 
				matlines(y = dbc[,c((1:length(indc))[-sel.c],sel.c)], type = "l",col=col.bcline[c((1:length(indc))[-sel.c],sel.c)],lty = rep(1, nrow(dset)) , lwd = lwd.line[c((1:length(indc))[-sel.c],sel.c)], pch = 1)
			}
			
			
			legend("topright",c( "Bic genes","Outbic genes"), col=c("red","black"),lty=c(1,1))
			if(BClabel){
				
				if(is.null(highlight.c)){
					axis(side=2,at=d[1,],labels=cnams[indc],las=1,cex.axis=0.7)
					
				}
				else{
					axis(side=2,at=d[1,][-sel.c],labels=cnams[indc][-sel.c],las=2,cex.axis=0.7,col.axis="grey")
					axis(side=2,at=d[1,][sel.c],labels=cnams[indc][sel.c],las=2,cex.axis=0.7,col.axis="black")
				}
				
				if(is.null(highlight.g)){
					axis(side=1,at=c(1:length(indg)),labels=rownames(dbc),las=2,cex.axis=0.7)
					
				}
				else{
					axis(side=1,at=c(1:length(indg))[-sel.g],labels=rownames(dbc)[-sel.g],las=1,cex.axis=0.7,col.axis="grey")
					axis(side=1,at=c(1:length(indg))[sel.g],labels=rownames(dbc)[sel.g],las=1,cex.axis=0.7,col.axis="black")
				}
				
#			print(par()$mar)
				par(mar=mar.temp)
#			print(par()$mar)
			}
#		box()
		}
		
		if(bplot=="boxplot"){
		#boxplot
		boxplot.matrix(t(d),col="black",main="",xlab="",axes=T)
		lines(t(dbc),col="red",xlab="")
		legend("topright",c( "Bic genes","Outbic genes"), col=c("red","black"),lty=c(1,1))
		box()
		}
			
		if(bplot=="histogram"){
		#histogram
		hist(d,col="black",xlab="",main="")
		hist(dbc,col="red",add=T,xlab="",main="")
		legend("topright",c( "Bic genes","Outbic genes"), col=c("red","black"),lty=c(1,1))
		box()
		}
			
		if(bplot=="3D"){
			d1<-c(1:nrow(d))
			d2<-c(1:ncol(d))
			fill <- matrix("black", nrow = nrow(d), ncol = ncol(d))
			fill[(grp==2),] <- "red";fill<-sort(fill, decreasing=T)
			persp(d2,d1,t(d),theta = teta, phi = ph, expand = 0.5, col = fill,
				 ltheta = 120, shade = 0.75, ticktype = "detailed",
				,xlab="Condtions",ylab="Genes",zlab="Gene expression")
			legend("topright",c( "Bic genes","Outbic genes"), col=c("red","black"),lty=c(1,1))
			box()
		}
		if(bplot=="all"){
			#par(mfrow=c(2,2))
			.checkcurrentgrid(4,2,2)
			profileAll(dset,indg,indc,grp,gby)
		}
	}
}
