gseaplot <-
function(x,cat=NULL){
	if(is.null(cat)) stop("BOG : Category needs to be specifiled.")

    stat=x[[1]]
	if(is.null(stat$gsea)) stop("BOG : gsea was not set to be TRUE in creating BOG object.")
	cat=toupper(cat)
    par(mfrow=c(1,1))
    map=stat$plot$map
    delta=stat$plot$delta
    S=stat$plot$S_location
    index=values(map,keys=cat,USE.NAMES=FALSE)
    delta=stat$plot$delta[,index]
    pval=stat$gsea$adj[which(stat$gsea$COG==cat)]
    n.index=length(delta)
    max.y=max(delta)
    max.x=which(delta==max.y)
		
    horiz.y=rep(max.y,n.index)
    horiz.x=c(1:n.index)
		
    ver.y=seq(0,max.y+0.2,length.out=1000)
    ver.x=rep(max.x,1000)
		
    if(pval<0.00001){
        label=paste("Ranked genes - COG ",cat," genes in red"," : p-value < 0.00001",sep="")
		}else{
			label=paste("Ranked genes - COG ",cat," genes in red : p-value = ",sep="")
			label=paste(label,round(pval,digits=5))
		}

		plot(delta,ylim=c(0,1),xlab=label,ylab="Running GSEA Score",type="l",col="brown",cex.lab=1.4, cex.axis=1.5, cex.main=1.5, cex.sub=1.5)
		points(horiz.x,horiz.y,type="l",col="red",lty=3)
		points(ver.x,ver.y,type="l",col="red",lty=3)
		points(max.x,max.y,col="red",pch=18)
		axis(4,max.y,labels=round(max.y,digits=3),lty=3,las=0,cex.axis=1.5)
		axis(3,max.x,labels=max.x,lty=3,las=0,cex.axis=1.5)
		
		Spoints.x=which(S[,index]!=1)
		Spoints.ys=rep(0,length(Spoints.x))
		Spoints.ye=rep(0.1,length(Spoints.x))
		segments(Spoints.x,Spoints.ys,Spoints.x,Spoints.ye,col="green")
		
		Spoints.x=which(S[,index]==1)
		Spoints.ys=rep(0,length(Spoints.x))
		Spoints.ye=rep(0.1,length(Spoints.x))
		segments(Spoints.x,Spoints.ys,Spoints.x,Spoints.ye,col="red")	
}
