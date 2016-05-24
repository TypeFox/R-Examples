#prettyScree <- function(eigs,retain.col="mediumorchid4",dismiss.col="gray",n.comps=NULL,perc.exp=1.0,comp.labels=NULL,show.only.good=FALSE,full.lim=FALSE,main=""){
prettyScree <- function(eigs,retain.col="mediumorchid4",dismiss.col="gray",perc.exp=1.0,n.comps=NULL,broken.stick=TRUE,kaiser=TRUE,main=""){

	##stolen from http://lamages.blogspot.com/2013/04/how-to-change-alpha-value-of-colours-in.html
	add.alpha <- function(col, alpha=1){

	  apply(sapply(col, col2rgb)/255, 2, 
                     function(x) 
                       rgb(x[1], x[2], x[3], alpha=alpha))  
	}
	
	
	eig.length <- length(eigs)
	mean.eig <- mean(eigs)
	eigs.round <- round(eigs,digits=3)
	exp.var <- eigs/sum(eigs)*100
	exp.var.round <- round(exp.var,digits=2)

	# if(is.null(comp.labels)){
		# #comp.labels <- paste(exp.var
	# }
	if(is.null(n.comps)){
		n.comps <- eig.length
	}
	if(n.comps > eig.length || n.comps < 1){
		n.comps <- eig.length
	}
	
	
	perc.exp.comps <- cumsum(exp.var) < (perc.exp * 100)
	perc.exp.comps[head(which(!(perc.exp.comps)),n=1)] <- TRUE
	
	keep.n.comps <- rep(FALSE,eig.length)
	keep.n.comps[1:n.comps] <- rep(TRUE,length(1:n.comps))
	
	comps.tests <- rbind(perc.exp.comps,keep.n.comps)
	
	if(broken.stick){
		broken.stick.distribution <- unlist(lapply(X=1:eig.length,FUN=function(x,n){return(tail((cumsum(1/x:n))/n,n=1))},n=eig.length)) * 100
		broken.stick.comps <- exp.var > broken.stick.distribution
		broken.stick.comps[head(which(!broken.stick.comps),n=1):eig.length] <- rep(FALSE,length(head(which(!broken.stick.comps),n=1):eig.length))
		comps.tests <- rbind(comps.tests,broken.stick.comps)	
	}
	if(kaiser){
		kaiser.mean.comps <- eigs > mean.eig
		comps.tests <- rbind(comps.tests,kaiser.mean.comps)			
	}
	
	comp.sums <- colSums(comps.tests)
	alpha.map <- 1/abs((comp.sums-(nrow(comps.tests)+1)))
	color.map <- rep(retain.col,eig.length)
	for(i in 1:eig.length){color.map[i] <- add.alpha(color.map[i],alpha.map[i])}
	color.map[which(comp.sums==0)] <- rep(dismiss.col,sum(comp.sums==0))
	
	
	dev.new()
	par(mar=c(5, 5, 4, 5) + 0.1)
	these.sizes <- ((log(exp.var) + abs(min(log(exp.var))))/max(log(exp.var) + abs(min(log(exp.var))))+0.1) * 3
	plot(exp.var,axes=FALSE,ylab="",xlab="Components",type="l",main="",ylim=c(-1,max(exp.var)))
	points(exp.var,cex= these.sizes,pch=20,col=dismiss.col)
	points(exp.var,cex= these.sizes,pch=21,bg=color.map)		
	box()
	axis(2,at= exp.var,labels=eigs.round,las=2,lwd=3,cex.axis=.65)
	mtext("Eigenvalues",2,line=4)
	axis(4,at= exp.var,labels=paste(exp.var.round,"%",sep=""),las=2,lwd=3,cex.axis=.65)
	mtext("Explained Variance",4,line=4)
	axis(1,at=1:length(exp.var),lwd=3)	

	return(comps.tests)

}