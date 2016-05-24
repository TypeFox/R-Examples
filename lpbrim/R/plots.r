#' @title Plot the modules as a network
#' @description Plot the modules as a network
#' @export
#' @param mod The output of findModules
plotModules <- function(mod)
{
	opar <- graphics::par(no.readonly=TRUE)
	##
	#TX <- spread(x[x>0],0.5,3)
 	x <- mod$M
	x[x>0] <- 1
	##
 	graphics::par(mar=c(0,0,0,0),xaxt='n',yaxt='n')
	##
 	S <- mod$S
 	Modules <- numeric(nrow(S))
 	names(Modules) <- rownames(S)
 	for(i in 1:nrow(S))
 	{
 		Modules[i] <- as.numeric(colnames(S)[S[i,]==1])
 	}
 	ModNum <- as.numeric(as.factor(Modules))
 	names(ModNum) <- names(Modules)
 	Tm <- ModNum[rownames(x)]
 	Bm <- ModNum[colnames(x)]
 	## Initate Plot
 	CommColor <- grDevices::colorRampPalette(RColorBrewer::brewer.pal(11,'Spectral'))(ncol(S))
 	graphics::plot(0,xlim=c(-1,2),ylim=c(0,1),pch=NA,bty='n')
 	##
 	SeqT <- rank(Tm,ties.method='random')
 	SeqB <- rank(Bm,ties.method='random')
 	## Within each module
 	for(cmod in unique(ModNum))
 	{
 		# degree-ranking of the HTL species
 		Tcurr <- SeqT[Tm==cmod]
   		if(length(Tcurr)>1)
   		{
   			Gen <- plyr::aaply(x[names(Tcurr),],1,function(x) sum(x>0))
	   		WCrank <- rank(Gen,ties.method='random')
	   		Tcurr <- (min(Tcurr)-1)+WCrank
	   		SeqT[Tm==cmod] <- Tcurr
   		}

   		# degree-ranking of the HTL species
 		Bcurr <- SeqB[Bm==cmod]
   		if(length(Bcurr)>1)
   		{
	 		Gen <- plyr::aaply(x[,names(Bcurr)],2,function(x) sum(x>0))
	   		WCrank <- rank(Gen,ties.method='random')
	   		Bcurr <- (min(Bcurr)-1)+WCrank
	   		SeqB[Bm==cmod] <- Bcurr
   		}
 	}

 	SeqT <- spread(SeqT)
 	SeqB <- spread(SeqB)

 	for(ts in 1:nrow(x)) for(bs in 1:ncol(x)) if(x[ts,bs]>0)
 	{
 		ccol <- ifelse(Tm[ts]==Bm[bs],CommColor[Tm[ts]],'grey')
 		clwd <- ifelse(Tm[ts]==Bm[bs],2,1)
 		graphics::segments(0,SeqB[bs],1,SeqT[ts],lwd=clwd,col=ccol)
 	}

 	## Plot Top Species
 	graphics::points(rep(1,nrow(x)),SeqT,col=CommColor[Tm],pch=19,cex=1.4)
 	graphics::text(rep(1.1,nrow(x)),SeqT,rownames(x),cex=0.8,adj=c(0,0.5))
 	## Plot Bottom Species
 	graphics::points(rep(0,ncol(x)),SeqB,col=CommColor[Bm],pch=19,cex=1.4)
	graphics::text(rep(-0.1,ncol(x)),SeqB,colnames(x),cex=0.8,adj=c(1,0.5))
 	##
 	graphics::par(opar)
}

#' @title Plot the modules as a matrix
#' @description Plot the modules as a network
#' @export
#' @param mod The output of findModules
#' @param mode The type of plot, either blocks (different modules are colored), frames (modules are outlined), or both
plotMatrixModules <- function (mod,mode='blocks') {

	if(!(mode%in%c('blocks','frames','both'))) warning('Plot mode should be one of blocks, frames, or both')

	x <- mod$M
 	S <- mod$S
 	Modules <- numeric(nrow(S))
 	names(Modules) <- rownames(S)
 	for(i in 1:nrow(S))
 	{
 		Modules[i] <- as.numeric(colnames(S)[S[i,]==1])
 	}
 	ModNum <- as.numeric(as.factor(Modules))
 	names(ModNum) <- names(Modules)
 	Tm <- ModNum[rownames(x)]
 	Bm <- ModNum[colnames(x)]

	opar <- graphics::par(no.readonly=TRUE)
 	graphics::par(xaxt='n',yaxt='n',mar=c(0,0,0,0),mai=c(0.2,0.2,0.2,0.2))

 	CommColor <- grDevices::colorRampPalette(RColorBrewer::brewer.pal(11,'Paired'))(ncol(S))

 	SeqT <- rank(Tm,ties.method='random')
 	SeqB <- rank(Bm,ties.method='random')

 	## Within each module
 	for(cmod in unique(ModNum))
 	{
 		## degree-ranking of the HTL species
 		Tcurr <- SeqT[Tm==cmod]
   		if(length(Tcurr)>1)
   		{
   			Gen <- plyr::aaply(x[names(Tcurr),],1,function(x) sum(x>0))
	   		WCrank <- rank(Gen,ties.method='random')
	   		Tcurr <- (min(Tcurr)-1)+WCrank
	   		SeqT[Tm==cmod] <- Tcurr
   		}

   		## degree-ranking of the HTL species
 		Bcurr <- SeqB[Bm==cmod]
   		if(length(Bcurr)>1)
   		{
	 		Gen <- plyr::aaply(x[,names(Bcurr)],2,function(x) sum(x>0))
	   		WCrank <- rank(Gen,ties.method='random')
	   		Bcurr <- (min(Bcurr)-1)+WCrank
	   		SeqB[Bm==cmod] <- Bcurr
   		}
 	}

 	Tm <- Tm[names(SeqT)]
 	Bm <- Bm[names(SeqB)]

graphics::plot(0,pch=NA,xlim=c(0.9,ncol(x)+0.1),ylim=c(0.9,nrow(x)+0.1),asp=1,xlab='',ylab='',bty='n')
	graphics::rect(0.2,0.2,ncol(x)+0.8,nrow(x)+0.8)
   	## FRAMES
   	if(mode%in%c('frames','both'))
   	{
   		Ucom <- unique(Tm)
   		for(cc in 1:length(Ucom))
   		{
   			Ccom <- Ucom[cc]
   			Tpos <- range(SeqT[Tm==Ccom])
	  		Bpos <- range(SeqB[Bm==Ccom])
	  		Tpos <- Tpos + c(-0.45,0.45)
	  		Bpos <- Bpos + c(-0.45,0.45)
	  		## Module
	  		graphics::rect(Bpos[1],Tpos[1],Bpos[2],Tpos[2],lwd=2,border=ifelse(mode=='both',CommColor[Ccom],'black'))
   		}
   	}
   	## Squares
	for(i in 1:length(SeqT))
	{
	  for(j in 1:length(SeqB))
	  {
	  	Tcoord <- SeqT[i]
		Bcoord <- SeqB[j]
	  	ccol <- 'grey'
	  	if(mode!='frames')
	  	{
	  		ccol <- ifelse(Tm[i]==Bm[j],CommColor[Tm[i]],'grey')
	  	} else {
	  		ccol <- ifelse(Tm[i]==Bm[j],'black','grey')
	  	}
	  	if(x[i,j]>0) graphics::symbols(Bcoord,Tcoord,squares=0.7,add=TRUE,inches=FALSE,bg=ccol,fg=NA)
	  }
	}
	graphics::par(opar)
}
