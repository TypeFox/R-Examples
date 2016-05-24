
chrom.isoreg <- function(A, tresh=0.2)
{

	A <- as.matrix(A)
	if(all(is.na(A))) return(A)
	if(ncol(A)==0) return(A)
	if(sum(A)==0) return(A)
	#max.ind <- apply(abs(A),2,which.max)
	#max.sign <- sweep(A,2,max.ind)
	#mul.sign <- diag(as.matrix(apply(as.matrix(A[max.ind,]),2,sign)))
	#A <- sweep(as.matrix(A),2,mul.sign,"*")
	#x <- A[,2]
	
	A <- apply(as.matrix(A),2,function(x)
	{			
		if(any(diff(x)>0))
		{
		if(which.max(x)==1) {x[1:which(diff(x)>0)[1]] <- 0}
		xr <- x[length(x):1]
		if(which.max(x)==1)
		{
			xr[1:which(diff(xr)>0)[1]] <- 0
			x <- xr[length(x):1]
		}
		}	
		x			
	})	
	
		
	Ar <- apply(as.matrix(A),2,function(x)
	{		
		x[x<0] <- 0
		max.x <- max(x, na.rm=T)
		x <- normalize(x)
		x.rigth <- x[which.max(x):1]
		if(which.max(x)<length(x)) {x.left <- x[(which.max(x)+1):length(x)]}else{x.left <- NA}
		
		if(length(x.rigth)>2)
		{
			flag.isoreg <- F
			for(i in 2:length(x.rigth))
			{
				if(x.rigth[i]<tresh) flag.isoreg=T
				if(flag.isoreg) 
					if(x.rigth[i]>x.rigth[i-1]) x.rigth[i] <- x.rigth[i-1]
			}
		}
		if(length(x.left)>2)
		{
			flag.isoreg <- F
			for(i in 2:length(x.left))
			{
				if(x.left[i]<tresh) flag.isoreg=T
				if(flag.isoreg) 
					if(x.left[i]>x.left[i-1]) x.left[i] <- x.left[i-1]
			}
		}
		x <- c(x.rigth[length(x.rigth):1],x.left)
		del.na <- which(is.na(x)==T)
		if(length(del.na)!=0) x <- x[-del.na]
		#length(x)
		#plot(x, type="l")
		x <- x*max.x
		x
	})
		
	xmax.vect <- apply(Ar,2,max)
	Ar.diff <- normalize(apply(Ar,2,diff))
	Ar.n <- apply(Ar.diff,2, function(x){	
		x[which(abs(x[1:which.max(x)])<0.01)] <- 0
		x[(which.min(x)-1)+ which(abs(x[which.min(x):length(x)])<0.01)] <- 0
		x <- cumsum(c(0,x))
		x[x<0] <- 0
		x[normalize(x)<0.01] <- 0
		x
	})
	
	Ar.n <- sweep(normalize(Ar.n),2,xmax.vect,"*")
	Ar.n	
}


getS.OSD <- function(mod.c, D, ref.response=NULL, beta=2,cutoff=0.05)
{	
		#beta <- 2
		thr.det <- cutoff
		an.win <- D[which(normalize(mod.c)>thr.det),]
		an.win[an.win<0] <- 0
		if(dim(as.matrix(an.win))[1]<3 | dim(as.matrix(an.win))[2]<3) return(rep(0,ncol(D)))
		mod.l <- mod.c[which(normalize(mod.c)>thr.det)]	
		comp.m <- an.win^(1/beta)
		if(any(is.na(comp.m))) comp.m[is.na(comp.m)] <- 0
		pr.an <- try(prcomp(comp.m), silent=T)
		if(class(pr.an)=="try-error") {return(rep(0,ncol(D)))}
		
		if(is.null(ref.response))
		{
			mat.cor <- cor(mod.l, pr.an$x)
			cor.vect <- abs(mat.cor)
			cor.vect[summary(pr.an)$importance[2,]<0.005] <- 0
			i.comp <- order(cor.vect, decreasing=T)
		}else{
			mat.cor <- cor(ref.response, pr.an$rotation)
			cor.vect <- abs(mat.cor)
			cor.vect[summary(pr.an)$importance[2,]<0.005] <- 0
			i.comp <- order(cor.vect, decreasing=T)	
			}	
			
		spc <- pr.an$rotation[,i.comp[1]]*sign(mat.cor)[i.comp[1]]  #*pr.an$sdev[1]
	
		spc[spc<0] <- 0 
		spc <- spc^beta
		spc <- refine.extraction(D, mod.c, spc)
		spc <- normalize(spc)#*apex.h		
		spc
}

refine.extraction <- function(m.data, c.model, s.model)
{
	
	neg.sign.flag <- any(sign(normalize(m.data[which.max(c.model),])-normalize(s.model))==(-1))
	its <- 1
	while(neg.sign.flag)
	{
		substr <- normalize(m.data[which.max(c.model),]) - normalize(s.model)
		substr[substr>0] <- 0
		s.model <- normalize(s.model) + substr
		neg.sign.flag <- any(sign(normalize(m.data[which.max(c.model),])-normalize(s.model))==(-1))
		its <- its + 1
		if(its>10) break
	}
	s.model[s.model<0] <- 0
		
	s.model	
}

getC.tP <- function(winD, target.s)
{
	winD[,which(normalize(target.s)<0.01)] <- 0
	winDU <- winD^0.5
	target.sU <- target.s^0.5
	winDU[is.na(winDU)] <- 0
	target.sU[is.na(target.sU)] <- 0
	
	Cps <- as.vector(winDU %*% target.sU)^2
	Cps <- normalize(Cps)
	target.s <- target.s/max(target.s)
	
	max.I <- t(as.matrix(target.s)^2) %*% as.matrix(winD[which.max(Cps),]^2)
	if(is.na(max.I)) max.I <- 0
	max.I <- sqrt(max.I)
	if(max.I>max(winD[which.max(Cps),])) max.I <- max(winD[which.max(Cps),])
	Cps <- Cps*max.I
	
	Cps		
}

getC.rq <- function(winD, target.s)
{
	#winD <- Cmp.Matrix*hanning(nrow(Cmp.Matrix))
	#target.s <- Target.S
	
	target.s <- normalize(target.s)
	winD[,which(target.s<0.005)] <- 0

	#target.s <- target.s/max(target.s)	
	#winD[,which(target.s<0.005)] <- 0

	Cps <- apply(winD, 1, function(x) as.numeric(coefficients(rq(x~target.s))[2]))
	Cps[which(Cps<0)] <- 0
	Cps[which(is.na(Cps))] <- 0
	Cps <- normalize(Cps)
	
	max.I <- t(as.matrix(target.s)^2) %*% as.matrix(winD[which.max(Cps),]^2)
	if(is.na(max.I)) max.I <- 0
	max.I <- sqrt(max.I)
	if(max.I>max(winD[which.max(Cps),])) max.I <- max(winD[which.max(Cps),])
	Cps <- Cps*max.I	
	
	Cps			
}


get.factor.list <- function(sampleRD, analysis.window, plotting=F)
{
	from.s <- analysis.window[1]*sampleRD@scans.per.second*60 - sampleRD@start.time*sampleRD@scans.per.second
	if(trunc(from.s)==0) from.s <- 1
	to.s <- analysis.window[2]*sampleRD@scans.per.second*60 - sampleRD@start.time*sampleRD@scans.per.second
	
	if(length(intersect(c(trunc(from.s),trunc(to.s)),1:nrow(sampleRD@data)))!=2)
	{
		min.t <- sampleRD@start.time/60
		max.t <- round(nrow(sampleRD@data)/sampleRD@scans.per.second/60 + min.t,2)
		error.txt <- paste("Analysis window out of limits, time must be between ", min.t ," and ", max.t ,sep="")
		stop(error.txt)
	}
	
	## Pre-processing
	
	sampleRD@data <- sampleRD@data[trunc(from.s):trunc(to.s),]
	sampleRD@data <- pre.process(sampleRD@data, sampleRD@min.peak.width)
		
	## Parameters:

	Chrm.RawData <- sampleRD@data
	sigma.scans <-  sampleRD@min.peak.width#*0.75
	#sigma.scans <- (1/60)*5*60*0.75
	#sampleRD@min.peak.height <- 5000
	
	Scan.K <- sigma.scans	
	ksp <- vector()
	sewq <- seq(0.1,1000,0.01)
	for( i in 1:length(sewq))
	{
		krn <- normalize(dnorm(1:nrow(Chrm.RawData),trunc(nrow(Chrm.RawData)/2),sewq[i]))
		#krn <- normalize(dnorm(1:length(BIC.signal),trunc(length(BIC.signal)/2),sewq[i]))
		k.span <- trunc((length(which(normalize(krn)>0.01))/2)+0.5)	 	
		ksp[i] <- k.span
		if(k.span - Scan.K >0) break
	}
	sigma.model <- sewq[(i-1)]
	
	## Match Filter:

	krn <- normalize(dnorm(1:nrow(Chrm.RawData),trunc(nrow(Chrm.RawData)/2),sigma.model))
	k.span <- trunc((length(which(normalize(krn)>0.01))/2)+0.5)	 	
	local.kernel.o <- krn[(trunc(nrow(Chrm.RawData)/2)-k.span):(trunc(nrow(Chrm.RawData)/2)+(k.span-0))]
	local.kernel <- local.kernel.o

	ZeroPaddingMat <- matrix(0,nrow=(k.span+1), ncol=ncol(Chrm.RawData))
	Chrm.RawData.a <- rbind(ZeroPaddingMat,Chrm.RawData,ZeroPaddingMat)
	mz.index.vector <- (k.span+1):(nrow(Chrm.RawData)-k.span)
		
		pb <- txtProgressBar(min=1,max=length(mz.index.vector)*3, width=50, style=3)
		mf <- vector()
		for(i in (k.span+1):(nrow(Chrm.RawData)-k.span))
		{
			local.kernel <- local.kernel.o
			local.matrix <- Chrm.RawData.a[(i-k.span):(i+k.span),]
			local.matrix <- local.matrix#*hanning(nrow(local.matrix))
		
			cov.m <- cov(t(local.matrix))   
			is.zeroitems <- unique(c(which(apply(cov.m,1,max)==0),which(apply(cov.m,2,max)==0)))

			if(length(is.zeroitems)!=0) 
			{
				cov.m <- cov.m[-is.zeroitems,-is.zeroitems]
				local.kernel <- local.kernel[-is.zeroitems]	
				local.matrix <- local.matrix[-is.zeroitems,]
			}
		
			solv.m <- try(solve(cov.m), silent=T)
			if(class(solv.m)=="try-error") {
				mf[i] <- 0
				next
			}
				
	 		mf[i] <- sum((t(as.matrix(local.kernel)) %*% solv.m %*% as.matrix(local.matrix)) )
			setTxtProgressBar(pb, i)
		}
		
	Cmp.SetPoints <- which(diff(sign(diff(mf, na.pad = FALSE)), na.pad = FALSE) > 0) + 1
	noise.points <- which(apply(Chrm.RawData.a,1,max)[Cmp.SetPoints]<=sampleRD@min.peak.height)
	if(length(noise.points)!=0) Cmp.SetPoints <- Cmp.SetPoints[-noise.points]

	## Match Filter Correlation Filtering:
	
	pk.out <- vector()
	indx <- 1
	for(k in Cmp.SetPoints)
	{
		span.len <- k.span*2
		window.span <- (k-span.len):(k+span.len)
		if(min(window.span)<1) window.span <- window.span[-which(window.span<1)]
		if(max(window.span)>nrow(Chrm.RawData.a)) window.span <- window.span[-which(window.span>ncol(Chrm.RawData.a))]
		Cmp.Matrix <- Chrm.RawData.a[window.span,] #*hamming(length(window.span))
					
		krn <- normalize(dnorm(1:nrow(Chrm.RawData.a),k,sigma.model))
		
		SubModel <- (krn[window.span])
		CrM <- suppressWarnings(cor(SubModel,Cmp.Matrix))
		
		Kmax <- which.max(CrM)
		if(CrM[1,Kmax]<0.5)
		{
			pk.out <- c(pk.out,indx) 
		}else{
			peakMax <- (k-span.len) + (which.max(Cmp.Matrix[,Kmax]) - 1)
			Cmp.SetPoints[indx] <- peakMax
		}
		indx <- indx + 1	
	}
	if(length(pk.out)!=0) Cmp.SetPoints <- Cmp.SetPoints[-pk.out]	
	
	## MOSD:	
	
	C.matrix <- matrix(0,nrow=nrow(Chrm.RawData.a), ncol=0)
	S.matrix <- matrix(0,nrow=ncol(Chrm.RawData.a), ncol=0)
	
	#ini.pb.val <- length(Cmp.SetPoints)/2
	iteration <- length(Cmp.SetPoints)
	pb <- txtProgressBar(min=1,max=length(Cmp.SetPoints)*3, width=50, style=3)
	setTxtProgressBar(pb, iteration)
	
	for(k in Cmp.SetPoints)
	{
		span.len <- k.span*2
		window.span <- (k-span.len):(k+span.len) - 1
		if(min(window.span)<1) window.span <- window.span[-which(window.span<1)]
		if(max(window.span)>nrow(Chrm.RawData.a)) window.span <- window.span[-which(window.span>ncol(Chrm.RawData.a))]
		Cmp.Matrix <- Chrm.RawData.a[window.span,]#*hamming(length(window.span))
					
		## OSD:
		
		krn <- normalize(dnorm(1:nrow(Chrm.RawData.a),k,sigma.model))
		
		SubModel <- (krn[window.span])
		SubModel.Spectra <- getS.OSD(SubModel,Cmp.Matrix*hanning(nrow(Cmp.Matrix)),cutoff=0.15)
			
		S.cl <- SubModel.Spectra 
		#C.mod <- getC.tP(Cmp.Matrix, S.cl)
		C.mod <- try(getC.rq(Cmp.Matrix, S.cl), silent=T)
		if(class(C.mod)=="try-error") next
		
		#matplot(Cmp.Matrix, type="l", col="gray", lty=1)	
		#matplot(C.mod, type="l", lty=1, add=T)
		
		if(which.max(C.mod)==1) next
		if(which.max(C.mod)==length(C.mod)) next
		
		C.mod <- chrom.isoreg(C.mod)	
		#C.mod <- C.mod*hanning(nrow(C.mod))
			
		S.matrix <- cbind(S.matrix,as.matrix(SubModel.Spectra))
		C.matrix.aux <- matrix(0,nrow=nrow(Chrm.RawData.a), ncol=1)
		C.matrix.aux[window.span,] <- C.mod
		C.matrix <- cbind(C.matrix, C.matrix.aux)
		
		setTxtProgressBar(pb, iteration)
		iteration <- iteration + 1
	}
		
	## Duplicity Filter:

	preCouts <- which(apply(C.matrix,2,max)==0)
	if(length(preCouts)!=0)
	{
		C.matrix <- C.matrix[,-preCouts]
		S.matrix <- S.matrix[,-preCouts]
	}
	
	
	iteration <- ncol(C.matrix)*2
	pb <- txtProgressBar(min=1,max=ncol(C.matrix)*3.05, width=50, style=3)
	setTxtProgressBar(pb, iteration)
	
	CorGen <- suppressWarnings(fastCor(C.matrix))
	
	C.out.inds <- vector()
	for(j in 1:ncol(C.matrix))
	{	
		setTxtProgressBar(pb, iteration)
		iteration <- iteration + 1
			
		if(j %in% C.out.inds) next
		#if(j %in% C.ok.inds) next
		conf <- which(CorGen[,j]>0.85)
		if(length(conf)==1) next
		
		combi.mat <- sapply(1:(2^length(conf)), function(i) as.numeric(intToBits(i))[1:length(conf)])
		sel.inds <- which(normalize(rowSums(C.matrix[,conf]))>0.01)
		if(length(sel.inds)==1) 
		{
			C.out.inds <- c(C.out.inds,conf)
			next
		}
		SubSData <- Chrm.RawData.a[sel.inds,]
		maxFac <- max(SubSData)
		SubSData <- SubSData/maxFac
		SubCMod <- C.matrix[sel.inds,conf]/maxFac
		SubSMod <- S.matrix[,conf]

		recMatrix <- sapply(1:ncol(SubCMod),function(x) as.vector(SubCMod[,x] %*% t(SubSMod[,x])))
	 	
		cov.v <- vector()
		cov.e <- vector()
		for(i in 1:(ncol(combi.mat)-1))
		{
			cov.total.m <- rowSums(as.matrix(recMatrix[,which(combi.mat[,i]==1)]))	
			cov.v[i] <- sum((as.vector(SubSData) - cov.total.m)^2)	
		}
		local.ko.inds <- which(combi.mat[,which.min(cov.v)]!=1)
		C.out.inds <- c(C.out.inds,conf[local.ko.inds])
		
	}
	
	C.matrix <- C.matrix[-(1:(k.span+1)),]
	
	if(trunc(from.s)>1)
	{
		IniMat <- matrix(0,nrow=(trunc(from.s)-1), ncol=ncol(C.matrix))
		C.matrix <- rbind(IniMat,C.matrix)
	}
	
	
	## CREATING TABLE:
	
	window.features <- list(profile=C.matrix[,-C.out.inds], spectra=S.matrix[,-C.out.inds])
	
	## Peak Position
	iter.ind <- seq(1,ncol(window.features$profile))
	models.maximas <- unlist(apply(window.features$profile,2,which.max))
	#models.maximas <- models.maximas + (trunc(from.s) - 1)
	
	## Quantification
	#models.areas <- sapply(iter.ind,function(j){ sum(as.matrix(window.features$profile)[,j]%*%t(as.matrix(window.features$spectra)[,j])) })
	
	colSum.C <- colSums(window.features$profile)
	#models.areas <- sapply(iter.ind,function(j){ sum(colSum.C[j]*t(as.matrix(window.features$spectra)[,j])) })	
	models.areas <- colSum.C	
				
	## Peak Height
	models.peakheight <-  apply(as.matrix(window.features$profile),2, function(x) max(x, na.rm=T))
		
	## Spectra:
	models.spectra <- apply(as.matrix(window.features$spectra),2,function(spectra){
			spectra.index <- which(normalize(spectra)>0.005) 
			spectra.pos <- spectra.index + (sampleRD@min.mz - 1)
			spectra.int <- round(spectra[spectra.index]*1000)
			spectra.text <- paste(sweep(as.matrix(spectra.pos),1,as.matrix(spectra.int),"paste.sp"), collapse=" ")
			spectra.text
		})	

	## Profile:
	models.profile <- apply(as.matrix(window.features$profile),2,function(profile){
			profile.index <- which(normalize(profile)>0.0001) 
			profile.pos <- ((profile.index)/sampleRD@scans.per.second/60) + (sampleRD@start.time/60)
			profile.int <- profile[profile.index]/max(profile)
			profile.text <- paste(sweep(as.matrix(profile.pos),1,as.matrix(profile.int),"paste.sp"), collapse=" ")
			profile.text
		})	
	
	feature.list <- data.frame()
	
	feature.prelist <- apply(as.matrix(iter.ind),1,function(window.number){
		items <- length(models.maximas[window.number])
		subfeat.list <- matrix(0,ncol=6,nrow=items)		
		subfeat.list[,2] <- models.maximas[window.number]
		subfeat.list[,3] <- round(models.areas[window.number], digits=0)
		subfeat.list[,4] <- round(models.peakheight[window.number], digits=0)
		subfeat.list[,5] <- models.spectra[window.number]
		subfeat.list[,6] <- models.profile[window.number]
		subfeat.list
		})

	feature.list <- as.data.frame(t(feature.prelist))
	feature.list[,1] <- seq(1,nrow(feature.list))
	feature.list[,2] <- round((as.numeric(as.vector((feature.list[,2])))/sampleRD@scans.per.second)/60 + (sampleRD@start.time/60), digits=4)	
	colnames(feature.list) <- c("ID","RT","Area","Peak Height","Spectra","Profile")	
	
	setTxtProgressBar(pb, ncol(C.matrix)*3.05)
	
	feature.list


}	


