
recMissComp <- function(Experiment, min.samples, free.model=F)
{
	
	alParameters <- Experiment@Results@Parameters@Alignment
	al.par <- list(alignment.algorithm=alParameters$alignment.algorithm, min.spectra.cor=alParameters$min.spectra.cor, max.time.dist=alParameters$max.time.dist, mz.range=alParameters$mz.range, min.samples=min.samples)
	Experiment@Results@Parameters@Alignment <- al.par
	
	align.list <- as.matrix(alignList(Experiment))
	data.list <- as.matrix(align.list[,-c(1:4)])	
	
	## Delete the peaks below the minimum peak threshold
	delete.indexes <- which(as.numeric(as.vector(apply(data.list,1,function(x) max(as.numeric(as.vector(x))))))<Experiment@Data@Parameters$min.peak.height)		
	if(length(delete.indexes)!=0) {data.list <- data.list[-delete.indexes,]; align.list <- align.list[-delete.indexes,]}
	

	meta.data <- metaData(Experiment)
	sample.filenames <- as.vector(apply(as.matrix(colnames(data.list)),1,function(x) meta.data[which(x==meta.data),"filename"]))
	
	model.indexes <- which(as.numeric(as.vector(align.list[,"FoundIn"]))>=min.samples)

	#if(length(re.analyze.indexes)==0) stop("Model already fitted for the min.samples value selected") 
	times.list <-  align.list[model.indexes,c("AlignID","tmean")]
	
	total.iterations <- 0
	for(sample in 1:length(sample.filenames)) total.iterations <- total.iterations + length(data.list[model.indexes,sample])

	cat("\n")
	pb <- txtProgressBar(min=1,max=total.iterations, width=50, style=3)
	k.it <- 0
	for(sample in 1:length(sample.filenames))
	{
		## Error on sample 23
    
	  re.analyze.in <- which(as.numeric(data.list[model.indexes,sample])==0)
	  if(length(re.analyze.in)!=0) 
	  {
	    re.analyze <- model.indexes[re.analyze.in]
	    
    
	 		if(Experiment@MetaData@DataDirectory=="") {filename <- as.character(Experiment@MetaData@Instrumental$filename[sample])
			}else{filename <- paste(Experiment@MetaData@DataDirectory,"/",Experiment@MetaData@Instrumental$filename[sample], sep="")}
			
			sampleRD <- load.file(filename)
		
			sampleRD@scans.per.second <- Experiment@Data@Parameters$scans.per.second
			sampleRD@compression.coef <- Experiment@Data@Parameters$compression.coef
			sampleRD@avoid.processing.mz <- Experiment@Data@Parameters$avoid.processing.mz 
			sampleRD@min.peak.width <- Experiment@Data@Parameters$min.peak.width 
			sampleRD@noise.threshold <- Experiment@Data@Parameters$noise.threshold 
			sampleRD@min.peak.height <- Experiment@Data@Parameters$min.peak.height
			#sampleRD@factor.minimum.sd <- Experiment@Data@Parameters$factor.minimum.sd 
			
			avoid.mz <- sampleRD@avoid.processing.mz - (sampleRD@min.mz - 1)
			sampleRD@data[,avoid.mz] <- 0
	
			moving.maximas <- apply(sampleRD@data,2,function(x){max(x)})
			sampleRD@data[,moving.maximas<sampleRD@noise.threshold] <- 0


			#sigma.object <- get.SigmaModel(sampleRD)
			#sigma.model <- sigma.object$sigma
			#sigma.scans <-  sigma.object$sigma.scans
			
			sigma.scans <-Experiment@Data@Parameters$min.peak.width*Experiment@Data@Parameters$scans.per.second*60
			sampleRD@data <- pre.process(sampleRD@data, sigma.scans)
			#sampleRD@min.peak.width <- sigma.object$sigma
			
			analyzed.vector <- vector()

			for(j in 1:length(re.analyze))
			{
				k.it <- k.it+1

				lost.FactorID <-  as.numeric(as.vector(align.list[re.analyze[j],"AlignID"]))
				lost.factor <- fit.model.in.data.tosd(sampleRD, Experiment,lost.FactorID, free.model)

				if(is.null(lost.factor)) next
				
					free.factor.index <- which(Experiment@Data@FactorList[colnames(data.list)[sample]][[1]]$"AlignID"==lost.FactorID)
					if(length(free.factor.index)!=0) {Experiment@Data@FactorList[colnames(data.list)[sample]][[1]]$"AlignID"[free.factor.index] <- 0}
			
					new.list <- as.data.frame(matrix(0,ncol=7,nrow=1))
					colnames(new.list) <- c("ID","RT","Area","Peak Height","Spectra","Profile","AlignID")
						
					new.list$"ID" <- length(Experiment@Data@FactorList[colnames(data.list)[sample]][[1]]$"ID")+1
					new.list$"Peak Height" <- lost.factor$height
					new.list$"RT" <- as.numeric(lost.factor$rt)
					new.list$"Area" <- lost.factor$area
					new.list$"Spectra" <- lost.factor$lost.spectra
					new.list$"Profile" <- lost.factor$lost.profile
					new.list$"AlignID" <- lost.FactorID 
				
					Experiment@Data@FactorList[colnames(data.list)[sample]][[1]] <- data.frame(rbind(as.matrix(Experiment@Data@FactorList[colnames(data.list)[sample]][[1]]),as.matrix(new.list)))		
					Experiment@Data@FactorList[colnames(data.list)[sample]][[1]]$"ID" <- as.integer(Experiment@Data@FactorList[colnames(data.list)[sample]][[1]]$"ID")
					Experiment@Data@FactorList[colnames(data.list)[sample]][[1]]$"RT" <- as.numeric(as.vector(Experiment@Data@FactorList[colnames(data.list)[sample]][[1]]$"RT"))
					Experiment@Data@FactorList[colnames(data.list)[sample]][[1]]$"AlignID" <- as.numeric(as.vector(Experiment@Data@FactorList[colnames(data.list)[sample]][[1]]$"AlignID"))
					colnames(Experiment@Data@FactorList[colnames(data.list)[sample]][[1]]) <- c("ID", "RT", "Area", "Peak Height", "Spectra", "Profile", "AlignID")						
				
				setTxtProgressBar(pb, k.it)
			}
		}	
	}
	setTxtProgressBar(pb, total.iterations)
	cat("\n Updating alignment table... \n")
	Experiment@Results@Alignment <- create.factorlist.table(Experiment)
	cat("Model fitted! \n")
	Experiment	
	
}




fit.model.in.data.tosd <- function(sampleRD, Experiment, lost.factor.alignId, free.model)
{

	maxMZ <- sampleRD@max.mz
	minMZ <- sampleRD@min.mz
	
	lost.factor.index <- which(Experiment@Results@Alignment[,"AlignID"]==lost.factor.alignId)	
	lost.factor.position <- Experiment@Results@Alignment[lost.factor.index,"tmean"]*sampleRD@scans.per.second*60-sampleRD@start.time*sampleRD@scans.per.second
		
	min.peak.width <- sampleRD@min.peak.width
	max.time.dist <- Experiment@Results@Parameters@Alignment$max.time.dist*sampleRD@scans.per.second*60
	
	window.start <- as.integer(lost.factor.position - max.time.dist - min.peak.width*4)
	if(window.start<=0) window.start <- 1
	if(window.start>nrow(sampleRD@data)) window.start <- as.integer(nrow(sampleRD@data) - max.time.dist - min.peak.width*4)
	window.end <- as.integer(lost.factor.position + max.time.dist + min.peak.width*4)
	if(window.end>nrow(sampleRD@data)) window.end <- nrow(sampleRD@data)
	
	analysis.window <- sampleRD@data[window.start:window.end,]
	
	if(max(analysis.window)<Experiment@Data@Parameters$noise.threshold) return(NULL)
	
	## Lost Spectra:
	
	lost.spectra <- rep(0,(maxMZ-minMZ)+1)
	lost.spectra <- normalize(convertMSPspectra(Experiment@Results@Alignment[lost.factor.index,"Spectra"], maxMZ)[minMZ:maxMZ])

	## Lost Profile:
	
	#C.model <- getC.tOSD(analysis.window, lost.spectra)
	#C.model <- getC.tP(analysis.window, lost.spectra)
	#C.model <- getC.rq(analysis.window, lost.spectra)
	C.model <- try(getC.rq(analysis.window, lost.spectra), silent=T)
	if(class(C.model)=="try-error") C.model <- rep(0, nrow(analysis.window))

	C.model <- chrom.isoreg(C.model)
		
	Smod <- getS.OSD(C.model, analysis.window, ref.response=lost.spectra)
	# if(is.na(suppressWarnings(cor(Smod, lost.spectra))<0.75)) 
	# {
		# max.mz <- apply(analysis.window,2,max)
		# max.mz[max.mz==0] <- max(max.mz)
		# C.model <- analysis.window[,which.min(max.mz)]
	# }else{
		# if(suppressWarnings(cor(Smod, lost.spectra))<0.75) 
		# {
			# max.mz <- apply(analysis.window,2,max)
			# max.mz[max.mz==0] <- max(max.mz)
			# C.model <- analysis.window[,which.min(max.mz)]
		# }  
	# }
	
	if(free.model==T) lost.spectra <- Smod
			
	lost.profile <- C.model		
			
	## Results Formatting
	
	#t.vect <- (window.start:window.end)/sampleRD@scans.per.second/60 + sampleRD@start.time/60
	#matplot(t.vect, analysis.window, type="l", col="gray", lty=1)
	#matplot(t.vect, lost.profile, type="l", add=T, lty=1)

	return.list <- list()

		lost.rt <- ((window.start-1) + which.max(lost.profile))/sampleRD@scans.per.second/60 + sampleRD@start.time/60
		#lost.area <- sum(as.matrix(lost.profile) %*% t(as.matrix(lost.spectra)))
		lost.area <- sum(as.matrix(lost.profile))
		
		lost.spectra[normalize(lost.spectra)<0.005] <- 0
		spectra.index <- which(lost.spectra!=0) 
		spectra.pos <- spectra.index + (sampleRD@min.mz - 1)
		spectra.int <- round(lost.spectra[spectra.index]*1000)
		del.pos <- which(spectra.int==0)
		if(length(del.pos)!=0) {spectra.int <- spectra.int[-del.pos]; spectra.pos <- spectra.pos[-del.pos]}
		spectra.text <- paste(sweep(as.matrix(spectra.pos),1,as.matrix(spectra.int),"paste.sp"), collapse=" ")
		
		lost.height <- lost.profile[which.max(lost.profile)]
		lost.profile[normalize(lost.profile)<0.005] <- 0
		profile.index <- which(lost.profile!=0) 
		profile.pos <- ((profile.index + (window.start - 1))/sampleRD@scans.per.second/60) + (sampleRD@start.time/60)
		profile.int <- lost.profile[profile.index]/max(lost.profile)
		profile.text <- paste(sweep(as.matrix(profile.pos),1,as.matrix(profile.int),"paste.sp"), collapse=" ")

		return.list <- list(lost.profile=profile.text, lost.spectra=spectra.text, rt=lost.rt, area=lost.area, height=lost.height)	

	#lines(c(t.vect[lost.factor.index],t.vect[lost.factor.index]),c(0,max(analysis.window)), col="red")
	return.list	
}

