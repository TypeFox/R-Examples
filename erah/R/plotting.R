plotSpectra <- function(Experiment, AlignId, n.putative=1, compare=T, id.database=mslib, comp.db=NULL, return.spectra=F, draw.color="purple", xlim=NULL)
{
	if(length(AlignId)!=1) stop("Only one spectrum can be shown at once")
	if(compare==T) if(is.null(id.database)) stop("A database is needed for spectra comparison. Select a database or set 'compare' parameter to 'False'")

	index <- which(as.numeric(as.vector(Experiment@Results@Alignment[,"AlignID"]))==AlignId)
	MSP.spect.emp <- Experiment@Results@Alignment[index,"Spectra"]

	current.column <- paste("DB.Id.",n.putative, sep="")

	if(compare==T) 
	{
		if(is.null(comp.db)) if(nrow(Experiment@Results@Identification)==1) stop("Factors must be identified first")
		if(is.null(comp.db)) MSP.spect.db <- id.database@database[[as.numeric(as.vector(Experiment@Results@Identification[index,current.column]))]]$Spectra
		if(!is.null(comp.db)) MSP.spect.db <- id.database@database[[comp.db]]$Spectra
		
		maxMZ <- max(Experiment@Results@Parameters@Identification$compare.only.mz)
		empiric.name <- id.database@database[[as.numeric(as.vector(Experiment@Results@Identification[index,current.column]))]]$Name
		
		empiric.spectra <- convertMSPspectra(MSP.spect.emp,maxMZ)
	
		delete.mz <- 1:length(empiric.spectra)
		delete.mz <- delete.mz[-Experiment@Results@Parameters@Identification$compare.only.mz]
		if(length(delete.mz)!=0) empiric.spectra[delete.mz] <- 0
		
		empiric.spectra <- empiric.spectra/max(empiric.spectra)*1000
			
	}else{
		maxMZ <- max(Experiment@Results@Parameters@Alignment$mz.range)
		empiric.name <- paste("Factor #",AlignId, sep="")
		
		empiric.spectra <- convertMSPspectra(MSP.spect.emp,maxMZ)
		empiric.spectra <- empiric.spectra/max(empiric.spectra)*1000

	}
		
	mz.len <- 5
	if(length(which(empiric.spectra!=0))<5) mz.len <- length(which(empiric.spectra!=0))
	main_mz.empiric <- sort(empiric.spectra, decreasing=T, index.return=T)$ix[1:mz.len]

	if(compare==F)
	{	
		plot(empiric.spectra, type="h", main=paste(empiric.name, "\n (Empiric Spectra)"), xlab="Mz", ylab="Intensity",xlim=xlim)
		text(main_mz.empiric, empiric.spectra[main_mz.empiric], labels=main_mz.empiric, cex=0.7)	
	}else{
		
		db.spectra <- convertMSPspectra.dot(MSP.spect.db,maxMZ)
		delete.mz1 <- Experiment@Data@Parameters$avoid.processing
		delete.mz2 <- 1:length(empiric.spectra)
		delete.mz2 <- delete.mz2[-Experiment@Results@Parameters@Identification$compare.only.mz]
		delete.mz <- unique(c(delete.mz1,delete.mz2))
		if(length(delete.mz)!=0) db.spectra[delete.mz] <- 0
		db.spectra <- normalize(db.spectra)*(-1000)

		match.factor <- cor.sinus(empiric.spectra,abs(db.spectra))
		match.factor <- round(match.factor*100, digits=1)
		
		if(is.null(comp.db)) main.title <- paste(empiric.name, "\n Match Factor:",match.factor)
		if(!is.null(comp.db)) main.title <- paste("Align ID #", AlignId, " VS ", id.database@database[[comp.db]]$Name, "\n Match Factor:",match.factor, sep="")
		
		plot(empiric.spectra, type="h", main=main.title, ylim=c(-1000,1000), xlab="Mz", ylab="Intensity", xlim=xlim)
		text(main_mz.empiric, empiric.spectra[main_mz.empiric], labels=main_mz.empiric, cex=0.7)	
		lines(db.spectra, col=draw.color, type="h")
		main_mz.db <- sort(abs(db.spectra), decreasing=T, index.return=T)$ix[1:mz.len]
		text(main_mz.db, db.spectra[main_mz.db], labels=main_mz.db, cex=0.7, col=draw.color)	
		if(is.null(comp.db)) legend("topright", legend=c("Empiric","Database"), col=c("black",draw.color), pch=19)
		if(!is.null(comp.db)) legend("topright", legend=c(paste("Align ID #", AlignId, sep=""),id.database@database[[comp.db]]$Name), col=c("black",draw.color), pch=19)
	}
	if(return.spectra==T) return(empiric.spectra)
	
}

plotProfile <- function(Experiment,AlignId, per.class=T, xlim=NULL)
{	
	if(!(any(unlist(lapply(Experiment@Data@FactorList,function(x) {is.null(x$AlignID)} ))==FALSE))) stop("Factors must be aligned first")
	
	if(nrow(Experiment@MetaData@Phenotype)==0) 
		if(per.class==T) {per.class=F; warning("The experiment does not contain phenotypic metadata. Profiles are shown per sample.")}
	
	empty.samples <- which(lapply(Experiment@Data@FactorList,nrow)==0)
	if(length(empty.samples)!=0)
	{
		Experiment@Data@FactorList <- Experiment@Data@FactorList[-empty.samples]
		Experiment@MetaData@Phenotype <- Experiment@MetaData@Phenotype[-empty.samples,]
	}
	
	alignId <- lapply(Experiment@Data@FactorList,function(x){x$AlignID})
	N.groups <- max(unique(unlist(alignId)))	
	N.samples <- length(Experiment@Data@FactorList)
	
	samples.name <- names(Experiment@Data@FactorList)		
	for(i in 1:N.samples) samples.name[i] <- strsplit(as.character(samples.name[i]), split="\\.")[[1]][1]
	
	profile.list <- lapply(Experiment@Data@FactorList,function(x) {
			outp <- as.character(x[which(x$AlignID==AlignId),"Profile"])
			time <- NA
			int <- NA
			if(length(outp)!=0)
			{
				output <- sparse.to.vector(outp)
				time <- output$time
				int <- output$int*as.numeric(as.character(x[which(x$AlignID==AlignId),"Peak Height"]))
			}
			list(time=time,int=int)
			})
	
	profile.len <- max(unlist(lapply(unlist(profile.list, recursive=F),length)))
	
	profile.time = matrix(NA, nrow=profile.len,ncol=N.samples)
	profile.int = matrix(NA, nrow=profile.len,ncol=N.samples)

	for(i in 1:N.samples) {profile.time[1:length(profile.list[[i]]$time),i] <- profile.list[[i]]$time; profile.int[1:length(profile.list[[i]]$int),i] <- profile.list[[i]]$int}
	
	na.samples.i <- which(apply(apply(profile.int,2,is.na),2,all)==T)
	na.samples.t <- which(apply(apply(profile.time,2,is.na),2,all)==T)
	na.samples <- unique(na.samples.i,na.samples.t)
	
	if(length(na.samples)!=0)
	{ 
		samples.name <- samples.name[-na.samples]
		profile.int <- profile.int[,-na.samples]
		profile.time <- profile.time[,-na.samples]
	}
	
	vector.time <- diag(profile.time[apply(profile.int,2,which.max),])
	time.mean <- mean(vector.time)
	
	profile.time <- sweep(profile.time,2,(vector.time-time.mean),"-") 
	
	compound.name <- as.character(Experiment@Results@Identification[which(Experiment@Results@Identification$AlignID==AlignId),"Name"])
	
	if(per.class==F)
	{
		
		matplot(profile.time, profile.int, type="l", lty=1, col=(1:length(samples.name)), main=paste("Profile Comparison \n",compound.name), xlab="time (min)", ylab="Intensity", xlim=xlim)
		par(font=2)
		legend("topright",legend=samples.name, pch=19, col=(1:length(samples.name)), title="Samples")
		par(font=1)
	}else{
		pn <- Experiment@MetaData@Phenotype
		indx <- apply(as.matrix(samples.name),1,function(x) which(pn[,"sampleID"]==x))		
		class.names <- pn[indx,"class"]
		
		samples.class.type <- levels(pn$class)
			
		matplot(profile.time,profile.int, type="l", lty=1, col=class.names, main=paste("Profile Comparison \n",compound.name), xlab="time (min)", ylab="Intensity", xlim=xlim)
		par(font=2)
		legend("topright",legend=samples.class.type, pch=19, col=1:length(samples.class.type), title="Classes")
		par(font=1)

	}
	
}


plotAlign <- function(Experiment, AlignId, per.class=T, xlim=NULL)
{	
	if(!(any(unlist(lapply(Experiment@Data@FactorList,function(x) {is.null(x$AlignID)} ))==FALSE))) stop("Factors must be aligned first")
	
	if(nrow(Experiment@MetaData@Phenotype)==0) 
	{
		per.class=F
		warning("The experiment does not contain phenotypic metadata. Profiles are shown per sample.")	
	}
	
	empty.samples <- which(lapply(Experiment@Data@FactorList,nrow)==0)
	if(length(empty.samples)!=0)
	{
		Experiment@Data@FactorList <- Experiment@Data@FactorList[-empty.samples]
		Experiment@MetaData@Phenotype <- Experiment@MetaData@Phenotype[-empty.samples,]
	}
	
	alignId <- lapply(Experiment@Data@FactorList,function(x){x$AlignID})
	N.groups <- max(unique(unlist(alignId)))	
	N.samples <- length(Experiment@Data@FactorList)
	
	samples.name <- names(Experiment@Data@FactorList)		
	for(i in 1:N.samples) samples.name[i] <- strsplit(as.character(samples.name[i]), split="\\.")[[1]][1]
	
	profile.list <- lapply(Experiment@Data@FactorList,function(x) {
			outp <- as.character(x[which(x$AlignID==AlignId),"Profile"])
			time <- NA
			int <- NA
			if(length(outp)!=0)
			{
				output <- sparse.to.vector(outp)
				time <- output$time
				int <- output$int*as.numeric(as.character(x[which(x$AlignID==AlignId),"Peak Height"]))
			}
			list(time=time,int=int)
			})
	
	profile.len <- max(unlist(lapply(unlist(profile.list, recursive=F),length)))
	
	profile.time = matrix(NA, nrow=profile.len,ncol=N.samples)
	profile.int = matrix(NA, nrow=profile.len,ncol=N.samples)

	for(i in 1:N.samples) {profile.time[1:length(profile.list[[i]]$time),i] <- profile.list[[i]]$time; profile.int[1:length(profile.list[[i]]$int),i] <- profile.list[[i]]$int}
	
	na.samples <- which(apply(apply(profile.time,2,is.na),2,all)==T)
	if(length(na.samples)!=0)
	{ 
		samples.name <- samples.name[-na.samples]
		profile.int <- profile.int[,-na.samples]
		profile.time <- profile.time[,-na.samples]
	}
	
	vector.time <- diag(profile.time[apply(profile.int,2,which.max),])
	time.mean <- mean(vector.time)
	
	profile.un.time <- profile.time
	profile.time <- sweep(profile.time,2,(vector.time-time.mean),"-") 
	
	compound.name <- as.character(Experiment@Results@Identification[which(Experiment@Results@Identification$AlignID==AlignId),"Name"])
	
	if(per.class==F)
	{	
		par(mfrow=c(1,2))
		matplot(profile.un.time, profile.int, type="l", lty=1, col=(1:length(samples.name)), main=paste("Unaligned \n",compound.name), xlab="time (min)", ylab="Intensity", xlim=xlim)
		par(font=2)
		legend("topright",legend=samples.name, pch=19, col=(1:length(samples.name)), title="Samples")
		par(font=1)
		matplot(profile.time, profile.int, type="l", lty=1, col=(1:length(samples.name)), main=paste("Aligned \n",compound.name), xlab="time (min)", ylab="Intensity", xlim=xlim)
		par(font=2)
		legend("topright",legend=samples.name, pch=19, col=(1:length(samples.name)), title="Samples")
		par(font=1)
	}else{
		pn <- Experiment@MetaData@Phenotype
		indx <- apply(as.matrix(samples.name),1,function(x) which(pn[,"sampleID"]==x))		
		class.names <- pn[indx,"class"]
		
		samples.class.type <- levels(pn$class)
		
		par(mfrow=c(1,2))
		matplot(profile.un.time, profile.int, type="l", lty=1, col=class.names, main=paste("Unaligned \n",compound.name), xlab="time (min)", ylab="Intensity", xlim=xlim)
		par(font=2)
		legend("topright",legend=samples.class.type, pch=19, col=1:length(samples.class.type), title="Classes")
		par(font=1)
		matplot(profile.time, profile.int, type="l", lty=1, col=class.names, main=paste("Aligned \n",compound.name), xlab="time (min)", ylab="Intensity", xlim=xlim)
		par(font=2)
		legend("topright",legend=samples.class.type, pch=19, col=1:length(samples.class.type), title="Classes")
		par(font=1)
	}
	
}


plotChr <- function(Experiment, N.sample=1, type=c("BIC","TIC","EIC"), xlim=NULL, mz=NULL)
{
	type <- match.arg(type, c("BIC","TIC","EIC"), several.ok = FALSE)
	
	sampleRD <- load.file(paste(Experiment@MetaData@DataDirectory, Experiment@MetaData@Instrumental$filename[[N.sample]], sep="/"))
	
	max.rt <- (nrow(sampleRD@data)/(sampleRD@scans.per.second*60)) + sampleRD@start.time/60
	min.rt <- sampleRD@start.time/60
	vect.rt <- seq(min.rt, max.rt, length.out=nrow(sampleRD@data))
	
	if(type=="TIC") plot(vect.rt, rowSums(sampleRD@data), type="l", main="TIC", xlab="RT(min)", ylab="TIC", xlim)
	if(type=="BIC") plot(vect.rt, apply(sampleRD@data,1,max), type="l", main="BIC", xlab="RT(min)", ylab="BIC", xlim)	
	if(type=="EIC") 
	{
		if(is.null(mz)) stop("When plotting the EIC, please, specify the masses or range of masses to be plotted.")
		if(any(mz<sampleRD@min.mz) || any(mz>sampleRD@max.mz)) stop(paste("The adquisition m/z range selected is above or under the selected masses. For this sample, the adquistion range is from ", sampleRD@min.mz, " to ", sampleRD@max.mz, ". Please, change the m/z paramter accordingly, inside the limits.", sep=""))

		mz.rang <- mz - (sampleRD@min.mz - 1)
	    if(!is.null(xlim)) scan.range <- c(which.min(abs(vect.rt - min(xlim))):which.min(abs(vect.rt - max(xlim))))
        if(is.null(xlim)) scan.range <- 1:nrow(sampleRD@data)
		matplot(vect.rt[scan.range], sampleRD@data[scan.range,mz.rang], type="l", lty=1, main="EIC", xlab="RT(min)", ylab="EIC")	
		legend("topright",legend=mz, pch=19, col=1:length(mz), title="m/z")

	}
}



# plotBoxplot <- function(Experiment, AlignId, classes.to.compare=NULL, outline = TRUE, log = "", horizontal = FALSE)
# {
	# if(!(any(unlist(lapply(Experiment@Data@FactorList,function(x) {is.null(x$AlignID)} ))==FALSE))) stop("Factors must be aligned first")
	# if(nrow(Experiment@MetaData@Phenotype)==0) stop("No Phenotype data has been attached to this experiment.")
	
	# align.list <- NULL
	# align.list <- alignList(Experiment)
	# if(is.null(align.list)) return()
	
	# empty.samples <- which(lapply(Experiment@Data@FactorList,nrow)==0)
	# if(length(empty.samples)!=0)
	# {
		# Experiment@Data@FactorList <- Experiment@Data@FactorList[-empty.samples]
		# Experiment@MetaData@Phenotype <- Experiment@MetaData@Phenotype[-empty.samples,]
	# }
	
	# pn <- Experiment@MetaData@Phenotype
	# samples.name <- names(Experiment@Data@FactorList)		
	# indx <- apply(as.matrix(samples.name),1,function(x) which(pn[,"sampleID"]==x))		
	# class.names <- as.vector(pn[indx,"class"])
		
	# samples.class.type <- levels(as.factor(class.names))
	
	# data.list <- align.list[,-c(1:4)]

	# if(!is.null(classes.to.compare))
	# {
		# classes.selected <- unique(unlist(apply(as.matrix(classes.to.compare),1,function(x) which(x==samples.class.type ))))
		# if(length(classes.selected)!=length(classes.to.compare)) stop("Invalid class selected. The class name is invalid or no sample of this class has been processed. Use expClasses() function to obtain more information about the existing classes in this experiment")
		# sel.classes.index <- unlist(apply(as.matrix(classes.to.compare),1,function(x) which(x==class.names)))
		# class.names <- class.names[sel.classes.index]
		# data.list <- data.list[,sel.classes.index]
	# }
	
	# p.value <- round(summary(aov(as.numeric(as.vector(as.matrix(data.list[which(align.list[,"AlignID"]==AlignId),])))~as.factor(class.names)))[[1]][1,5], digits=4)
	# main.title <- paste(align.list[which(align.list[,"AlignID"]==AlignId),"Factor"], "\n p-value = ", round(p.value, digits=4), sep="")
	
	# boxplot(as.numeric(as.vector(as.matrix(data.list[which(align.list[,"AlignID"]==AlignId),])))~as.factor(class.names), outline = outline, log = log, horizontal = horizontal, main=main.title, ylab="Peak Height")

# }

# plotPca <- function(Experiment)
# {
	# if(nrow(Experiment@MetaData@Phenotype)==0) stop("No Phenotype data has been attached to this experiment.")
	
	# align.list <- NULL
	# align.list <- alignList(Experiment)
	# if(is.null(align.list)) return()
	
	# empty.samples <- which(lapply(Experiment@Data@FactorList,nrow)==0)
	# if(length(empty.samples)!=0)
	# {
		# Experiment@Data@FactorList <- Experiment@Data@FactorList[-empty.samples]
		# Experiment@MetaData@Phenotype <- Experiment@MetaData@Phenotype[-empty.samples,]
	# }
	
	# pn <- Experiment@MetaData@Phenotype
	# samples.name <- names(Experiment@Data@FactorList)		
	# indx <- apply(as.matrix(samples.name),1,function(x) which(pn[,"sampleID"]==x))		
	# class.names <- as.vector(pn[indx,"class"])
		
	# samples.class.type <- levels(as.factor(class.names))
	
	# data.list <- align.list[,-c(1:3)]
	
	# PCA.analysis <- prcomp(matrix(as.numeric(as.matrix(data.list)), ncol=ncol(data.list)), scale=F, center=T)
	
	# x.lab <- paste("Scores on PC 1 (", round(summary(PCA.analysis)$importance[2,1]*100, digits=1), "%)", sep="")
	# y.lab <- paste("Scores on PC 2 (", round(summary(PCA.analysis)$importance[2,2]*100, digits=1), "%)", sep="")
	# plot(PCA.analysis$rotation[,1],PCA.analysis$rotation[,2], pch=20, col=as.factor(class.names), xlab=x.lab, ylab=y.lab)
	# legend("topright",legend=samples.class.type, pch=19, col=1:length(samples.class.type), title="Classes")
# }