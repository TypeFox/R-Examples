idList <- function(object, id.database=mslib) {
	#if(!(any(unlist(lapply(object@Data@FactorList,function(x) {is.null(x$AlignID)} ))==FALSE))) stop("Factors must be aligned and identified first")
	if(nrow(object@Results@Identification)==0) stop("Factors must be identified first")
	if(object@Results@Parameters@Identification$database.name!=id.database@name)
	{
		error.msg <- paste("This experiment was not processed with the database selected. Please use ", object@Results@Parameters@Identification$database.name, sep="")
		stop(error.msg)
	} 

	n.putative <- object@Results@Parameters@Identification$n.putative
	
	id.general <- object@Results@Identification[,c("AlignID","tmean","FoundIn")]	
	id.list.up <- NULL	
	
	for(current.putative in 1:n.putative)
	{
		current.column <- paste("DB.Id.",current.putative, sep="")
		current.column.match <- paste("MatchFactor.",current.putative, sep="")
		current.column.cas <- paste("CAS.",current.putative, sep="")
		current.column.form	<- paste("Formula.",current.putative, sep="")
		current.column.name	<- paste("Name.",current.putative, sep="")

		id.found <- as.numeric(as.vector(object@Results@Identification[,current.column]))
		met.name <- unlist(lapply(id.database@database[c(id.found)], function(x) x$Name))
		met.cas <- unlist(lapply(id.database@database[c(id.found)], function(x) x$CAS))
		met.form <- unlist(lapply(id.database@database[c(id.found)], function(x) x$Formula))
				
		id.list <- object@Results@Identification[,c(current.column.match,current.column)]
		del.na <- which(is.na(id.list[,1]))
		if(length(del.na)!=0) id.list <- id.list[-del.na,]
		id.list <- cbind(id.list,met.name,met.cas,met.form)
		colnames(id.list)[ncol(id.list):(ncol(id.list)-2)] <- c(current.column.form,current.column.cas,current.column.name)
		cas.no.na <- apply(id.list,1,function(x) {if(as.vector(x[current.column.cas])=="NA"){
		 	x[current.column.cas]=""		
		} else{
			x[current.column.cas]	
		}
		})
		id.list[,current.column.cas] <- cas.no.na
		id.list <- id.list[,c(current.column.name,current.column.match,current.column,current.column.cas,current.column.form)]
		id.list.up <- cbind(id.list.up,as.matrix(id.list))
	}
	if(length(del.na)!=0) id.general <- id.general[-del.na,]
	id.list <- cbind(id.general,id.list.up)	
	if(!is.null(object@Results@Parameters@Alignment$min.samples)){
		id.list <- id.list[which(id.list$FoundIn>=object@Results@Parameters@Alignment$min.samples),]
	}
		
	return(as.data.frame(id.list[order(as.numeric(as.vector(id.list[,"tmean"]))),], row.names=1:nrow(id.list)))
}

# factorList <- function(object, sample)
# {
	# if(is.null(object@Data@FactorList[[sample]]$AlignID)) 
	# {
		# return(as.data.frame(object@Data@FactorList[[sample]][,c("ID","RT","Peak Height","Area")]))
	# }else{
		# return(as.data.frame(object@Data@FactorList[[sample]][,c("ID","AlignID","RT","Peak Height","Area")]))
	# }
# }

alignList <- function(object, by.area=TRUE) {
		
	if(!(any(unlist(lapply(object@Data@FactorList,function(x) {is.null(x$AlignID)} ))==FALSE))) stop("Factors must be aligned first")
	height <- FALSE
	if(by.area==FALSE) height <- TRUE
	if(height==TRUE)
	{
		del.in <- which(colnames(object@Results@Alignment)=="Spectra")
		alList <- object@Results@Alignment[,-del.in]
		if(!is.null(object@Results@Parameters@Alignment$min.samples)){
			alList  <- alList[which(alList$FoundIn>=object@Results@Parameters@Alignment$min.samples),]
		}		
		return(alList)
	}else{
		# align.inds <- as.numeric(as.vector(object@Results@Alignment[,"AlignID"]))
		# align.area <- lapply(object@Data@FactorList,function(x) {
			# #search.for <- which(x$"AlignID" %in% align.inds)
			# search.for <- unlist(sapply(align.inds, function(i) which(as.numeric(as.vector(x$"AlignID"))==i)))
			# as.numeric(as.vector(x$"Area"[search.for]))
		# })
		# if(!is.matrix(align.area))
		# {
			# del.list <- as.vector(which(unlist(lapply(align.area,length))==0))
			# if(length(del.list)!=0) align.area <- align.area[-del.list]
			# align.area <- unlist(align.area)
			# dim(align.area) <- c(length(align.inds), length(align.area)/length(align.inds))
		# }
		# del.in <- which(colnames(object@Results@Alignment)=="Spectra")
		# area.list <- cbind(object@Results@Alignment[,c("AlignID","Factor","tmean","FoundIn")],align.area)
		# colnames(area.list) <- colnames(object@Results@Alignment[,-del.in])
		# return(area.list)
		
		align.inds <- as.numeric(as.vector(object@Results@Alignment[,"AlignID"]))
		align.area <- lapply(object@Data@FactorList,function(x) {
			#search.for <- which(x$"AlignID" %in% align.inds)
			search.for <- (sapply(align.inds, function(i) which(as.numeric(as.vector(x$"AlignID"))==i)))
			if(class(search.for)=="list") 
			{
				fill.inds <- unlist(lapply(search.for, length))
				fill.vector <- rep(0,length(align.inds))
				fill.vector[fill.inds==1] <- as.numeric(as.vector(x$"Area"[unlist(search.for)]))
			}else{
				fill.vector <- as.numeric(as.vector(x$"Area"[search.for]))
			}
			fill.vector
		})
		
		if(!is.matrix(align.area))
		{
			del.list <- as.vector(which(unlist(lapply(align.area,length))==0))
			if(length(del.list)!=0) align.area <- align.area[-del.list]
			align.area <- unlist(align.area)
			dim(align.area) <- c(length(align.inds), length(align.area)/length(align.inds))
		}
		del.in <- which(colnames(object@Results@Alignment)=="Spectra")
		area.list <- cbind(object@Results@Alignment[,c("AlignID","Factor","tmean","FoundIn")],align.area)
		colnames(area.list) <- colnames(object@Results@Alignment[,-del.in])
		
		if(!is.null(object@Results@Parameters@Alignment$min.samples)){
			area.list  <- area.list[which(area.list$FoundIn>=object@Results@Parameters@Alignment$min.samples),]
		}		
		
		return(area.list)
	}
}

expClasses <- function(object)
{
	if(nrow(object@MetaData@Phenotype)==0) stop("No Phenotype data has been attached to this experiment.")
	
	pn <- object@MetaData@Phenotype
	samples.name <- names(object@Data@FactorList)		
	indx <- apply(as.matrix(samples.name),1,function(x) which(pn[,"sampleID"]==x))		
	class.names <- as.vector(pn[indx,"class"])
	samples.class.type <- levels(as.factor(class.names))
	empty.samples <- which(lapply(object@Data@FactorList,nrow)==0)
	
	classes.list <- matrix(c(samples.name, class.names, rep("Processed", length(samples.name))), ncol=3)
	classes.list[empty.samples,3] <- "Not processed"
	colnames(classes.list) <- c("Sample ID", "Class Type", "Processing Status")
	
	classes.summary <- as.data.frame(classes.list, row.names=1:nrow(classes.list))
	
	classes.string <- paste(samples.class.type, collapse=", ")
	cat("Experiment containing ", nrow(classes.summary), " samples in ", length(samples.class.type), " different type of classes named: ",classes.string, ". \n \n", sep="")
	print(classes.summary)
	
	#return(new("expClasses", classes.type=samples.class.type, classes.summary=classes.summary))
}

dataList <- function(Experiment, id.database=mslib, by.area=TRUE)
{
	ID.table <- idList(Experiment, id.database)	
	Al.table <- alignList(Experiment, by.area)
	Al.table <- Al.table[, !(names(Al.table) %in% c("FoundIn","tmean","Factor"))]	
	data.table <- merge(ID.table, Al.table, by="AlignID")
	data.table
}

sampleInfo <- function(Experiment, N.sample=1)
{
	sampleRD <- load.file(paste(Experiment@MetaData@DataDirectory, Experiment@MetaData@Instrumental$filename[[N.sample]], sep="/"))
	
	max.rt <- (nrow(sampleRD@data)/(sampleRD@scans.per.second*60)) + sampleRD@start.time/60
	min.rt <- sampleRD@start.time/60

	cat(" Name: \t", as.vector(Experiment@MetaData@Instrumental$filename[[N.sample]]), "\n", "Start Time: \t", min.rt, "\n", "End Time: \t", max.rt, "\n", "Min MZ: \t", sampleRD@min.mz,  "\n", "Max MZ: \t", sampleRD@max.mz, "\n", "Scans/sec: \t", sampleRD@scans.per.second)	
	
}



