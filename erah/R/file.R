load.file <- function(filename)
{
	sampleRD <- NULL
	file.extension <- strsplit(filename, split="\\.")[[1]]
	file.type <- file.extension[length(file.extension)]

	if(file.type=="cdf") sampleRD <- load.ncdf4(filename)
	if(file.type=="mzXML" || file.type=="xml") sampleRD <- load.xml(filename)
	if(file.type=="erahrd") sampleRD <- load.erahrd(filename)
	if(is.null(sampleRD)) stop("File extension not recognized. Avalible extensions are: .cdf, .mzXML and .xml")
	
	sampleRD
}

load.erahrd <- function(filename)
{
	sampleRD <- 0
	samp.file <- load(filename)
	sampleObject <- sampleRD
	remove(sampleRD)
	sampleObject
}

load.xml <- function(filename)
{	
	
	if (requireNamespace("mzR", quietly = TRUE)) {
	  	xmlO <- mzR::openMSfile(filename)
		
		metadata <- mzR::runInfo(xmlO)
		
		scans <- metadata$scanCount
		lowMZ <- round(metadata$lowMz+0.5)
		highMZ <- round(metadata$highMz+0.5)
		StartTime <- metadata$dStartTime
		ScansPerSecond <- 1/((metadata$dEndTime - metadata$dStartTime)/metadata$scanCount)
		
		raw.data <- mzR::get3Dmap(object=xmlO, scans=1:scans, lowMz=lowMZ, highMz=highMZ, resMz=1)
		
		sampleRD <- new("RawDataParameters", data = raw.data, min.mz = lowMZ, max.mz = highMZ, start.time = StartTime, mz.resolution = 1, scans.per.second=ScansPerSecond)
		return(sampleRD)

  
   } else {
		msg <- c("mzR is not installed. eRah can operate withouth mzR, unless you want to process .mzXML files (as in this case). To install the mzR package and be able to read mzXML files, please visit its bioconductor website: http://bioconductor.org/packages/release/bioc/html/mzR.html
Or, alternatively, execute the following R code:
		
		## try http:// if https:// URLs are not supported 
		source('https://bioconductor.org/biocLite.R')
		biocLite('mzR')")
		    
		warning(msg)   

   }
}


load.ncdf4 <- function(filename)
{	
	
	measurement = nc_open(filename)	
	
	mass_values <- ncvar_get(measurement,"mass_values")   
	mass_intensities <- ncvar_get(measurement,"intensity_values") 
	scan_indexes <- ncvar_get(measurement,"scan_index")
	min_mz <- as.numeric(ncvar_get(measurement,"mass_range_min", count=1))
	max_mz <- as.numeric(ncvar_get(measurement,"mass_range_max", count=1))
	start_time <- as.numeric(ncvar_get(measurement,"scan_acquisition_time", count=1))
	scans_per_second <- as.numeric((1/ncvar_get(measurement,"scan_duration", count=1)))
				
	full.matrix <- matrix(0,length(scan_indexes),((max_mz-min_mz)+1))
	mass_values <- mass_values - (min_mz - 1)
	
	for(i in 1:(length(scan_indexes)-1))
	{
		full.matrix[i,mass_values[(scan_indexes[i]+1):scan_indexes[i+1]]] <- mass_intensities[(scan_indexes[i]+1):scan_indexes[i+1]]
	}
	
	sampleRD <- new("RawDataParameters", data = full.matrix, min.mz = min_mz, max.mz = max_mz, start.time = start_time, mz.resolution = 1, scans.per.second=scans_per_second)
	sampleRD
}


createdt <- function(path)
{
	#path <- "Valli/GenCond"
	#path.dir <- list.files(path)
	path.name <- strsplit(path, "/")[[1]]
	path.name <- path.name[length(path.name)]

	#dirs.c <- unlist(apply(as.matrix(path.dir), 1, function(x) rep(x,length(list.files(paste(path, x, sep="/"))))))
	#path.dir.c <- apply(as.matrix(path.dir),1, function(x) paste(path, x, sep="/"))
	#files.c <- list.files(path.dir.c)
	
	#files.name <- apply(as.matrix(1:length(dirs.c)),1, function(x) paste(dirs.c[x],files.c[x], sep="/"))
	#files.name <- list.files(path.dir.c, full.name=T)
	#files.ID <- apply(as.matrix(files.c), 1, function(x) strsplit(x, "\\.")[[1]][1])

	files.name <- list.files(path, recursive=T)
	
	if(any(apply(as.matrix(files.name), 1, function(x) length(strsplit(x, "/")[[1]]))==1)) stop("There are files without directory in the selected path. Remove all the files in the path, only folders are allowed")
	
	files.class <- apply(as.matrix(files.name), 1, function(x) strsplit(x, "/")[[1]][1])
	files.ID <- apply(as.matrix(files.name), 1, function(x) {
		out.s <- strsplit(x, "/")[[1]]
		strsplit(out.s[length(out.s)], "\\.")[[1]][1] 
		})

	files.path <- apply(as.matrix(files.name),1, function(x) paste(path, x, sep="/"))	
	files.cdate <- apply(as.matrix(files.path), 1, function(x) as.character(file.info(x)$mtime))
	files.date <- apply(as.matrix(files.cdate),1, function(x) strsplit(x, " ")[[1]][1])
	files.time <- apply(as.matrix(files.cdate),1, function(x) strsplit(x, " ")[[1]][2])
		
	inst.table <- matrix(0, ncol=4, nrow=length(files.ID))
	colnames(inst.table) <- c("sampleID", "filename", "date", "time") 
	
	inst.table[,1] <- files.ID
	inst.table[,2] <- files.name
	inst.table[,3] <- files.date
	inst.table[,4] <- files.time
	
	meta.table <- matrix(0, ncol=2, nrow=length(files.ID))
	colnames(meta.table) <- c("sampleID", "class") 
	
	meta.table[,1] <- files.ID
	meta.table[,2] <- files.class
	
	inst.file <- paste(path, "/", path.name, "_inst.csv", sep="")
	meta.file <- paste(path, "/", path.name, "_pheno.csv", sep="")

	write.table(inst.table, file=inst.file, sep=";", row.names=FALSE, eol="\n", quote=F)
	write.table(meta.table, file=meta.file, sep=";", row.names=FALSE, eol="\n", quote=F)	

}

