checkUpdates <- function(con,GSEid=NULL){

	gseData <- GSEinDB(con,GSEid)

	dir <- path.package("compendiumdb")

	scriptLoc <- paste(dir,"/scripts/R/downloadSoftFile.R",sep="")
	scriptLoc <- gsub("^","\"",scriptLoc)
	scriptLoc <- gsub("$","\"",scriptLoc) 

	# Download soft file from GEO and retrieve the last updated date
	updateDates <- unlist(lapply(unique(gseData$Experiment),function(x){
				system(paste("Rscript",scriptLoc,x,"gse.new 0",sep=" "))
				data=read.delim("gse.new.soft",sep="=",header=F)
				c(x,as.character(data[data[,1]=="!Series_last_update_date ",2]))
			}))
	updateDates <- matrix(updateDates,ncol=2,byrow=T)
	updateDates <- as.data.frame(updateDates)
	colnames(updateDates)=c("Experiment","LastUpdateDate")

	result <- unique(merge(updateDates,gseData[,c("Experiment","DateLoaded")],by="Experiment"))
	diff <- as.numeric(as.Date(gsub(" .*","",result$DateLoaded),"%Y-%m-%d")-as.Date(result$LastUpdateDate," %b %d %Y"))

	## GSEs updated on GEO
	fRes <- result[which(ifelse(sign(diff)<0,TRUE,FALSE)),]
  if (nrow(fRes)==0){
    cat("None of the specified GSEs loaded in the compendium database have been updated on GEO\n")
  }else{
    rownames(fRes) <- 1:nrow(fRes)
	  fRes
  }
}