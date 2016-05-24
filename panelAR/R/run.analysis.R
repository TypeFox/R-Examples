getRuns <- function(x,times){
	runs <- rle(x)
	runs.length <- runs$lengths[runs$values==TRUE]
	time.end.run <- names(runs$values[runs$values==TRUE])
	time.start.run <- times[which(times %in% time.end.run)-runs.length+1]
	run.table <- data.frame(Start=as.integer(time.start.run),End=as.integer(time.end.run),stringsAsFactors=FALSE)
	return(run.table)
}

run.analysis <- function(object){
	if(class(object)!="panelAR"){
		stop("object must be of class 'panelAR'.",call.=FALSE)
	}
	obs.mat <- object$panelStructure$obs.mat
	times <- colnames(obs.mat)
	units <- rownames(obs.mat)
	
	run.list <- apply(obs.mat,MARGIN=1,function(x) getRuns(x,times=times))
	run.count <- sapply(run.list,function(x) nrow(x))
	run.table <- do.call(rbind,run.list)
	run.table$Length <- run.table$End-run.table$Start+1
	run.table <- as.matrix(run.table)
	rownames(run.table) <- rep(units,run.count)
	

	out <- list(run.count=run.count, runs=run.table,rho=object$panelStructure$rho)
	class(out) <- "panelAR.runs"
	out
	}

print.panelAR.runs <- function(x,...){
	count.table <- data.frame(Unit=names(x$run.count),Runs=x$run.count)
	
	if(any(x$run.count>1) & !is.null(x)){
		names <- paste(names(x$run.count[which(x$run.count>1)]), collapse=", ")
		cat(paste("Calculation of autocorrelation coefficient restarted for each run for the following panels: ",names,"\n\n" ,sep=""))
	}
	cat("Run Counts:\n")
	print(count.table, row.names=FALSE)
}