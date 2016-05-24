setClass(
		Class = "NG_Visualization",
		representation = representation(
				graph = "character",     # graph name
				data = "character",      # data name
				from = "character",      # bullet comes from node
				to = "character",        # bullet goes to node
				varList = "character"    # list all variables in the data
		),
		contains = "VIRTUAL"
)


## show method
setMethod(f = "show",
		signature = "NG_Visualization",
		definition = function(object){
			if(is(object,"NG_Viz2DAxis")){
				cat("2D Axis Plot: nd_2d_myplot()\n")
				cat(paste("Plot function:",object@FUN,'\n'))
			}else if(is(object,"NG_Viztk2d")){
				cat("tk2d scatter plot: ng_2d()\n")
			}else {
				cat("Custom plot\n")
			}
			cat(paste("Graph:", object@graph,'\n'))
			cat(paste("Data:", object@data,'\n'))
		}
)


## check whether the names in the graph nodes come from the data
## names or shortnames
vizVarNames <- function(graph,data) {
	parsedNames <- unlist(strsplit(nodes(graph@graph),graph@sep, fixed = TRUE ))
	isLong <- all(parsedNames %in% names(data@data))
	
	if(length(shortnames(data))>0) {
		isShort <- all(parsedNames %in% data@shortnames)
	}else {
		isShort <- FALSE
	}
	
	if(isShort && !isLong) {
		varNames <- shortnames(data)
	}else if(isLong && !isShort) {
		varNames <- names(data)
	}else if (all(names(data) == shortnames(data))) {
		varNames <- names(data)
	}else {
		stop("[ng.sp] Can not figure out whether to associate shortnames or names with the graph nodes.")
	}
	
	## check if sep is consistant
	if(any(sapply(varNames,function(name)grepl(graph@sep, name, fixed = TRUE)))){
		stop(paste("[ng_sp] your NG_data object has shortnames or names which include sep (",graph@sep,")\n",sep = ""))
	}
	
	return(varNames)
}

## initialize rotation matrix mat and axis history matrix axisMat
initRotation <- function(viz,ngEnv) {
	viz@mat <- viz@mat*0
	ifrom <- match(unlist(strsplit(viz@from, split = ngEnv$graph@sep, fixed = TRUE)), viz@varList)
	
	viz@mat[ifrom[1],1] <- 1 ## x-axis
	viz@mat[ifrom[2],2] <- 1 ## y-axis
	
	## Reset History Matrix
	viz@axisMat <- matrix(rep(0,4),ncol = 2)
	viz@axisMat[1,1] <- ifrom[1]
	viz@axisMat[2,1] <- ifrom[2]
	
	ngEnv$bulletState$to <- ""
	
	return(viz)
}

