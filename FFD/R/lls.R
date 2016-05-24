## Ian Kopacka
## 2011-05-31
##
## Function: lls
## 
## Function works like ls but it returns the names of the objects 
## in the workspace, together with class, dimension and size in the
## form of a data frame.

lls <- function(name, pos = -1, envir = as.environment(pos), all.names = FALSE, 
		pattern, classFilter, sort = "size")
{ 
	if( !missing(name) )
	{
		envir <- name
	} else
	{
		envir <- parent.frame()
	}	
	lsVec <- ls(name, pos, envir, all.names, pattern)  
	classVec <- character(length(lsVec))
	sizeVec <- character(length(lsVec))
	dimVec <- character(length(lsVec))
	
	for(ii in seq(along = lsVec))
	{
		item <- lsVec[ii]
		realItem <- eval(parse(text = item), envir)
		classItem <- class(realItem)
		if (length(classItem) > 0) classItem <- Reduce(function(x,y) paste(x,y,sep = ", "), classItem)
		classVec[ii] <- classItem
		sizeVec[ii] <- as.numeric(object.size(realItem))
		dimVec[ii] <- paste(dim(realItem), collapse = "x" )
		if( dimVec[ii] == "" )
		{
			dimVec[ii] <- length(realItem)
		}
	} 
	sizeVec <- as.numeric(sizeVec)
	sizeUnits <- ifelse(sizeVec < 1024, paste(sizeVec, "B"), 
		ifelse(sizeVec < 1024*1024, paste(round(sizeVec/1024, 1), "KB"), 
				paste(round(sizeVec/1024/1024, 1 ),"MB")))
	lsDf <- data.frame(Name = lsVec, Class = classVec, Dimension = dimVec, 
		Size_Bytes = sizeVec, Size_Unit = sizeUnits)
	if (sort == "size") lsDf <- lsDf[order(lsDf$Size_Bytes, decreasing = TRUE),]
	if (!missing(classFilter)) lsDf <- lsDf[lsDf$Class == classFilter,]
	rownames(lsDf) <- NULL
	lsDf
}
