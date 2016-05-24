generate.names <-
function(original.names,regex.codes,standard.names=NULL,set.na=FALSE,suggest=FALSE,verbose=FALSE) {
  
  #Get both classified and unclassified locations returned
  obtained.values <- find.names(original.names,regex.codes,standard.names,set.na=set.na,suggest=suggest,verbose=verbose)
  original.locations <- obtained.values[[1]]
  unclassified.values <- obtained.values[[2]]
  suggested.values <- unclassified.values
  #Order names according to original order
  ordered.names <- original.names[original.locations]
  #Get all regex codes
  full.codes <- regex.codes[,2]
  
  #If no standard names provided, will consider the input the standard naming system
  #Will return the standard names and corresponding codes
  if(dim(unclassified.values)[1] > 0) {
    unclassified.values <- unclassified.values[,1:2]
    suggested.values <- na.omit(suggested.values[,2:3])
  }
  if (!is.null(unclassified.values)) {
    full.codes <- c(full.codes,as.character(unclassified.values[,1]))
  }
  if (is.null(standard.names)){
    return(data.frame(code=full.codes,name=ordered.names))
  }
  else if (!is.null(unclassified.values)) {
    standard.names <- rbind(standard.names,unclassified.values)
  }
  
  #Find locations in standard names which correspond to input names
  classified.locations <-!is.na(ordered.names)
  #Get ordered regex codes for these locations
  standard.codes <- full.codes[classified.locations]
  if (sum(is.na(standard.names))>0) {
    fallback.standard.names <- generate.names(original.names,regex.codes)
  }
  fallback.names.used <- NULL
  
  #Function to convert codes to standard names
  obtain.standard <- function(x) {
    loop.location <- which(standard.names[,1]==x)
    loop.name <- as.character(standard.names[loop.location,2])
    if (is.na(loop.name)) {
      loop.location <- which(fallback.standard.names[,1]==x)
      loop.name <- as.character(fallback.standard.names[loop.location,2])
      fallback.names.used <<- c(fallback.names.used,loop.name)
      if(set.na) {
        loop.name <- NA
      }
    }
    return(loop.name)
  }
  
  #Change codes to standardized names
  complete.standard <- as.vector(sapply(standard.codes,obtain.standard))
  #Print unrecognized names
  if (!is.null(fallback.names.used) && length(fallback.names.used)>0) {
    if (set.na) {
      action.taken <- "removed"
    } else {
      action.taken <- "left unchanged"
    }
    if(verbose || set.na) {
      cat(sprintf("\nThe following names were not in the standard set and %s:\n",action.taken))
      print(fallback.names.used)
    } else {
      cat(sprintf("\nNote: %d names were not in the standard set and %s.\n",length(fallback.names.used),action.taken))
    }
  }
  
  #Original ordering of the names
  original.locations.nona <- as.numeric(na.omit(original.locations))
  #Order standardized names by orginal ordering
  output.names <- alternate.order(complete.standard,original.locations.nona)
  
  #Returns
  if (suggest && dim(suggested.values)[1] > 0) {
    #Change suggestion codes to standard names and return with ordered standard names
    suggestions.standard <- as.vector(sapply(as.character(suggested.values[,2]),obtain.standard))
    suggested.values <- data.frame(Original=suggested.values[,1],Suggested=suggestions.standard)
    return(list(output.names,suggested.values))
  }else{
    #Return ordered standardized names
    return(output.names)
  }
}
