# -- TNRS FUNCTIONS -- #
# Make a function for dealing with names seperate from tnrs functions

ResolveNames <- function(names, max.per.call=100, verbose=TRUE) {
  max.per.call <- 100
  verbose <- FALSE
  # takes a list of names and sends it to the iPlant TNRS site
  # http://tnrs.iplantcollaborative.org/
  # names <- c("zea mays","acacia","solanum","saltea","rosa_rugoso")
  # returnedNames <- ResolveNames(names)
  # print(returnedNames)
  names <- sapply(names, sub, pattern="_", replacement=" ", USE.NAMES=FALSE)
  names <- sapply(names, URLencode, USE.NAMES = FALSE)
  names <- sapply(names, sub, pattern="=", replacement="", USE.NAMES=FALSE)
  names <- sapply(names, sub, pattern="&", replacement="", USE.NAMES=FALSE)
  call.base <- 'http://tnrs.iplantc.org/tnrsm-svc/matchNames?retrieve=best&names='
  new.names <- rep(NA, length(names))
  names.in.call <- 0
  actual.call <- call.base
  starting.position <- 1
  
  for (name.index in sequence(length(names))) {
    names.in.call <- names.in.call + 1
    actual.call <- paste(actual.call, names[name.index], ",", sep="")
    if (names.in.call == max.per.call || name.index == length(names)) {
      returned.values <- suppressWarnings(fromJSON(file=actual.call)$items)
      for (return.index in sequence(length(returned.values)))
        new.names[starting.position + return.index - 
                  1] <- returned.values[[return.index]]$nameScientific
      if (verbose) 
        message(paste("finished ", name.index, "of ", length(names), "names")) 
      starting.position <- name.index + 1
      names.in.call <- 0
      actual.call <- call.base
    }
  }
  if (length(new.names) != length(names)) 
    warning(paste("the input name list was", length(names), 
                  "long but the new one is ", length(new.names), "long"))
  new.names <- sapply(new.names, sub, pattern=" ", replacement="_", USE.NAMES=F)
  return(new.names)
}

GetPhylotasticToken <- function(names, max.per.call=100, verbose=TRUE) {
  max.per.call <- 100
  verbose <- FALSE
  # takes a list of names and sends it to the phylotastic TNRS site
  # http://api.phylotastic.org/tnrs
  # names <- c("zea mays","acacia","solanum","saltea","rosa_rugoso")
  # returns a token that will check names later
  for(i in sequence(length(names))) {
  	if (length(which(strsplit(names[i], split="")[[1]] == " ")) > 1){
  	  #second space in a name is where to cut name sequence
      WhereToCut <- which(strsplit(names[i], split="")[[1]] == " ")[2]-1   
      names[i] <- paste(paste(strsplit(names[i], split="")[[1]][1:WhereToCut], 
                      collapse=""))
  	}
  }
  names <- sapply(names, sub, pattern="_", replacement="+", USE.NAMES=FALSE)
  names <- sapply(names, sub, pattern=" ", replacement="+", USE.NAMES=FALSE)
  names <- sapply(names, URLencode, USE.NAMES = FALSE)
  names <- sapply(names, sub, pattern="=", replacement="", USE.NAMES=FALSE)
  call.base <- "http://www.taxosaurus.org/submit?query="
  new.names <- rep(NA, length(names))
  names.in.call <- 0
  starting.position <- 1
  name.call <- paste(call.base, names[1], sep="")
  for (name.index in 2:length(names)) {
    names.in.call <- names.in.call + 1
    name.call <- paste(name.call, names[name.index], collapse="", sep="%0A")
      if (names.in.call == max.per.call || name.index == length(names)) {
        res <- suppressWarnings(fromJSON(postForm(name.call)))
    }
  }
    cat(res$message, "\n\nOR use the rPlant function RetrieveTNRSNames to pull directly into R")
    return(res$token)
}

RetrieveTNRSNames <- function(names, token, source=c("iPlant_TNRS", "NCBI"), 
                              match.threshold=0.5, verbose=FALSE) {
  web <- "http://www.taxosaurus.org/retrieve"
  res <- suppressWarnings(fromJSON(getURL(paste(web, token, sep="/"))))

  #make results a subsettable matrix
  TNRSnames <- matrix(nrow=length(names), ncol=2)  
  rownames(TNRSnames) <- names 
  colnames(TNRSnames) <- c("TNRS name", "match score") 
  for (i in sequence(length(res$names))) {
    for (j in sequence(length(res$names[[i]]$matches))) {
      if (res$names[[i]]$matches[[j]]$sourceId == source && 
          as.numeric(res$names[[i]]$matches[[j]]$score) > match.threshold) {  
        TNRSnames[which(rownames(TNRSnames) == res$names[[i]]$submittedName),
                  1] <- res$names[[i]]$matches[[j]]$matchedName
        TNRSnames[which(rownames(TNRSnames) == res$names[[i]]$submittedName),
                  2] <- res$names[[i]]$matches[[j]]$score
      }
    }
  }
  #if no TNRS name is returned, then return submitted name
  TNRSnames[which(is.na(TNRSnames[,1])),1] <- rownames(TNRSnames)[which(
            is.na(TNRSnames[,1]))] 
  if (verbose)
    return(TNRSnames)
  return(as.vector(TNRSnames[,1]))
}

CompareNames <- function(old.names, new.names, verbose=TRUE) {
# takes a list of old.names taxonomic names(same ones given as "names in 
# ResolveNames) and compares to the returned names from TNRS
# note that names are changed back to include an "_" instead of the " " they
# come with out of TNRS first, so that they do not count as taxonomic name changes
  taxa.changed <- 0
  comp <- cbind(old.names, new.names)
  for (i in 1: dim(comp)[1]){
    if (comp[i,1] != comp[i,2]) {
      taxa.changed <- taxa.changed + 1
      if (verbose)
        message(paste0(comp[i,1], " was changed to ", comp[i,2], cat("\n")))
    }
  }
  return(paste(taxa.changed, "taxa changed names according to TNRS"))
}
# -- END --#
