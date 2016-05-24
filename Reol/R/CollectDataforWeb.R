PageProcessing <- function(MyEOL, ...) {
  res <- xmlToList(xmlRoot(xmlParse(MyEOL, getDTD=FALSE), ...), simplify=FALSE)
  if(!is.null(res$error)) {
    system(paste("mv", MyEOL, "../TRASH/"))
    stop(paste("Bad file", MyEOL, "has been purged"))
    }
  return(res)
}

FirstTwo <- function(name) {
  if(is.null(name))
    name <- NA
  if(!is.na(name)) {
    name <- trim(name)
    name <- paste(strsplit(name," ")[[1]][1:2],sep=" ",collapse=" ")
  }
  return(name)
}



#common names counts
CNCount <- function(res) {
  whichCNs <- which(names(res$taxonConcept) == "commonName")
  languages <- NULL
  for(j in sequence(length(whichCNs))) {
    languages <- c(languages, as.character(res$taxonConcept[[whichCNs[j]]]$.attr[which(names(res$taxonConcept[[whichCNs[j]]]$.attr)=="lang")])) #not tidy, but effective for multiple entries
  }
  langCounts <- rep(NA,length(unique(languages)))
  for(k in sequence(length(unique(languages)))) {
  	langCounts[k] <- paste(unique(languages)[k], length(which(languages == unique(languages)[k])), sep=".", collapse="")
  }
  return(c(length(unique(languages)), paste(langCounts, collapse="_")))
}

#data object counts
DOCount <- function(res) {
  whichDOs <- which(names(res) == "dataObject")
  dataTypes <- NULL
  IUCNstat <- NA
  for(j in sequence(length(whichDOs))) {
    dataTypes <- c(dataTypes, res[[whichDOs[j]]]$mimeType)
    if(!is.null(res[[whichDOs[j]]]$title)) {
      if(res[[whichDOs[j]]]$title == "IUCNConservationStatus")
        IUCNstat <- res[[whichDOs[j]]]$description
    }
  }
  typeCounts <- rep(NA, length(unique(dataTypes)))
  for(k in sequence(length(unique(dataTypes)))) {
  	typeCounts[k] <- paste(unique(dataTypes)[k], length(which(dataTypes == unique(dataTypes)[k])), sep=".", collapse="")
  }
  return(c(length(unique(dataTypes)), paste(typeCounts, collapse="_"), IUCNstat))
}

#provider numbers
providerCount <- function(res) {
  whichProviders <- which(names(res$taxonConcept$additionalInformation) == "taxon")
  providerTypes <- NULL
  providerIDs <- NULL
  for(j in sequence(length(whichProviders))) {
    providerTypes <- c(providerTypes, res$taxonConcept$additionalInformation[[whichProviders[j]]]$nameAccordingTo)
    providerID <- paste(res$taxonConcept$additionalInformation[[whichProviders[j]]]$nameAccordingTo, res$taxonConcept$additionalInformation[[whichProviders[j]]]$taxonID, sep=".")
    providerIDs <- c(providerIDs, providerID)
  }
  return(c(length(unique(providerTypes)), paste(providerIDs, collapse="_")))
}


#gather vector of information
DataProcessing <- function(res, do.higher.taxonomy) {
  if(!is.null(res$taxonConcept)) {
    taxonData <- c(FirstTwo(res$taxonConcept$ScientificName), res$taxonConcept$ScientificName, res$taxonConcept$taxonConceptID)
    richness <- res$taxonConcept$additionalInformation$richness_score
    refCounts <- length(which(names(res$taxonConcept) == "reference"))
    CNs <- CNCount(res)
    providers <- providerCount(res)
    DOs <- DOCount(res)
    pageLength <- sum(nchar(unlist(res, use.names=FALSE)))
    higher.taxonomy <- ""
    if (do.higher.taxonomy) {
    	try(higher.taxonomy<-paste(MakeTreeData (DownloadHierarchy(res, to.file=FALSE)), collapse="/"))
    	if(nchar(higher.taxonomy)==0) {
    		try(higher.taxonomy<-paste(MakeTreeData (DownloadHierarchy(DownloadEOLpages(res$taxonConcept$taxonConceptID, to.file=FALSE), to.file=FALSE)), collapse="/")) #sometimes the taxon name does not match to the hierarchy, but the eol ID does
    	}	
    }
  }	
  return(matrix(c(taxonData, richness, refCounts, CNs, providers , DOs, pageLength, higher.taxonomy), nrow=1))
}

CollectDataforWeb <- function(MyEOL, do.higher.taxonomy=FALSE) {
#  higher.taxonomy<-""
#  if(do.higher.taxonomy) {
#    try(higher.taxonomy<-paste(MakeTreeData (DownloadHierarchy(MyEOL, to.file=FALSE)), collapse="/"))
#  }
  res <- PageProcessing(MyEOL)
  return(DataProcessing(res, do.higher.taxonomy))
}












