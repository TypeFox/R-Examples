GatherDataObjectInformation <- function(MyEOL) {
  #this function works for one EOL file only.  It will return information about all of the different data objects associated with each taxon.  
  res <- PageProcessing(MyEOL)  
  whichDataObjects <- which(names(res) == "dataObject") 
  NumberOfDataObjects <- length(whichDataObjects) 
  DataObjectInfo <- data.frame(matrix(nrow=NumberOfDataObjects, ncol=1), stringsAsFactors=F)
  scientificName  <- res$taxonConcept[[which(names(res$taxonConcept) == grep("ScientificName", names(res$taxonConcept), ignore.case=TRUE, value=T))]] #because some are cap and some are not
  taxon <-FirstTwo(scientificName)
  if (is.null(taxon)) 
    taxon <- NA
  eolID <- try(res$taxonConcept$taxonConceptID, silent=T)
  if (is.null(eolID)) 
    eolID <- NA
  DataObjectInfo <- data.frame(rep(taxon, NumberOfDataObjects), rep(eolID, NumberOfDataObjects), stringsAsFactors=F)  #initialize dataframe
  colnames(DataObjectInfo) <- c("Taxon", "eolID") 

  #add each data object one by one.  
  for(i in sequence(NumberOfDataObjects)) {
    DO <- unlist(res[[whichDataObjects[i]]])
    for(j in sequence(length(DO))) {
      nameOfColumn <- names(DO)[j]
      if(!any(grepl(paste(nameOfColumn,'*', sep=""), colnames(DataObjectInfo)))) {  #add new column if data doesn't exist
        DataObjectInfo <- cbind(DataObjectInfo, rep(NA, NumberOfDataObjects))
        colnames(DataObjectInfo) <- c(colnames(DataObjectInfo[-length(colnames(DataObjectInfo))]), nameOfColumn) #add new colname
      }
      column <- which(colnames(DataObjectInfo) == nameOfColumn)
      #DataObjectInfo[i,column] <- paste("DO$", nameOfColumn, sep="")
      DataObjectInfo[i,column] <- DO[[j]][1] #Note: sometimes there are multiple elements. Typically, this is is just a main element (like a name of a provider) followed by attributes (URL), but we are just taking the first element regardless.

    }
  }
  if(dim(DataObjectInfo)[1] == 0)
    DataObjectInfo[1,] <- c(taxon, eolID)	
  return(DataObjectInfo)
}


CombineDataObjectInformation <- function(MyEOLs, verbose=TRUE) {
  #Next: subset to Trusted Information only
  #this function works for multiple EOL files.  It will return information about all of the different data objects associated with each taxon.  
  CombinedDOI <- GatherDataObjectInformation(MyEOLs[1])
  for (i in 2:length(MyEOLs)){
    if(verbose)
      print(paste("combined", i, "files"))
    DOI <- GatherDataObjectInformation(MyEOLs[i])
    if(any(!colnames(DOI) %in% colnames(CombinedDOI))) { #check that all new data coming in will match existing data
      ColumnsToAdd <- which(!colnames(DOI) %in% colnames(CombinedDOI))
      for(j in sequence(length(ColumnsToAdd))) {
        CombinedDOI <- cbind(CombinedDOI, rep(NA, dim(CombinedDOI)[1]))
        colnames(CombinedDOI) <- c(colnames(CombinedDOI[-length(colnames(CombinedDOI))]), colnames(DOI)[ColumnsToAdd[j]]) #add new colname
      }
    }
    if(any(!colnames(CombinedDOI) %in% colnames(DOI))) { #check that all new data coming in will match existing data
      ColumnsToAdd <- which(!colnames(CombinedDOI) %in% colnames(DOI))
      for(j in sequence(length(ColumnsToAdd))) {
        DOI <- cbind(DOI, rep(NA, dim(DOI)[1]))
        colnames(DOI) <- c(colnames(DOI[-length(colnames(DOI))]), colnames(CombinedDOI)[ColumnsToAdd[j]]) #add new colname
        #print(paste("added column to DOI", colnames(CombinedDOI)[ColumnsToAdd[j]]))
      }
    }

    ColMatches <- match(colnames(CombinedDOI), colnames(DOI))
    #two ways of sorting dataframes by columns seem to work 1) b[,c(2,1,3)] and 2) subset(b, select=c(2,1,3))
    #rearrange colnames in DOI then add to CombinedDOI
    DOI <- DOI[,ColMatches]
    CombinedDOI <- rbind(CombinedDOI, DOI)
  }
  return(CombinedDOI)
}


DataObjectOverview <- function(MyEOLs, verbose=TRUE){
  MyEOLs <- RemoveNAFiles(MyEOLs)
  if(length(MyEOLs) == 1)
    cDOI <- GatherDataObjectInformation(MyEOLs)

  else
    cDOI <- CombineDataObjectInformation(MyEOLs, verbose=verbose)  
    UniqueTaxa <- unique(cDOI[,which(names(cDOI) == "Taxon")])
    UniqueDataTypes <- unique(cDOI[,which(names(cDOI) == "mimeType")])
    overview <- matrix(nrow=length(UniqueTaxa), ncol=2+length(UniqueDataTypes))
    colnames(overview) <- c("Taxon", "eolID", UniqueDataTypes)
  for(h in sequence(length(UniqueTaxa))) {
    taxonInfo <- NULL
    cDOIsubset <- cDOI[which(cDOI[,1] == UniqueTaxa[h]),]
    taxonName <- as.character(cDOIsubset[1,1])
    eolID <- as.character(cDOIsubset[1,2])
    for(i in sequence(length(UniqueDataTypes))){
      taxonInfo <- append(taxonInfo, length(which(cDOIsubset[,6] == UniqueDataTypes[i])))  
    }
  taxonInfo <- c(taxonName, eolID, taxonInfo)
  overview[h,] <- taxonInfo
  }  
  overview <- data.frame(overview, stringsAsFactors=FALSE)
  if(any(colnames(overview) == "NA.") ||  any(colnames(overview) == "NA")) {
    whichColToDelete <- c(which(colnames(overview)=="NA."), which(colnames(overview)=="NA"))
    overview <- overview[,-whichColToDelete]
  }
  return(overview)
}
