"pdb.annotate" <- function(ids, anno.terms=NULL, unique=FALSE, verbose=FALSE) {
  
  oopsa <- requireNamespace("XML", quietly = TRUE)
  oopsb <- requireNamespace("RCurl", quietly = TRUE)
  if(!all(c(oopsa, oopsb)))
    stop("Please install the XML and RCurl package from CRAN")
  
  if(!is.vector(ids)) {
    stop("Input argument 'ids' should be a vector of PDB identifiers/accession codes")
  }
  
  ## All available annotation terms (note 'citation' is a meta term)
  anno.allterms <- c("structureId", "experimentalTechnique", "resolution", "chainId", "ligandId",
                     "ligandName", "source", "scopDomain", "classification", "compound", "title",
                     "citation", "citationAuthor", "journalName", "publicationYear",
                     "structureTitle","depositionDate","structureMolecularWeight","macromoleculeType",
                     "chainId","entityId","sequence","chainLength","db_id","db_name")
                     ##"molecularWeight","secondaryStructure","entityMacromoleculeType")
    
  if(is.null(anno.terms)) {
    anno.terms <- anno.allterms
  } else {
    ## Check and exclude invalid annotation terms
    term.found <- (anno.terms %in% anno.allterms)
    if( any(!term.found) ) {
      warning( paste("Requested annotation term not available:",
                     paste(anno.terms[!term.found], collapse=", "),
                     "\n  Available terms are:\n\t ",
                     paste(anno.allterms, collapse=", ")) )
    }
    anno.terms <- match.arg(anno.terms, anno.allterms, several.ok=TRUE)
  }
  
  ## Check if we have any valid terms remaining
  if( length(anno.terms) == 0 ) {
    stop( paste("No valid anno.terms specified. Please select from:\n\t ",
                paste(anno.allterms, collapse=", ")) )
  }
  
  ## force the structureId and chainId term to be present
  if(any(anno.terms=="citation"))
    req.terms <- c("structureId", "chainId", "ligandId",
                   "citationAuthor", "journalName", "publicationYear")
  else
    req.terms <- c("structureId", "chainId", "ligandId")
  
  anno.terms.input <- anno.terms
  inds <- req.terms %in% anno.terms
  if(!all(inds))
    anno.terms <- c(req.terms[!inds], anno.terms)

  if (missing(ids)) 
    stop("please specify PDB ids for annotating")
  
  if (any(nchar(ids) != 4)) {
    warning("ids should be standard 4 character PDB-IDs: trying first 4 characters...")

    if(unique)
      ids <- unique(substr(basename(ids), 1, 4))

    ## first 4 chars should be upper
    ## any chainId should remain untouched - see e.g. PDB ID 3R1C
    mysplit <- function(x) {
      str <- unlist(strsplit(x, "_"))
      if(length(str)>1)
        paste(toupper(str[1]), "_", str[2], sep="")
      else
        toupper(str[1])
    }
    
    ids <- unlist(lapply(ids, mysplit))
  } else {
    ids <- toupper(ids)
  }
  ids.short <- substr(basename(ids), 1, 4)
  ids.string <- paste(unique(ids.short), collapse=",")
  
  ## prepare query
  query1 = paste(anno.terms, collapse=",")
  url <- paste('http://www.rcsb.org/pdb/rest/customReport')

  curl.opts <- list(httpheader = "Expect:",
                    httpheader = "Accept:text/xml",
                    verbose = verbose,
                    followlocation = TRUE
                    )
  
  curl <- RCurl::postForm(url, pdbids=ids.string, customReportColumns=query1,
                   ssa='n', primaryOnly=1,
                   style = "POST",
                   .opts = curl.opts,
                   .contentEncodeFun=RCurl::curlPercentEncode, .checkParams=TRUE )

  ## parse XML
  xml <- XML::xmlParse(curl)
  data <- XML::xmlToDataFrame(XML::getNodeSet(xml, "/dataset/record"), stringsAsFactors=FALSE)
  
  if(nrow(data)==0)
    stop("Retrieving data from the PDB failed")
  
  ## change colnames (e.g. dimEntity.structureId -> structureId)
  for( i in 1:ncol(data)) {
    a <- unlist(strsplit(colnames(data)[i], "\\."))
    colnames(data)[i] <- a[2]
  }
  
  ## merge data for unique structureId_chainId entries
  if(! (is.null(data$ligandId) & is.null(data$ligandName)) ) {
    if(unique)
      pdbc.ids <- data$structureId
    else
      pdbc.ids <- paste(data$structureId, data$chainId, sep="_")
    unq.ids <- unique(pdbc.ids)
    
    excl.inds <- NULL
    for( i in 1:length(unq.ids) ) {
      inds <- which(pdbc.ids==unq.ids[i])

      if(!is.null(data$ligandId)) {
        tmp.ligId   <- paste(unique(data$ligandId[inds]), collapse=",")
        data$ligandId[inds]   <- tmp.ligId
      }
      if(!is.null(data$ligandName)) {
        tmp.ligName <- paste(unique(data$ligandName[inds]), collapse=",")
        data$ligandName[inds] <- tmp.ligName
      }
      if(!is.null(data$chainId)) {
        tmp.chainId <- paste(unique(data$chainId[inds]), collapse=",")
        data$chainId[inds] <- tmp.chainId
      }
      
      excl.inds <- c(excl.inds, inds[-1])
    }
  }

  if(length(excl.inds)>0)
    data <- data[-excl.inds,]

  if(unique)
    rownames(data) <- data$structureId
  else
    rownames(data) <- paste(data$structureId, data$chainId, sep="_")

  ## include only requested "structureId_chainID"
  row.inds <- unique(unlist(lapply(ids, function(x) grep(x, rownames(data)))))
  data <- data[row.inds,, drop=FALSE]
  
  ## Format citation information
  if (any(anno.terms == "citation") ) {
    citation <- NULL
    lig.auth <- data[,"citationAuthor"]
    lig.year <- data[,"publicationYear"]
    lig.jnal <- data[,"journalName"]
    
    for(i in 1:length(lig.auth)) {
      citation <- c(citation, paste( unlist(strsplit(lig.auth[[i]], ","))[1],
                                    " et al. ", lig.jnal[i], " (", lig.year[i],")",sep=""))
    }
    data <- cbind(data, citation)
  }

  ## include only requested terms
  col.inds <- which(colnames(data) %in% anno.terms.input)
  data <- data[, col.inds, drop=FALSE]

  if(any(data=="null")) {
    inds <- which(data=="null", arr.ind=TRUE)
    for(i in 1:nrow(inds))
      data[ inds[i,"row"], inds[i,"col"] ] = NA
  }

  ## check for missing entries
  mygrep <- function(x, y) {
    inds <- grep(x, y)
    if(length(inds)==0)
      return(NA)
    else
      return(inds)
  }
  
  collected.ids <- rownames(data)
  requested.ids <- ids
  missing <- is.na(unlist(lapply(requested.ids, mygrep, collected.ids)))

  if(any(missing)) {
    missing.str <- paste(requested.ids[which(missing)], collapse=", ")
    
    warning(paste("Annotation data could not be found for PDB ids:\n  ",
                  missing.str))
  }
  
  return(data)
}
