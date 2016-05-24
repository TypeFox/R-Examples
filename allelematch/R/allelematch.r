## allelematch R Package
## v2.0
## allelematch:  Pairwise matching and identification of unique multilocus genotypes
##
## by Paul Galpern
## License: GPL-2
##
## Last update:  12 September 2011
##
## N.B. Renamed from MicroSatMatch as of v1.1
##
## Functions:
## amDataset()             Produces an input dataset object for allelematch routines
## print.amDataset()       Print method for amDataset objects
## amMatrix()              Produce a dissimilarity matrix
## amPairwise()            Pairwise matching of genotypes
## summary.amPairwise()    Summary method for amPairwise objects
## amCluster()             Clustering of genotypes  
## summary.amCluster()     Summary method for amCluster objects
## amUnique()              Identification of unique genotypes
## summary.amUnique()      Summary method for amUnique objects
## amUniqueProfile()       Utility to find optimal parameters for amUnique()
## amAlleleFreq()          Produces allele frequencies from an amDataset object
##
## Requires:  dynamicTreeCut
##
## Please see R documentation for full description of function parameters
##
## amDataset()
amDataset <- function(multilocusDataset, missingCode="-99", indexColumn=NULL, metaDataColumn=NULL, ignoreColumn=NULL) {
    
                
    ## Create amDataset object
    newDataset <- list()
    class(newDataset) <- "amDataset"

    ## Check function call variables for validity
    if (is.null(dim(multilocusDataset))) stop("allelematch:  multilocusDataset must be a matrix or a data.frame", call.=FALSE)
    if ((is.character(indexColumn) || is.character(metaDataColumn) || is.character(ignoreColumn)) && is.null(dimnames(multilocusDataset)[[2]])) {
        stop("allelematch:  multilocusDataset does not have dimnames()[[2]] set; use integer column indices or set dimnames()", call.=FALSE)
    }
    if (!is.null(indexColumn)) {
        if (length(indexColumn) > 1) stop("allelematch:  only one indexColumn permitted", call.=FALSE)
        if (is.character(indexColumn)) {
            indexColumnWhich <- which(indexColumn == dimnames(multilocusDataset)[[2]])
            if (length(indexColumnWhich) == 0) stop("allelematch:  indexColumn does not exist in multilocusDataset", call.=FALSE)
        }
        else {
            indexColumnWhich <- indexColumn
            if (length(indexColumn))
            if (indexColumnWhich > ncol(multilocusDataset)) stop("allelematch:  indexColumn does not exist in multilocusDataset", call.=FALSE)
        }
    }
    else {
        indexColumnWhich <- 0
    }
    if (!is.null(metaDataColumn)) {
        if (length(metaDataColumn) > 1) stop("allelematch:  only one metaDataColumn permitted", call.=FALSE)
        if (is.character(metaDataColumn)) {
            metaDataColumnWhich <- which(metaDataColumn == dimnames(multilocusDataset)[[2]])
            if (length(metaDataColumnWhich) == 0) stop("allelematch:  metaDataColumn does not exist in multilocusDataset", call.=FALSE)
        }
        else {
            metaDataColumnWhich <- metaDataColumn
            if (length(metaDataColumn))
            if (metaDataColumnWhich > ncol(multilocusDataset)) stop("allelematch:  metaDataColumn does not exist in multilocusDataset", call.=FALSE)
        }
    }
    else {
        metaDataColumnWhich <- 0
    }
    if (!is.null(ignoreColumn)) {
        if (is.character(ignoreColumn)) {
            ignoreColumnWhich <- which(as.logical(rowSums(sapply(ignoreColumn, function(x) x == dimnames(multilocusDataset)[[2]]))))
            if (length(ignoreColumnWhich) != length(ignoreColumn)) stop("allelematch:  one or more ignoreColumn does not exist in multilocusDataset", call.=FALSE)
        }
        else {
            ignoreColumnWhich <- ignoreColumn
            if (length(ignoreColumn))
            if (any(ignoreColumnWhich > ncol(multilocusDataset))) stop("allelematch:  one or more ignoreColumn does not exist in multilocusDataset", call.=FALSE)
        }
    }
    else {
        ignoreColumnWhich <- 0
    }
    
    
    
    ## Prepare multilocusDataset
    columnDataset <- dimnames(multilocusDataset)[[2]]
    multilocusDataset <- t(apply(multilocusDataset, 1, as.character))

    ## Change NA data to the missingCode
    if (sum(is.na(multilocusDataset)) > 0) {
        multilocusDataset[is.na(multilocusDataset)] <- missingCode
        cat("allelematch:  NA data converted to", missingCode, "\n")
    }
    newDataset$index <- multilocusDataset[, indexColumnWhich]
    newDataset$metaData <- multilocusDataset[, metaDataColumnWhich]
    keepTheseColumns <- 1:ncol(multilocusDataset) %in% c(indexColumnWhich, metaDataColumnWhich, ignoreColumnWhich)
    if (sum(!keepTheseColumns) < 3) stop("allelematch:  at least three data columns are required for allelematch", call.=FALSE)
    newDataset$multilocus <- multilocusDataset[, !keepTheseColumns]
    ## Remove spaces from the multilocus and index columns (overcomes problems caused by space padding that may crop up)
    newDataset$multilocus <- t(apply(newDataset$multilocus, 1, function(x) gsub(" ", "", x)))
    newDataset$index <- gsub(" ", "", newDataset$index)
    columnDataset <- columnDataset[!keepTheseColumns]

    
   
    ## Assign index and multilocus column names if not given
    if (length(newDataset$index)==0) {
        if (nrow(newDataset$multilocus) > 17576) {
            ## N.B. 17576 is length of labelRepository (below)
            stop("allelematch:  too many samples for automatic assignment of index;  please provide an index column", call.=FALSE)
        }
        ## Produce a vector of possible index or column name labels 
        labelRepository <- as.character(apply(cbind(rep(1:26, each=26*26),
                                        do.call(rbind,lapply(1:26, function(x) cbind(rep(1:26, each=26),rep(1:26, times=26))))), 1, function(x)
                                        paste(LETTERS[x[1]], LETTERS[x[2]], LETTERS[x[3]], sep="")))
        
        newDataset$index <- labelRepository[1:nrow(newDataset$multilocus)]
    }
    if (length(newDataset$metaData)==0) {
        newDataset$metaData <- NULL
    }
    if (is.null(columnDataset)) {
        dimnames(newDataset$multilocus)[[2]] <- paste("loc", 1:ncol(newDataset$multilocus), sep="") 
    }
    else {
        dimnames(newDataset$multilocus)[[2]] <- columnDataset
    }
    
    if (sum(duplicated(newDataset$index)) > 0) stop("allelematch:  index column should contain a unique identifier for each sample", call.=FALSE)
    
    if (is.na(missingCode)) {
        newDataset$multilocus[is.na(newDataset$multilocus)] <- "NA"
        newDataset$missingCode <- "NA"
    }
    else {
        newDataset$missingCode <- missingCode
    }
    
    return(newDataset)
}
##
## print.amDataset()
print.amDataset <- function(x, ...) {
    
    cat("allelematch\namDataset object\n")
    if (!is.null(x$metaData)) {
        xPretty <- paste(format(x$index), format(x$metaData), sep="  ")    
    }
    else {
        xPretty <- format(x$index)
    }
    
    xPrint <- data.frame(x$multilocus, row.names=xPretty)
    
   print(xPrint)
}

##
## amMatrix()
amMatrix <- function(amDatasetFocal, missingMethod=2) {
  
  ## Check function call variables for validity
  if (class(amDatasetFocal) != "amDataset") stop("allelematch:  amDatasetFocal must be an object of class \"amDataset\";  use amDataset() first", call.=FALSE)
  if (!(missingMethod %in% c(1,2))) stop("allelematch:  missingMethod must equal 1 or 2", call.=FALSE)
  
  ## Create variables from amDatasetFocal object
  numGenotypes <- nrow(amDatasetFocal$multilocus)
  genotypes <- amDatasetFocal$multilocus
  ## Set missing data to NA for convenience
  genotypes[genotypes==amDatasetFocal$missingCode] <- NA
  
  
  ## Empty data structure to store results
  simMatrix <- matrix(, nrow=numGenotypes, ncol=numGenotypes)

  ## Determine allele similarity score
  for (i in 1:numGenotypes) {
    if (missingMethod > 0) {
        
        if (missingMethod==1) missingMultiplier <- 0.25
        else missingMultiplier <- 0.5
    
        ## Treat missing data in either the focal or comparison genotype as a partial match except:
        ##      When missingMethod=2 missing data matches perfectly with missing data
        ##      When missingMethod=1 missing data matches partially with missing data
        simMatrix[i,] <- as.double(rowSums(genotypes[rep(i, numGenotypes),]==genotypes, na.rm=TRUE) +
                                                          rowSums(is.na(genotypes) * missingMultiplier) + sum(is.na(genotypes[i,]) * missingMultiplier))

      }
    else {
        ## Treat missing data in the focal genotype, or missing data in the comparison genotype, or missing data matching in both as a match
        ##simMatrix[i,] <- as.double(rowSums(genotypes[rep(i, numGenotypes),]==genotypes, na.rm=TRUE) +
        ##                                                  rowSums(is.na(genotypes) | is.na(genotypes[i,])))
        stop("allelematch:  missingMethod=0 is currently not implemented", call.=FALSE)
      }
  }
  
  ## Turn matches into a percent and make it a dissimilarity
  dissimMatrix <- 1-(simMatrix/ncol(genotypes))
  
  ## Add labels to dissimilarity matrix
  dimnames(dissimMatrix) <- list(amDatasetFocal$index, amDatasetFocal$index)
  
  ## Note that the matrix produced has a non-zero diagonal when missingMethod=0 or 1
  ## This doesn't affect analyses downstream
  class(dissimMatrix) <- "amMatrix"
  return(dissimMatrix)
}
##
## amPairwise()
amPairwise <- function(amDatasetFocal, amDatasetComparison=amDatasetFocal, alleleMismatch=NULL, matchThreshold=NULL, missingMethod=2) {
    
    ## Check function call variables for validity
    if ((class(amDatasetFocal) != "amDataset") || (class(amDatasetComparison) != "amDataset")) {
          stop("allelematch:  amDatasetFocal and amDatasetComparison must be an object of class \"amDataset\";  use amDataset() first", call.=FALSE)
    }
    if (!(missingMethod %in% c(1,2))) stop("allelematch:  missingMethod must equal 1 or 2", call.=FALSE)
    
    ## More checking of input parameters
    if (sum(!(c(is.null(alleleMismatch), is.null(matchThreshold)))) != 1) {
        stop("allelematch:  please specify alleleMismatch OR matchThreshold.", call.=FALSE)
    }
    if (length(c(alleleMismatch, matchThreshold))>1) {
            stop("allelematch:  please provide a single parameter value for alleleMismatch OR matchThreshold.", call.=FALSE)
        }
    
    if (!is.null(matchThreshold)) {
        if ((matchThreshold < 0) || (matchThreshold > 1)) {
                stop("allelematch:  matchThreshold must be between 0 and 1", call.=FALSE)
            }
        alleleMismatch <- round((1-matchThreshold)*ncol(amDatasetFocal$multilocus),2)
    }
    else if (!is.null(alleleMismatch)) {
        matchThreshold <- 1-(alleleMismatch/ncol(amDatasetFocal$multilocus))
        }
  
    
    ## Put amDataset object into convenience variables
    focalGenotypes <- amDatasetFocal$multilocus
    comparisonGenotypes <- amDatasetComparison$multilocus
    indexFocal <- amDatasetFocal$index
    indexComparison <- amDatasetComparison$index
    metaDataFocal <- amDatasetFocal$metaData
    metaDataComparison <- amDatasetComparison$metaData
    columnNames <- dimnames(amDatasetFocal$multilocus)[[2]]
    
    ## Check focal and comparison datasets have the same number of loci, column names, and missing cdodes
    if (ncol(focalGenotypes) != ncol(comparisonGenotypes))  {
      stop("allelematch:  amDatasetFocal and amDatasetComparison must have the same number of columns / loci", call.=FALSE)
    }
    if (any(dimnames(amDatasetFocal$multilocus)[[2]] != dimnames(amDatasetComparison$multilocus)[[2]])) {
      stop("allelematch:  amDatasetFocal and amDatasetComparison must have columns with the same names", call.=FALSE)
    }
    if (amDatasetFocal$missingCode != amDatasetComparison$missingCode) {
      stop("allelematch:  amDatasetFocal and amDatasetComparison must have same missingCode", call.=FALSE)
    }
    
    ## Set missingCodes to NA
    focalGenotypes[focalGenotypes==amDatasetFocal$missingCode] <- NA
    comparisonGenotypes[comparisonGenotypes==amDatasetComparison$missingCode] <- NA
     
    numFocalGenotypes <- nrow(focalGenotypes)
    numComparisonGenotypes <- nrow(comparisonGenotypes)
    
    
    ## Empty data structures to store results
    simMatrix <- matrix(, nrow=numFocalGenotypes, ncol=numComparisonGenotypes)
    pairwiseMatches <- vector("list", numFocalGenotypes)
  
    ## Determine allele similarity score
    for (i in 1:numFocalGenotypes) {
      if (missingMethod > 0) {
          if (missingMethod==1) missingMultiplier = 0.25
          else missingMultiplier <- 0.5
      
          ## Treat missing data in either the focal or comparison genotype as a partial match except:
          ##      When missingMethod=2 missing data matches perfectly with missing data
          ##      When missingMethod=1 missing data matches partially with missing data
          simMatrix[i,] <- as.double(rowSums(focalGenotypes[rep(i, numComparisonGenotypes),]==comparisonGenotypes, na.rm=TRUE) +
                                                          rowSums(is.na(comparisonGenotypes) * missingMultiplier) + sum(is.na(focalGenotypes[i,]) * missingMultiplier))
      }
      else {
          ## Treat missing data in the focal genotype, or missing data in the comparison genotype, or missing data matching in both as a match
          ##simMatrix[i,] <- as.double(rowSums(focalGenotypes[rep(i, numComparisonGenotypes),]==comparisonGenotypes, na.rm=TRUE) +
          ##                                                rowSums(is.na(comparisonGenotypes) | is.na(focalGenotypes[i,])))
          stop("allelematch:  missingMethod=0 is currently not implemented", call.=FALSE)
      }
      
      ## Determine which comparison genotypes meet the threshold
      simMatrix <- simMatrix/ncol(focalGenotypes)
      pairwiseMatchesWhich <- which(simMatrix[i, ]  >= matchThreshold)
      pairwiseMatchesScores <- signif(simMatrix[i, pairwiseMatchesWhich],2)
      
      ## Focal genotype
      pairwiseMatches[[i]]$focal <- list(index=indexFocal[i], metaData=metaDataFocal[i], multilocus=t(focalGenotypes[i,]))
      ## Set column names on multilocus
      dimnames(pairwiseMatches[[i]]$focal$multilocus)[[2]] <- columnNames
      ## Flag missing genotypes
      pairwiseMatches[[i]]$focal$flags <- matrix(as.integer(is.na(pairwiseMatches[[i]]$focal$multilocus))+1, 1, ncol=length(columnNames))
      ## Return missing codes to multilocus
      pairwiseMatches[[i]]$focal$multilocus[is.na(pairwiseMatches[[i]]$focal$multilocus)] <- amDatasetFocal$missingCode
      
      ## Comparison or "matching" genotype
      pairwiseMatches[[i]]$match <- list(index=NULL, metaData=NULL, multilocus=NULL, score=NULL)
      ## No matches
      if (length(pairwiseMatchesWhich) == 0) {
          pairwiseMatches[[i]]$match$index <- "None"
          if (!is.null(metaDataFocal[i])) {
              pairwiseMatches[[i]]$match$metaData <- ""
          }
          else {
              pairwiseMatches[[i]]$match$metaData <- NULL
          }
          pairwiseMatches[[i]]$match$multilocus <- matrix(, 1, length(columnNames))      
          pairwiseMatches[[i]]$match$multilocus[1, ] <- rep("", length(columnNames))
          pairwiseMatches[[i]]$match$score <- ""
          pairwiseMatches[[i]]$match$flags <- matrix(1, nrow=1, ncol=length(columnNames))
          pairwiseMatches[[i]]$match$perfect <- 0
          pairwiseMatches[[i]]$match$partial <- 0
      }
      ## One or more matches
      else {
          pairwiseMatches[[i]]$match$multilocus <- matrix(, length(pairwiseMatchesWhich), length(columnNames))      
          for (j in 1:length(pairwiseMatchesWhich)) {
              pairwiseMatches[[i]]$match$index[j] <- indexComparison[pairwiseMatchesWhich[j]]
              pairwiseMatches[[i]]$match$metaData[j] <- metaDataComparison[pairwiseMatchesWhich[j]]
              pairwiseMatches[[i]]$match$multilocus[j, ] <- comparisonGenotypes[pairwiseMatchesWhich[j], ]
              pairwiseMatches[[i]]$match$score[j] <- format(signif(pairwiseMatchesScores[j],2))
          }
          ## Sort in decreasing order if two or more matches
          if (nrow(pairwiseMatches[[i]]$match$multilocus) > 1) {
              pairwiseMatches[[i]]$match$index <- pairwiseMatches[[i]]$match$index[order(pairwiseMatchesScores, decreasing=TRUE)]
              pairwiseMatches[[i]]$match$metaData <- pairwiseMatches[[i]]$match$metaData[order(pairwiseMatchesScores, decreasing=TRUE)]
              pairwiseMatches[[i]]$match$multilocus <- pairwiseMatches[[i]]$match$multilocus[order(pairwiseMatchesScores, decreasing=TRUE),]
              pairwiseMatches[[i]]$match$score <- pairwiseMatches[[i]]$match$score[order(pairwiseMatchesScores, decreasing=TRUE)]
          }
          ## Create match flags for display purposes
          ## 0 = alleles do not match
          ## 1 = alleles match
          ## 2 = alleles are missing or alleles are missing in focal but present in comparison
          
          ## Set all flags to matching
          pairwiseMatches[[i]]$match$flags <- matrix(1, nrow(pairwiseMatches[[i]]$match$multilocus), ncol(pairwiseMatches[[i]]$match$multilocus))
          ## Set flags that mismatch
          pairwiseMatches[[i]]$match$flags[matrix(pairwiseMatches[[i]]$focal$multilocus,
                                                 nrow(pairwiseMatches[[i]]$match$multilocus), ncol(pairwiseMatches[[i]]$match$multilocus), byrow=TRUE)
                                          != pairwiseMatches[[i]]$match$multilocus] <- 0
          ## Set flags that are missing
          pairwiseMatches[[i]]$match$flags[is.na(pairwiseMatches[[i]]$match$multilocus)] <- 2
          pairwiseMatches[[i]]$match$flags[, which(pairwiseMatches[[i]]$focal$flags==2)] <- 2
          
          ## Summarize the matches
          pairwiseMatches[[i]]$match$perfect <- sum(pairwiseMatchesScores==1)
          pairwiseMatches[[i]]$match$partial <- sum(pairwiseMatchesScores<1)
      }
      
      ## Set column names on multilocus
      dimnames(pairwiseMatches[[i]]$match$multilocus)[[2]] <- columnNames
      
      ## Return missing codes to multilocus
      pairwiseMatches[[i]]$match$multilocus[is.na(pairwiseMatches[[i]]$match$multilocus)] <- amDatasetFocal$missingCode
      }
    
    ## Add some meta-data to the amPairwise object
    amPairwise <- list()
    amPairwise$pairwise <- pairwiseMatches
    amPairwise$missingCode <- amDatasetFocal$missingCode
    amPairwise$matchThreshold <- matchThreshold
    amPairwise$alleleMismatch <- alleleMismatch
    amPairwise$missingMethod <- missingMethod
    amPairwise$focalDatasetN <- nrow(amDatasetFocal$multilocus)
    amPairwise$comparisonDatasetN <- nrow(amDatasetComparison$multilocus)
    if ((dim(amDatasetFocal$multilocus)==dim(amDatasetComparison$multilocus)) && (all(amDatasetFocal$multilocus==amDatasetComparison$multilocus))) {
      amPairwise$focalIsComparison <- TRUE
    }
    else {
      amPairwise$focalIsComparison <- FALSE
    }
  
    class(amPairwise) <- "amPairwise"
    return(amPairwise) 
  }


##
## summary.amPairwise()
summary.amPairwise <- function(object, html=NULL, csv=NULL, ...) {
    
    if (class(object) != "amPairwise") {
        stop("allelematch:  this function requires an \"amPairwise\" object", call.=FALSE)
    }
    
    if (!is.null(html)) {
        if (is.logical(html) && (html==TRUE)) amHTML.amPairwise(object)
        else if (is.logical(html) && (html==FALSE)) html <- NULL
        else amHTML.amPairwise(object, htmlFile=html)
    }
    if (!is.null(csv)) {
        amCSV.amPairwise(object, csvFile=csv)
    }
    if (is.null(html) && is.null(csv)) {
        cat("allelematch\npairwise analysis\n")
        cat("\nfocal dataset N=", object$focalDatasetN, "\n", sep="")
        if (object$focalIsComparison) cat("focal dataset compared against itself\n")
        else cat("comparison dataset N=", object$comparisonDatasetN, "\n", sep="")
        cat("missing data represented by: ", object$missingCode, "\n", sep="")
        cat("missing data matching method: ", object$missingMethod, "\n", sep="")
        cat("alleleMismatch (m-hat; maximum number of mismatching alleles): ", object$alleleMismatch, "\n", sep="")
        cat("matchThreshold (s-hat; lowest matching score returned): ", object$matchThreshold, "\n\n", sep="")
        cat("score flags:\n")
        cat("*101 allele does not match\n")
        cat("+101 allele is missing\n\n")
        y <- object$pairwise
        for (i in 1:length(y)) {
            cat("(", i, " of ", length(y), ")\n", sep="")
            if (is.null(y[[i]]$focal$metaData)) y[[i]]$focal$metaData <- ""
            if (is.null(y[[i]]$match$metaData)) y[[i]]$match$metaData <- rep("", length(y[[i]]$match$index))
            
            showFocalFlags <- y[[i]]$focal$multilocus
            showFocalFlags[y[[i]]$focal$flags==2] <- "+"
            showFocalFlags[y[[i]]$focal$flags==0] <- "*"
            showFocalFlags[y[[i]]$focal$flags==1] <- ""
            showFocal <- matrix(paste(showFocalFlags, y[[i]]$focal$multilocus, sep=""), nrow(showFocalFlags), ncol(showFocalFlags))
            dimnames(showFocal) <- dimnames(showFocalFlags)
            
            showMatchFlags <- y[[i]]$match$multilocus
            showMatchFlags[y[[i]]$match$flags==2] <- "+"
            showMatchFlags[y[[i]]$match$flags==0] <- "*"
            showMatchFlags[y[[i]]$match$flags==1] <- ""
            showMatch <- matrix(paste(showMatchFlags, y[[i]]$match$multilocus, sep=""), nrow(showMatchFlags), ncol(showMatchFlags))
            dimnames(showMatch) <- dimnames(showMatchFlags)
            
            print(data.frame(rbind(data.frame(showFocal, score=""), data.frame(showMatch, score=y[[i]]$match$score)),
                                   row.names=paste(c("FOCAL  ", rep("MATCH  ", length(y[[i]]$match$index))),
                                                   format(c(y[[i]]$focal$index, y[[i]]$match$index)), "  ",
                                                   format(c(y[[i]]$focal$metaData, y[[i]]$match$metaData)), sep="")))
            cat(y[[i]]$match$perfect, " perfect matches found.  ", y[[i]]$match$partial, " partial matches found.\n", sep="")
            cat("\n\n")
        }
    }   
}


##
## summary.amUnique()
summary.amUnique <- function(object, html=NULL, csv=NULL, ...) {
    
    if (class(object) != "amUnique") {
        stop("allelematch:  this function requires an \"amUnique\" object", call.=FALSE)
    }
    
    if (!is.null(html)) {
        if (is.logical(html) && (html==TRUE)) amHTML.amUnique(object)
        else if (is.logical(html) && (html==FALSE)) html <- NULL
        else amHTML.amUnique(object, htmlFile=html)
    }
    if (!is.null(csv)) {
        amCSV.amUnique(object, csvFile=csv, ...)
    }
    if (is.null(html) && is.null(csv)) {
        cat("allelematch:  Console summary is not available for \"amUnique\" objects.  Please use summary(x, html=TRUE) or summary(x, csv=\"file.csv\") options.\n")
    }
}   


##
## amCluster()
amCluster <- function(amDatasetFocal, runUntilSingletons=TRUE, cutHeight=0.3, missingMethod=2, consensusMethod=1, clusterMethod = "complete") {
  
  ## Check function call variables for validity
  if (!(class(amDatasetFocal) %in% c("amDataset", "amInterpolate", "amCluster"))) {
        stop("allelematch:  amDatasetFocal must be an object of class \"amDataset\" or for recursive use, class \"amCluster\"", call.=FALSE)
  }
  if (class(amDatasetFocal)=="amCluster") {
    amDatasetFocal <- amDatasetFocal$unique
  }
  
  if (class(amDatasetFocal)=="amInterpolate") {
    reClass <- amDatasetFocal
    amDatasetFocal <- list()
    amDatasetFocal$index <- reClass$index
    amDatasetFocal$metaData <- reClass$metaData
    amDatasetFocal$multilocus <- reClass$multilocus
    amDatasetFocal$missingCode <- reClass$missingCode
    class(amDatasetFocal) <- "amDataset"
  }
  
  originalFocalDatasetN <- nrow(amDatasetFocal$multilocus)
  
  if (!(missingMethod %in% c(1,2))) stop("allelematch:  missingMethod must equal 1 or 2", call.=FALSE)
  if (!(consensusMethod %in% c(1,2,3,4))) stop("allelematch:  consensusMethod must equal 1, 2, 3 or 4", call.=FALSE)
  if (!(tolower(clusterMethod) %in% c("complete", "average", "single"))) stop("allelematch:  clusterMethod must be \"complete\" or \"average\"", call.=FALSE)

  totalRuns <- 0
  
  repeat {
    totalRuns <- totalRuns + 1
    
    ## On recursive runs with extreme (and useless) datasets sometimes there will only be one unique individual
    ## In this case the function will return the previous run, which will include clusters
    if (is.null(dim(amDatasetFocal$multilocus))) break
    
    ## Produce dissimilarity matrix
    dissimMatrix <- amMatrix(amDatasetFocal, missingMethod=missingMethod)
    
    ## Do agglomerative hierarchical clustering 
    tryCatch(dendro <- hclust(as.dist(dissimMatrix), method=clusterMethod),
             error=function(x) stop("allelematch:  error in clustering, try a different clusterMethod", call.=FALSE))
    
    ## Do dynamic tree cutting, hybrid method
    tryCatch(uniqueIndex <- cutreeHybrid(dendro, dissimMatrix, cutHeight=cutHeight, minClusterSize=1, verbose=0)$labels+1,
             error=function(x) stop("allelematch:  error in dynamic tree cutting; several causes for this error; have you installed dynamicTreeCut package?  install.packages(\"dynamicTreeCut\")", call.=FALSE))
    numLabels <- length(unique(uniqueIndex))
    
    
    
    ## Build clusterAnalysis object where each $cluster$match are the clusters of genotypes
    ## and $cluster$focal are the consensus genotypes for each cluster
    clusterAnalysis <- list()
    
    ## Create data object to contain clusters (it will be a list of length(0) if there are non)
    clusterAnalysis$cluster <- vector("list", numLabels-1)
  
  
    
    j <- 0
    for (i in 1:numLabels) {
  
      thisGenotype <- matrix(amDatasetFocal$multilocus[uniqueIndex==i, ], sum(uniqueIndex==i), ncol(amDatasetFocal$multilocus), byrow=FALSE)
      dimnames(thisGenotype) <- list(NULL, dimnames(amDatasetFocal$multilocus)[[2]])
      thisIndex <- amDatasetFocal$index[uniqueIndex==i]
      thisMetaData <- amDatasetFocal$metaData[uniqueIndex==i]
      
      ## Produce a list of singletons (label=1) and determine their maximum similarity to any other
      ## in the focalGenotype set.  These values should all be below 50% if the cutHeight was
      ## chosen well.  Then create singletons paired with this next closest match.
      if ((i==1) && (length(thisGenotype) > 1)) {
       
        modifiedDissimMatrix <- dissimMatrix
        diag(modifiedDissimMatrix) <- 1
  
        maxScoreSingletons <- apply(matrix(1 - modifiedDissimMatrix[uniqueIndex==1,], sum(uniqueIndex==1), ncol(modifiedDissimMatrix), byrow=FALSE), 1, max)
  
        maxScoreSingletonsWhich <- apply(matrix(1 - modifiedDissimMatrix[uniqueIndex==1,], sum(uniqueIndex==1), ncol(modifiedDissimMatrix), byrow=FALSE), 1, which.max)
        clusterAnalysis$singletons <- vector("list", sum(uniqueIndex==1))
        
        if (sum(uniqueIndex==1) > 1) {
          reorderSingletons <- order(maxScoreSingletons, decreasing=TRUE)
          maxScoreSingletons <- maxScoreSingletons[reorderSingletons]
          maxScoreSingletonsWhich <- maxScoreSingletonsWhich[reorderSingletons]
          thisGenotype <- thisGenotype[reorderSingletons,]
          thisIndex <- thisIndex[reorderSingletons]
          thisMetaData <- thisMetaData[reorderSingletons]
        }
  
        for (k in 1:sum(uniqueIndex==1)) {
          ## Singleton, focal genotype
          clusterAnalysis$singletons[[k]]$focal$index <- thisIndex[k]
          clusterAnalysis$singletons[[k]]$focal$metaData <- thisMetaData[k]
          clusterAnalysis$singletons[[k]]$focal$multilocus <- t(thisGenotype[k,])
          ## Set all flags to matching
          clusterAnalysis$singletons[[k]]$focal$flags <- matrix(1, nrow(clusterAnalysis$singletons[[k]]$focal$multilocus), ncol(clusterAnalysis$singletons[[k]]$focal$multilocus))
          ## Set flags that are missing
          clusterAnalysis$singletons[[k]]$focal$flags[t(apply(clusterAnalysis$singletons[[k]]$focal$multilocus,1, function(x) x==amDatasetFocal$missingCode))] <- 2
          
          ## Singleton, matching genotype
          clusterAnalysis$singletons[[k]]$match$index <- amDatasetFocal$index[maxScoreSingletonsWhich[k]]
          clusterAnalysis$singletons[[k]]$match$metaData <- amDatasetFocal$metaData[maxScoreSingletonsWhich[k]]
          #browser()
          clusterAnalysis$singletons[[k]]$match$multilocus <- t(amDatasetFocal$multilocus[maxScoreSingletonsWhich[k],])
          clusterAnalysis$singletons[[k]]$match$score <- as.character(signif(maxScoreSingletons[k],2))
          ## Set all flags to matching
          clusterAnalysis$singletons[[k]]$match$flags <- matrix(1, nrow(clusterAnalysis$singletons[[k]]$match$multilocus), ncol(clusterAnalysis$singletons[[k]]$match$multilocus))
          ## Set flags that mismatch
          clusterAnalysis$singletons[[k]]$match$flags[matrix(clusterAnalysis$singletons[[k]]$focal$multilocus,
                                                 nrow(clusterAnalysis$singletons[[k]]$match$multilocus), ncol(clusterAnalysis$singletons[[k]]$match$multilocus), byrow=TRUE)
                                          != clusterAnalysis$singletons[[k]]$match$multilocus] <- 0
          ## Set flags that are missing
          clusterAnalysis$singletons[[k]]$match$flags[t(apply(clusterAnalysis$singletons[[k]]$match$multilocus,1, function(x) x==amDatasetFocal$missingCode))] <- 2
          clusterAnalysis$singletons[[k]]$match$flags[, which(clusterAnalysis$singletons[[k]]$focal$flags==2)] <- 2
          
        }
  
      }
      ## No singletons
      else if (i==1) {
        clusterAnalysis$singletons <- list()
      }
      ## Clusters of individuals
      else {
         j <- j+1
       
        ## Recreate the dissimilarity matrix for simplicity so that only those dissimilarity
        ## scores of for this matching genotype are there.  Assign false index labels so that
        ## amMatrix does not use CPU time creating these unnecessarily.
   
        score <- amMatrix(amDataset(cbind(1:nrow(thisGenotype), thisGenotype), indexColumn=1,
                                      missingCode=amDatasetFocal$missingCode), missingMethod=missingMethod)
         
        ## consensusMethod=1
        ## Find the genotype that has the highest similarity to other genotypes in the cluster and make it the consensus (focal)
        if (consensusMethod==1) {
          lowerTriScore <- score
          lowerTriScore[upper.tri(score, diag=TRUE)] <- 0
          consensusIndex <- which.min(rowSums(lowerTriScore))
          ## Clusters, focal genotype
          clusterAnalysis$cluster[[j]]$focal$index <- thisIndex[consensusIndex]
          clusterAnalysis$cluster[[j]]$focal$metaData <- thisMetaData[consensusIndex]
          clusterAnalysis$cluster[[j]]$focal$multilocus <- t(thisGenotype[consensusIndex, ])
        }
        
        ## consensusMethod=2
        ## Find the genotype that has the highest similarity to other genotypes in the cluster and make it the consensus (focal)
        ## and fill in missing alleles with the mode non-missing allele in each column (not done locus-wise, but column-wise)
        if (consensusMethod==2) {
          lowerTriScore <- score
          lowerTriScore[upper.tri(score, diag=TRUE)] <- 0
          consensusIndex <- which.min(rowSums(lowerTriScore))
          ## Clusters, focal genotype
          clusterAnalysis$cluster[[j]]$focal$index <- thisIndex[consensusIndex]
          clusterAnalysis$cluster[[j]]$focal$metaData <- thisMetaData[consensusIndex]
          clusterAnalysis$cluster[[j]]$focal$multilocus <- t(thisGenotype[consensusIndex, ])
          
          ## Do missing data interpolation for consensus genotype only if required
          if (sum(clusterAnalysis$cluster[[j]]$focal$multilocus==amDatasetFocal$missingCode) > 0) {
  
              ## Determine most frequent (mode) non-missing allele for all columns
              modes <- apply(thisGenotype, 2, function(x) {
                                                        modeAllele <- attr(sort(table(x), decreasing=TRUE), "name")
                                                      
                                                        if (length(modeAllele) > 1) {
                                                          if ((modeAllele[1] == amDatasetFocal$missingCode) && (modeAllele[2] != amDatasetFocal$missingCode)) {
                                                            returnVal <- modeAllele[2]
                                                          }
                                                          else {
                                                            returnVal <- modeAllele[1]
                                                          }
                                                        }
                                                        else {
                                                          returnVal <- modeAllele[1]
                                                        }
                                                        return(returnVal)
                                                      })
              
              ## Reassign missing alleles only with most frequent value
              reassignThese <- clusterAnalysis$cluster[[j]]$focal$multilocus==amDatasetFocal$missingCode
              clusterAnalysis$cluster[[j]]$focal$multilocus[reassignThese] <- modes[reassignThese]
              interpolated <- reassignThese & (modes != amDatasetFocal$missingCode)
  
          }
        }
        
        
        ## consensusMethod=3
        ## Find the genotype that has the least amount of missing data and make it the consensus (focal)
        else if (consensusMethod==3) {
            
          ## Determine counts of missing alleles for all genotypes
          missingAlleles <- rowSums(thisGenotype==amDatasetFocal$missingCode)
          consensusIndex <- which.min(missingAlleles)
          ## Clusters, focal genotype
          clusterAnalysis$cluster[[j]]$focal$index <- thisIndex[consensusIndex]
          clusterAnalysis$cluster[[j]]$focal$metaData <- thisMetaData[consensusIndex]
          clusterAnalysis$cluster[[j]]$focal$multilocus <- t(thisGenotype[consensusIndex, ])
  
          
        }
        
        ## consensusMethod=4
        ## Find the genotype with the least amount of missing data, and fill in missing alleles with
        ## the mode non-missing allele in each column (not done locus-wise, but column-wise)
        else if (consensusMethod==4) {
        
          ## Determine counts of missing alleles for all genotypes
          missingAlleles <- rowSums(thisGenotype==amDatasetFocal$missingCode)
          consensusIndex <- which.min(missingAlleles)
          ## Clusters, focal genotype
          clusterAnalysis$cluster[[j]]$focal$index <- thisIndex[consensusIndex]
          clusterAnalysis$cluster[[j]]$focal$metaData <- thisMetaData[consensusIndex]
          clusterAnalysis$cluster[[j]]$focal$multilocus <- t(thisGenotype[consensusIndex, ])
  
          ## Do missing data interpolation for consensus genotype only if required
          if (sum(clusterAnalysis$cluster[[j]]$focal$multilocus==amDatasetFocal$missingCode) > 0) {
  
              ## Determine most frequent (mode) non-missing allele for all columns
              modes <- apply(thisGenotype, 2, function(x) {
                                                        modeAllele <- attr(sort(table(x), decreasing=TRUE), "name")
                                                      
                                                        if (length(modeAllele) > 1) {
                                                          if ((modeAllele[1] == amDatasetFocal$missingCode) && (modeAllele[2] != amDatasetFocal$missingCode)) {
                                                            returnVal <- modeAllele[2]
                                                          }
                                                          else {
                                                            returnVal <- modeAllele[1]
                                                          }
                                                        }
                                                        else {
                                                          returnVal <- modeAllele[1]
                                                        }
                                                        return(returnVal)
                                                      })
              
              ## Reassign missing alleles only with most frequent value
              reassignThese <- clusterAnalysis$cluster[[j]]$focal$multilocus==amDatasetFocal$missingCode
              clusterAnalysis$cluster[[j]]$focal$multilocus[reassignThese] <- modes[reassignThese]
              interpolated <- reassignThese & (modes != amDatasetFocal$missingCode)
  
          }
        }
    
      
      ## Set all focal flags to matching
      clusterAnalysis$cluster[[j]]$focal$flags <- matrix(1, nrow(clusterAnalysis$cluster[[j]]$focal$multilocus), ncol(clusterAnalysis$cluster[[j]]$focal$multilocus))
      ## Set all focal flags that are missing
      clusterAnalysis$cluster[[j]]$focal$flags[t(apply(clusterAnalysis$cluster[[j]]$focal$multilocus,1, function(x) x==amDatasetFocal$missingCode))] <- 2
      ## Set all focal flags that are interpolated
      if (exists("interpolated")) clusterAnalysis$cluster[[j]]$focal$flags[interpolated] <- 3
   
      ## Clusters, matching genotypes
      reorderMatch <- order(1-score[consensusIndex,], decreasing=TRUE)
      clusterAnalysis$cluster[[j]]$match$index <- thisIndex[reorderMatch]
      clusterAnalysis$cluster[[j]]$match$metaData <- thisMetaData[reorderMatch]
      clusterAnalysis$cluster[[j]]$match$multilocus <- thisGenotype[reorderMatch, ]
      clusterAnalysis$cluster[[j]]$match$score <- as.character(signif(1-score[consensusIndex, ], 2)[reorderMatch])
      ## Set all flags to matching
      clusterAnalysis$cluster[[j]]$match$flags <- matrix(1, nrow(clusterAnalysis$cluster[[j]]$match$multilocus), ncol(clusterAnalysis$cluster[[j]]$match$multilocus))
      ## Set flags that mismatch
      clusterAnalysis$cluster[[j]]$match$flags[matrix(clusterAnalysis$cluster[[j]]$focal$multilocus,
                                             nrow(clusterAnalysis$cluster[[j]]$match$multilocus), ncol(clusterAnalysis$cluster[[j]]$match$multilocus), byrow=TRUE)
                                      != clusterAnalysis$cluster[[j]]$match$multilocus] <- 0
      ## Set flags that are missing
      clusterAnalysis$cluster[[j]]$match$flags[t(apply(clusterAnalysis$cluster[[j]]$match$multilocus,1, function(x) x==amDatasetFocal$missingCode))] <- 2
      clusterAnalysis$cluster[[j]]$match$flags[, which(clusterAnalysis$cluster[[j]]$focal$flags==2)] <- 2
      
      ## Summarize the matches
      clusterAnalysis$cluster[[j]]$match$perfect <- sum(apply(clusterAnalysis$cluster[[j]]$match$flags, 1, function(x) all(x==1)))
      clusterAnalysis$cluster[[j]]$match$partial <- nrow(clusterAnalysis$cluster[[j]]$match$flags) - clusterAnalysis$cluster[[j]]$match$perfect
      }
  
    }
    
    ## Produce a summary of unique genotypes
    ## No clusters
    if (length(clusterAnalysis$cluster) == 0) {
      clusterAnalysis$unique$index <- do.call(c,lapply(clusterAnalysis$singletons, function(x) x$focal$index))
      clusterAnalysis$unique$metaData <- do.call(c,lapply(clusterAnalysis$singletons, function(x) x$focal$metaData))
      clusterAnalysis$unique$multilocus <- do.call(rbind,lapply(clusterAnalysis$singletons, function(x) x$focal$multilocus))
      clusterAnalysis$unique$uniqueType <-  rep("SINGLETON", length(clusterAnalysis$singletons))
    }
    ## No singletons
    else if (length(clusterAnalysis$singletons) == 0) {
      clusterAnalysis$unique$index <- do.call(c,lapply(clusterAnalysis$cluster, function(x) x$focal$index))
      clusterAnalysis$unique$metaData <- do.call(c,lapply(clusterAnalysis$cluster, function(x) x$focal$metaData))
      clusterAnalysis$unique$multilocus <- do.call(rbind,lapply(clusterAnalysis$cluster, function(x) x$focal$multilocus))
      clusterAnalysis$unique$uniqueType <-  rep("CONSENSUS", length(clusterAnalysis$cluster))
    }
    ## Singletons and clusters
    else {
      clusterAnalysis$unique$index <- c(do.call(c,lapply(clusterAnalysis$cluster, function(x) x$focal$index)),
                                          do.call(c,lapply(clusterAnalysis$singletons, function(x) x$focal$index)))
      clusterAnalysis$unique$metaData <- c(do.call(c,lapply(clusterAnalysis$cluster, function(x) x$focal$metaData)),
                                          do.call(c,lapply(clusterAnalysis$singletons, function(x) x$focal$metaData)))
      clusterAnalysis$unique$multilocus <- rbind(do.call(rbind,lapply(clusterAnalysis$cluster, function(x) x$focal$multilocus)),
                                               do.call(rbind,lapply(clusterAnalysis$singletons, function(x) x$focal$multilocus)))
      clusterAnalysis$unique$uniqueType <-  c(rep("CONSENSUS", length(clusterAnalysis$cluster)),
                                                 rep("SINGLETON", length(clusterAnalysis$singletons)))
      
    }

    
   
    ## Sort by the index column
    ## First check to see if first column should be sorted as a character or as a numeric value
    if (sum(is.na(suppressWarnings(as.numeric(clusterAnalysis$unique$index)))) > 0) {
      orderUnique <- order(clusterAnalysis$unique$index)
    }
    else {
      orderUnique <- order(as.numeric(clusterAnalysis$unique$index))
    }
    clusterAnalysis$unique$index <- clusterAnalysis$unique$index[orderUnique]
    clusterAnalysis$unique$metaData <- clusterAnalysis$unique$metaData[orderUnique]
    clusterAnalysis$unique$multilocus <- clusterAnalysis$unique$multilocus[orderUnique, ]
    clusterAnalysis$unique$uniqueType <- clusterAnalysis$unique$uniqueType[orderUnique]
    clusterAnalysis$unique$missingCode <- amDatasetFocal$missingCode
    class(clusterAnalysis$unique) <- "amDataset"
    
    clusterAnalysis$cutHeight <- cutHeight
    clusterAnalysis$consensusMethod <- consensusMethod
    clusterAnalysis$missingMethod <- missingMethod
    clusterAnalysis$clusterMethod <- clusterMethod
    clusterAnalysis$missingCode <- amDatasetFocal$missingCode
    clusterAnalysis$focalDatasetN <- originalFocalDatasetN
    clusterAnalysis$totalRuns <- totalRuns
    clusterAnalysis$runUntilSingletons <- runUntilSingletons
    class(clusterAnalysis) <- "amCluster"
    
    
    if (runUntilSingletons) {
        if (length(clusterAnalysis$cluster)==0) break
        else amDatasetFocal <- clusterAnalysis$unique
    }
    else break;

  }
  ## Check that the multilocus output is satisfactory
  if (is.null(clusterAnalysis$unique$multilocus)) {
    stop("allelematch:  amCluster:  no clusters formed.  Please set cutHeight lower and run again.", call.=FALSE)
    } 
  if (is.null(dim(clusterAnalysis$unique$multilocus))) {
        tmpMultilocus <- matrix(clusterAnalysis$unique$multilocus, 1, length(clusterAnalysis$unique$multilocus),byrow=FALSE)
        dimnames(tmpMultilocus)[[2]] <- names(clusterAnalysis$unique$multilocus)
        clusterAnalysis$unique$multilocus <- tmpMultilocus
      }
  return(clusterAnalysis)
}
##
## summary.amCluster()
summary.amCluster <- function(object, html=NULL, csv=NULL, ...) {
    
    if (class(object) != "amCluster") {
        stop("allelematch:  this function requires an \"amCluster\" object", call.=FALSE)
    }
    
    if (!is.null(html)) {
        if (is.logical(html) && (html==TRUE)) amHTML.amCluster(object)
        else if (is.logical(html) && (html==FALSE)) html <- NULL
        else amHTML.amCluster(object, htmlFile=html)
        }
    if (!is.null(csv)) {
        amCSV.amCluster(object, csvFile=csv)
    }
    if (is.null(html) && is.null(csv)) {
        cat("allelematch\namCluster object\n")
        cat("\nFocal dataset N=", object$focalDatasetN, "\n", sep="")
        cat("Unique genotypes (total): ", nrow(object$unique$multilocus), "\n", sep="")
        cat("Unique genotypes (by cluster consensus): ", length(object$cluster), "\n", sep="")
        cat("Unique genotypes (singletons): ", length(object$singletons), "\n\n", sep="")
        cat("Missing data represented by: ", object$missingCode, "\n", sep="")
        cat("Missing data matching method: ", object$missingMethod, "\n", sep="")
        cat("Clustered genotypes consensus method: ", object$consensusMethod, "\n", sep="")
        cat("Hierarchical clustering method: ", object$clusterMethod, "\n", sep="")
        cat("Dynamic tree cutting height (cutHeight): ", object$cutHeight, "\n\n", sep="")
        cat("Run until only singletons: ", object$runUntilSingletons, "\n", sep="")
        cat("Runs: ", object$totalRuns, "\n\n", sep="")
        cat("Score flags:\n")
        cat("*101 Allele does not match\n")
        cat("+101 Allele is missing\n")
        cat("[101] Allele is interpolated in consensus from non-missing cluster members (consensusMethod = 2, 4 only)\n\n")
        
        ## Write clusters
        y <- object$cluster
        if (length(y) > 0) {
            for (i in 1:length(y)) {
                cat("(Cluster ", i, " of ", length(y), ")\n", sep="")
                if (is.null(y[[i]]$focal$metaData)) y[[i]]$focal$metaData <- ""
                if (is.null(y[[i]]$match$metaData)) y[[i]]$match$metaData <- rep("", length(y[[i]]$match$index))
                
                showFocalFlags <- y[[i]]$focal$multilocus
                showFocalFlags[y[[i]]$focal$flags==2] <- "+"
                showFocalFlags[y[[i]]$focal$flags==0] <- "*"
                showFocalFlags[y[[i]]$focal$flags==1] <- ""
                showFocalFlags[y[[i]]$focal$flags==3] <- "["
                showFocal <- matrix(paste(showFocalFlags, y[[i]]$focal$multilocus, sep=""), nrow(showFocalFlags), ncol(showFocalFlags))
                showFocal[, grepl("\\[", showFocal)] <- paste(showFocal[,grepl("\\[", showFocal)], "]", sep="")
                dimnames(showFocal) <- dimnames(showFocalFlags)
                
                showMatchFlags <- y[[i]]$match$multilocus
                showMatchFlags[y[[i]]$match$flags==2] <- "+"
                showMatchFlags[y[[i]]$match$flags==0] <- "*"
                showMatchFlags[y[[i]]$match$flags==1] <- ""
                showMatch <- matrix(paste(showMatchFlags, y[[i]]$match$multilocus, sep=""), nrow(showMatchFlags), ncol(showMatchFlags))
                dimnames(showMatch) <- dimnames(showMatchFlags)
                
                print(data.frame(rbind(data.frame(showFocal, score=""), data.frame(showMatch, score=y[[i]]$match$score)),
                                       row.names=paste(c("CONSENSUS  ", rep("MEMBER     ", length(y[[i]]$match$index))),
                                                       format(c(y[[i]]$focal$index, y[[i]]$match$index)), "  ",
                                                       format(c(y[[i]]$focal$metaData, y[[i]]$match$metaData)), sep="")))
                cat(y[[i]]$match$perfect, " perfect matches found.  ", y[[i]]$match$partial, " partial matches found.\n", sep="")
                cat("\n\n")
            }
            cat("(Clusters END)\n\n\n")
        }
    
        ## Write singletons
        y <- object$singletons
        if (length(y) > 0) {
            for (i in 1:length(y)) {
                cat("(Singleton ", i, " of ", length(y), ")\n", sep="")
                if (is.null(y[[i]]$focal$metaData)) y[[i]]$focal$metaData <- ""
                if (is.null(y[[i]]$match$metaData)) y[[i]]$match$metaData <- rep("", length(y[[i]]$match$index))
                
                showFocalFlags <- y[[i]]$focal$multilocus
                showFocalFlags[y[[i]]$focal$flags==2] <- "+"
                showFocalFlags[y[[i]]$focal$flags==0] <- "*"
                showFocalFlags[y[[i]]$focal$flags==1] <- ""
                showFocal <- matrix(paste(showFocalFlags, y[[i]]$focal$multilocus, sep=""), nrow(showFocalFlags), ncol(showFocalFlags))
                dimnames(showFocal) <- dimnames(showFocalFlags)
                
                showMatchFlags <- y[[i]]$match$multilocus
                showMatchFlags[y[[i]]$match$flags==2] <- "+"
                showMatchFlags[y[[i]]$match$flags==0] <- "*"
                showMatchFlags[y[[i]]$match$flags==1] <- ""
                showMatch <- matrix(paste(showMatchFlags, y[[i]]$match$multilocus, sep=""), nrow(showMatchFlags), ncol(showMatchFlags))
                dimnames(showMatch) <- dimnames(showMatchFlags)
                
                print(data.frame(rbind(data.frame(showFocal, score=""), data.frame(showMatch, score=y[[i]]$match$score)),
                                       row.names=paste(c("SINGLETON  ", rep("CLOSEST    ", length(y[[i]]$match$index))),
                                                       format(c(y[[i]]$focal$index, y[[i]]$match$index)), "  ",
                                                       format(c(y[[i]]$focal$metaData, y[[i]]$match$metaData)), sep="")))
                cat("\n\n")
            }
            cat("(Singletons END)\n\n\n")
        }
            
        ## Write unique
        y <- object$unique
        cat("(Unique genotypes)\n")
        if (is.null(y$metaData)) y$metaData <- ""     
        print(data.frame(y$multilocus, row.names=paste(y$uniqueType, "  ", format(y$index), "  ", format(y$metaData), sep="")))
    }    
        
}
    
    
    




#
## amAlleleFreq()
amAlleleFreq <- function(amDatasetFocal, multilocusMap=NULL) {
    
    if (class(amDatasetFocal) != "amDataset") {
        stop("allelematch:  amDatasetFocal must be an object of class \"amDataset\"", call.=FALSE)
    }
    
    ## Set multilocusMap to default if not given
    if (is.null(multilocusMap)) {
        if ((ncol(amDatasetFocal$multilocus)%%2) != 0) {
            stop("allelematch:  there are an odd number of genotype columns in amDatasetFocal; Please specify multilocusMap manually", call.=FALSE)
        }
        else cat("allelematch:  assuming genotype columns are in pairs, representing", ncol(amDatasetFocal$multilocus)/2, "loci\n")
        multilocusMap <- rep(1:(ncol(amDatasetFocal$multilocus)/2), each=2)
    }
    ## Check multilocusMap is the correct length
    else if (length(multilocusMap) != ncol(amDatasetFocal$multilocus))  {
        stop("allelematch:  multilocusMap must be a vector of integers or strings giving the mappings onto loci for all genotype columns in amDatasetFocal;
             Example: gender followed by 4 diploid loci in paired columns could be coded: mutlilocusMap=c(1,2,2,3,3,4,4,5,5)
             or as: multilocusMap=c(\"GENDER\",\"LOC1\",\"LOC1\",\"LOC2\",\"LOC2\",\"LOC3\",\"LOC3\",\"LOC4\",\"LOC4\")", call.=FALSE)
    }
    else if (sum(table(multilocusMap) > 2) > 0) {
        stop("allelematch:  multilocusMap indicates that a locus occurs in three or more columns;  this situation is not yet handled", call.=FALSE)
    }
    multilocusMap <- as.integer(as.factor(multilocusMap))
    
    alleleFreq <- list()
    alleleFreq$multilocusMap <- multilocusMap
    alleleFreq$loci <- vector("list", max(multilocusMap))
    ## Determine allele frequencies for each locus
    for (locus in 1:max(multilocusMap)) {

        thisLocusIndex <- which(multilocusMap==locus)
        alleleFreq$loci[[locus]]$name <- paste(dimnames(amDatasetFocal$multilocus)[[2]][thisLocusIndex], collapse="-")
        alleleFreq$loci[[locus]]$columnNames <- dimnames(amDatasetFocal$multilocus)[[2]][thisLocusIndex]
        thisLocus <- as.character(amDatasetFocal$multilocus[, thisLocusIndex])
        thisLocus[thisLocus==amDatasetFocal$missingCode] <- NA
        thisLocusUnique <- unique(thisLocus)[!is.na(unique(thisLocus))]
        alleleFreq$loci[[locus]]$alleleFreq <- sort(sapply(thisLocusUnique, function(x) sum(thisLocus==x, na.rm=TRUE)/length(thisLocus[!is.na(thisLocus)])), decreasing=TRUE)
        alleleFreq$loci[[locus]]$missingFreq <- sum(is.na(thisLocus))/length(thisLocus)
        alleleFreq$loci[[locus]]$numAlleles <- length(thisLocusUnique)
        #alleleFreq$loci[[locus]]$PIC <- 1 - sum(alleleFreq$loci[[locus]]$alleleFreq^2) - sum(apply(t(combn(1:length(alleleFreq$loci[[locus]]$alleleFreq), 2)), 1, function(x) prod(alleleFreq$loci[[locus]]$alleleFreq[x]^2)))
    }
    class(alleleFreq) <- "amAlleleFreq"

    return(alleleFreq)
}
#
## print.amAlleleFreq()
print.amAlleleFreq <- function(x, ...) {
    
    cat("allelematch\namAlleleFreq object\nFrequencies calculated after removal of missing data\n")
    y <- x$loci
    for (i in 1:length(y)) {
        cat("\n", locus=y[[i]]$name, " (", y[[i]]$numAlleles, " alleles)\n", sep="")
        for (j in 1:length(y[[i]]$alleleFreq)) {
            cat("\tAllele\t", names(y[[i]]$alleleFreq)[j], "\t", signif(y[[i]]$alleleFreq[j],3), "\n")
        }
    }
}
##
## amUnique()
amUnique <- function(amDatasetFocal, multilocusMap=NULL, alleleMismatch=NULL, matchThreshold=NULL, cutHeight=NULL, doPsib="missing", consensusMethod=1, verbose=FALSE) {
    
    if (!class(amDatasetFocal)=="amDataset") {
        stop("allelematch:  amDatasetFocal must be an object of class \"amDataset\"", call.=FALSE)
    }
    
    # Set multilocusMap to default if not given
    if (is.null(multilocusMap)) {
        if ((ncol(amDatasetFocal$multilocus)%%2) != 0) {
            stop("allelematch:  there are an odd number of genotype columns in amDatasetFocal; Please specify multilocusMap manually", call.=FALSE)
        }
        else cat("allelematch:  assuming genotype columns are in pairs, representing", ncol(amDatasetFocal$multilocus)/2, "loci\n")
        multilocusMap <- rep(1:(ncol(amDatasetFocal$multilocus)/2), each=2)
    }
    ## Check multilocusMap is the correct length
    else if (length(multilocusMap) != ncol(amDatasetFocal$multilocus))  {
        stop("allelematch:  multilocusMap must be a vector of integers or strings giving the mappings onto loci for all genotype columns in amDatasetFocal;
             Example: gender followed by 4 diploid loci in paired columns could be coded: mutlilocusMap=c(1,2,2,3,3,4,4,5,5)
             or as: multilocusMap=c(\"GENDER\",\"LOC1\",\"LOC1\",\"LOC2\",\"LOC2\",\"LOC3\",\"LOC3\",\"LOC4\",\"LOC4\")", call.=FALSE)
    }
    else if (sum(table(multilocusMap) > 2) > 0) {
        stop("allelematch:  multilocusMap indicates that a locus occurs in three or more columns;  this situation is not yet handled", call.=FALSE)
    }
    multilocusMap <- as.integer(as.factor(multilocusMap))
    
    ## More checking of input parameters
    if (sum(!(c(is.null(alleleMismatch), is.null(matchThreshold), is.null(cutHeight)))) != 1) {
        stop("allelematch:  please specify alleleMismatch OR matchThreshold OR cutHeight.", call.=FALSE)
    }
    
    if (length(c(alleleMismatch, matchThreshold, cutHeight))>1) {
            stop("allelematch:  please provide a single parameter value for alleleMismatch OR matchThreshold OR cutHeight.  Use amUniqueProfile() to examine a range of values", call.=FALSE)
        }
    
    if (!is.null(matchThreshold)) {
        if ((matchThreshold < 0) || (matchThreshold > 1)) {
                stop("allelematch:  matchThreshold must be between 0 and 1", call.=FALSE)
            }
        cutHeight <- 1-matchThreshold
        alleleMismatch <- round((1-matchThreshold)*length(multilocusMap),2)
    }
    else if (!is.null(alleleMismatch)) {
            matchThreshold <- 1-(alleleMismatch/length(multilocusMap))
            cutHeight <- 1-matchThreshold
        }
    else if (!is.null(cutHeight)) {
        if ((cutHeight < 0) || (cutHeight > 1)) {
                stop("allelematch:  cutHeight must be greater than 0 and less than 1", call.=FALSE)
            }
            matchThreshold <- 1-cutHeight
            alleleMismatch <- round((1-matchThreshold)*length(multilocusMap),2)
        }

    if (matchThreshold==1 && cutHeight==0) {
        if (verbose) cat("allelematch: cutHeight cannot be zero.  Setting cutHeight=0.00001.  This will return perfect matches.\n")
        cutHeight <- 0.00001
    }
    
    ## Run the required analyses
    if (verbose) cat("allelematch:  amUnique:  Clustering genotypes\n")
    clusterAnalysis <- amCluster(amDatasetFocal, cutHeight=cutHeight, runUntilSingletons=TRUE, consensusMethod=consensusMethod)
    
    if (verbose) cat("allelematch:  amUnique:  Comparing unique genotypes identified by clustering to all samples\n")
    clusterAnalysisPairwise <- amPairwise(clusterAnalysis$unique, amDatasetFocal, matchThreshold=matchThreshold)

    if (verbose) cat("allelematch:  amUnique:  Determining allele frequencies of unique genotypes identified by cluster\n")
    clusterAnalysisAlleleFreq <- amAlleleFreq(clusterAnalysis$unique, multilocusMap=multilocusMap)
    
    if (verbose) cat("allelematch:  amUnique:  Finding Psib\n")
    ## Prepare an output object of class amUnique
    uniqueAnalysis <- clusterAnalysisPairwise
    class(uniqueAnalysis) <- "amUnique"
    for (i in 1:length(uniqueAnalysis$pairwise)) {
        uniqueAnalysis$pairwise[[i]]$match$psib <- rep(NA, nrow(uniqueAnalysis$pairwise[[i]]$match$multilocus))
        uniqueAnalysis$pairwise[[i]]$match$rowFlag <- rep("", nrow(uniqueAnalysis$pairwise[[i]]$match$multilocus))
    }

    ## Find Psib 
    multilocusMap <- clusterAnalysisAlleleFreq$multilocusMap
    ## Loop over all pairwise comparison sets
    for (iPairwise in 1:length(uniqueAnalysis$pairwise)) {
        thisPairwise <- uniqueAnalysis$pairwise[[iPairwise]]
        ## Find all the genotype rows that have no mismatches if doPsib=="missing" 
        if (doPsib != "all") {
            doTheseGenotypes <- rowSums(thisPairwise$match$flags==0)==0
        }
        ## Or if doPsib=="all" return all genotype rows
        else {
            doTheseGenotypes <- rep(TRUE, nrow(thisPairwise$match$flags))
        }

        if (sum(doTheseGenotypes) > 0) {
            for (iGenotype in which(doTheseGenotypes)) {
                thisGenotype <- thisPairwise$match$multilocus[iGenotype, ]
                psib <- 1
                ## Loop over all loci for this multilocus genotype
                for (iLocus in 1:max(multilocusMap)) {
                    thisLocus <- thisGenotype[which(multilocusMap==iLocus)]
                    ## Continue if this locus is not missing or mismatched in match or in focal (i.e. flags are 1)
                    if (all(thisPairwise$match$flags[iGenotype, which(multilocusMap==iLocus)]==1)) {
                        ## Single column locus (e.g. gender)
                        if (length(thisLocus)==1) {
                            alleleA <- thisLocus   
                            pA <- clusterAnalysisAlleleFreq$loci[[iLocus]]$alleleFreq[alleleA]
                            psib <- psib * (1 + (2 * pA) + (pA*pA)) /4    
                        }
                        ## Two column locus
                        else {
                            alleleA <- thisLocus[1]
                            alleleB <- thisLocus[2]
                            ## Homozygote case
                            if (alleleA==alleleB) {
                                pA <- clusterAnalysisAlleleFreq$loci[[iLocus]]$alleleFreq[alleleA]
                                psib <- psib * (1 + (2 * pA) + (pA*pA)) /4
                            }   
                            ## Heterozygote case
                            else {
                                pA <- clusterAnalysisAlleleFreq$loci[[iLocus]]$alleleFreq[alleleA]
                                pB <- clusterAnalysisAlleleFreq$loci[[iLocus]]$alleleFreq[alleleB]
                                psib <- psib * (1 + pA + pB + (2*pA*pB)) / 4
                            }
                        }
                    }
                   # cat("iPairwise=", iPairwise, " iGenotype=", iGenotype, " iLocus=", iLocus, " psib=", psib, "\n", sep="")
                }
                ## Save the psib we calculated only if focal and match sample index are not the same
                #if (uniqueAnalysis$pairwise[[iPairwise]]$focal$index != uniqueAnalysis$pairwise[[iPairwise]]$match$index[iGenotype]) {
                #    uniqueAnalysis$pairwise[[iPairwise]]$match$psib[iGenotype] <- psib
                #}

                ## Save the psib we calculated
                uniqueAnalysis$pairwise[[iPairwise]]$match$psib[iGenotype] <- psib

                uniqueAnalysis$pairwise[[iPairwise]]$match$psibNotCalculable <- sum(!doTheseGenotypes)
            }
        }
    }
    
    ## Identify samples that do not appear either as focal or in pairwise matches at this matchThreshold
    indexUnclassified <- amDatasetFocal$index[!(amDatasetFocal$index %in% unique(do.call(c,lapply(clusterAnalysisPairwise$pairwise, function(x) c(x$focal$index, x$match$index)))))]
    uniqueAnalysis$numUnclassified <- length(indexUnclassified)
    ## If there are any samples not matched, add a section to the output object of
    ## class amDataset containing only these not matched or "missing" rows
    if (uniqueAnalysis$numUnclassified > 0) {
        uniqueAnalysis$unclassified <- amDatasetFocal
        unclassifiedDatasetFocal <- which(amDatasetFocal$index %in% indexUnclassified)
        uniqueAnalysis$unclassified$index <- uniqueAnalysis$unclassified$index[unclassifiedDatasetFocal]
        if (!is.null(uniqueAnalysis$unclassified$metaData)) {
            uniqueAnalysis$unclassified$metaData <- uniqueAnalysis$unclassified$metaData[unclassifiedDatasetFocal]
        }
        tmpMultilocus <- matrix(uniqueAnalysis$unclassified$multilocus[unclassifiedDatasetFocal, ], length(indexUnclassified), ncol(uniqueAnalysis$unclassified$multilocus),byrow=FALSE)
        dimnames(tmpMultilocus)[[2]] <- dimnames(uniqueAnalysis$unclassified$multilocus)[[2]]
        uniqueAnalysis$unclassified$multilocus <- tmpMultilocus
    }
    
    ## Identify samples that appear in more than one pairwise match
    ## by placing a flag indicating this in the rowFlag column
    ## and create an amDataset object to contain these for external use
    indexMatches <- do.call(c,lapply(clusterAnalysisPairwise$pairwise, function(x) x$match$index))
    indexMultipleMatches <- unique(indexMatches[duplicated(indexMatches)])
    uniqueAnalysis$numMultipleMatches <- length(indexMultipleMatches)
    for (i in 1:length(uniqueAnalysis$pairwise)) {
        theseMultipleMatches <- which(uniqueAnalysis$pairwise[[i]]$match$index %in% indexMultipleMatches)
        if (length(theseMultipleMatches != 0)) {
            uniqueAnalysis$pairwise[[i]]$match$rowFlag[theseMultipleMatches] <- "MULTIPLE_MATCH"
            uniqueAnalysis$pairwise[[i]]$focal$rowFlag <- "CHECK_UNIQUE"
        }
        else {
            uniqueAnalysis$pairwise[[i]]$focal$rowFlag <- "UNIQUE"
        }
    }
    if (length(indexMultipleMatches > 0)) {
        uniqueAnalysis$multipleMatches <- amDatasetFocal
        multipleMatchesDatasetFocal <- which(amDatasetFocal$index %in% indexMultipleMatches)
        uniqueAnalysis$multipleMatches$index <- uniqueAnalysis$multipleMatches$index[multipleMatchesDatasetFocal]
        if (!is.null(uniqueAnalysis$multipleMatches$metaData)) {
            uniqueAnalysis$multipleMatches$metaData <- uniqueAnalysis$multipleMatches$metaData[multipleMatchesDatasetFocal]
        }
        tmpMultilocus <- matrix(uniqueAnalysis$multipleMatches$multilocus[multipleMatchesDatasetFocal, ], length(indexMultipleMatches), ncol(uniqueAnalysis$multipleMatches$multilocus),byrow=FALSE)
        dimnames(tmpMultilocus)[[2]] <- dimnames(uniqueAnalysis$multipleMatches$multilocus)[[2]]
        uniqueAnalysis$multipleMatches$multilocus <- tmpMultilocus
    }

    uniqueAnalysis$unique <- clusterAnalysis$unique
    uniqueAnalysis$unique$psib <- unlist(lapply(uniqueAnalysis$pairwise, function(x) x$match$psib[1]))
    uniqueAnalysis$unique$rowFlag <- unlist(lapply(uniqueAnalysis$pairwise, function(x) x$focal$rowFlag))
    
    
    ## Add in other reference items involved in the analysis
    uniqueAnalysis$cutHeight <- cutHeight
    uniqueAnalysis$consensusMethod <- clusterAnalysis$consensusMethod
    uniqueAnalysis$alleleMismatch <- alleleMismatch
    uniqueAnalysis$doPsib <- doPsib
    uniqueAnalysis$alleleFreq <- clusterAnalysisAlleleFreq
    
    
    return(uniqueAnalysis)

    
    
}
##
## amUniqueProfile()
amUniqueProfile <- function(amDatasetFocal, multilocusMap=NULL, alleleMismatch=NULL, matchThreshold=NULL, cutHeight=NULL, guessOptimum=TRUE, doPlot=TRUE, consensusMethod=1, verbose=TRUE) {
    
    if (!class(amDatasetFocal)=="amDataset") {
        stop("allelematch:  amDatasetFocal must be an object of class \"amDataset\"", call.=FALSE)
    }
    
    # Set multilocusMap to default if not given
    if (is.null(multilocusMap)) {
        if ((ncol(amDatasetFocal$multilocus)%%2) != 0) {
            stop("allelematch:  there are an odd number of genotype columns in amDatasetFocal; Please specify multilocusMap manually", call.=FALSE)
        }
        else cat("allelematch:  assuming genotype columns are in pairs, representing", ncol(amDatasetFocal$multilocus)/2, "loci\n")
        multilocusMap <- rep(1:(ncol(amDatasetFocal$multilocus)/2), each=2)
    }
    ## Check multilocusMap is the correct length
    else if (length(multilocusMap) != ncol(amDatasetFocal$multilocus))  {
        stop("allelematch:  multilocusMap must be a vector of integers or strings giving the mappings onto loci for all genotype columns in amDatasetFocal;
             Example: gender followed by 4 diploid loci in paired columns could be coded: mutlilocusMap=c(1,2,2,3,3,4,4,5,5)
             or as: multilocusMap=c(\"GENDER\",\"LOC1\",\"LOC1\",\"LOC2\",\"LOC2\",\"LOC3\",\"LOC3\",\"LOC4\",\"LOC4\")", call.=FALSE)
    }
    else if (sum(table(multilocusMap) > 2) > 0) {
        stop("allelematch:  multilocusMap indicates that a locus occurs in three or more columns;  this situation is not yet handled", call.=FALSE)
    }
    multilocusMap <- as.integer(as.factor(multilocusMap))

    ## More checking of input parameters
    if (sum(!(c(is.null(alleleMismatch), is.null(matchThreshold), is.null(cutHeight)))) > 1) {
        stop("allelematch:  please specify alleleMismatch OR matchThreshold OR cutHeight.", call.=FALSE)
    }
    
    if (sum(!(c(is.null(alleleMismatch), is.null(matchThreshold), is.null(cutHeight))))==0) {
        alleleMismatch <- seq(0, floor(length(multilocusMap))*0.4, 1)
        matchThreshold <- 1-(alleleMismatch/length(multilocusMap))
        cutHeight <- 1-matchThreshold
        profileType <- "alleleMismatch"
    }
    else {
        if (length(c(alleleMismatch, matchThreshold, cutHeight))<2) {
            stop("allelematch:  please provide a range of parameter values for alleleMismatch OR matchThreshold OR cutHeight.  e.g. alleleMismatch=c(0,1,2,3,4,5,6,7,8)", call.=FALSE)
        }
    
        if (!is.null(alleleMismatch)) {
            if (length(alleleMismatch) <= 2) {
                stop("allelematch:  alleleMismatch must be a vector containing a sequence of three or more values", call.=FALSE)
            }
            matchThreshold <- 1-(alleleMismatch/length(multilocusMap))
            cutHeight <- 1-matchThreshold
            profileType <- "alleleMismatch"
        }
        else if (!is.null(cutHeight)) {
            if ((any(cutHeight < 0)) || (any(cutHeight > 1))) {
                stop("allelematch:  cutHeight must be between 0 and 1", call.=FALSE)
            }
            if(length(cutHeight <= 2)) {
                stop("allelematch:  cutHeight must be a vector containing a sequence of three or more values", call.=FALSE)    
            }
            matchThreshold <- 1-cutHeight
            alleleMismatch <- round((1-matchThreshold)*length(multilocusMap),2)
            profileType <- "cutHeight"
        }
        else {
            if ((any(matchThreshold < 0)) || (any(matchThreshold > 1))) {
                stop("allelematch:  matchThreshold must be between 0 and 1", call.=FALSE)
            }
            if(length(matchThreshold <= 2)) {
                stop("allelematch:  matchThreshold must be a vector containing a sequence of three or more values", call.=FALSE)    
            }
            cutHeight <- 1-matchThreshold
            alleleMismatch <- round((1-matchThreshold)*length(multilocusMap),2)
            profileType <- "matchThreshold"
        }
    }
    if (verbose) cat("allelematch:  running amUnique() at", length(matchThreshold), "different values of", profileType, "\n")
    
    profileResults <- data.frame(matchThreshold=NA, cutHeight=NA, alleleMismatch=NA, samples=NA, unique=NA, unclassified=NA, multipleMatch=NA, guessOptimum=NA)
    
    for (i in 1:length(matchThreshold)) {
        if (verbose) cat("allelematch:  ", i, " of ", length(matchThreshold), " (matchThreshold=", round(matchThreshold[i],2), ", cutHeight=", round(cutHeight[i],2), ", alleleMismatch=", alleleMismatch[i], ")\n", sep="")
        profileResults[i, "matchThreshold"] <- matchThreshold[i]
        profileResults[i, "cutHeight"] <- cutHeight[i]
        profileResults[i, "alleleMismatch"] <- alleleMismatch[i]
        amUniqueResult <- amUnique(amDatasetFocal, matchThreshold=matchThreshold[i], multilocusMap=multilocusMap, verbose=FALSE)
        profileResults[i, "unclassified"] <- amUniqueResult$numUnclassified
        profileResults[i, "multipleMatch"] <- amUniqueResult$numMultipleMatch
        profileResults[i, "unique"] <- amUniqueResult$focalDatasetN
        profileResults[i, "samples"] <- amUniqueResult$comparisonDatasetN
        profileResults[i, "guessOptimum"] <- NA
    }
    profileResults <- profileResults[order(matchThreshold, decreasing=TRUE), ]
    
    if (guessOptimum) {
        
        LeftTrimZero <- function(trimString) {
            oldString <- trimString
            repeat {
                newString <- sub("^0", "", oldString)
                if (newString == oldString) break   
                oldString <- newString
            }
            return(newString)
        }
    
        ## Guess the second minimum...if it exists
        i <- 0
        repeat {
            i <- i + 1
            minimum <- which.min(profileResults$multipleMatch[i:nrow(profileResults)])[1] + (i-1)
            if (minimum == 1) {
                slopeAtMinimum <- 0
            }
            else {
                slopeAtMinimum <- profileResults$multipleMatch[minimum]-profileResults$multipleMatch[minimum-1]
            }
            if ((slopeAtMinimum < 0) || (minimum==nrow(profileResults))) break
        }
        secondMinimum <- minimum
        
        ## Guess the morphology
        checkMultipleMatch <- LeftTrimZero(paste(profileResults$multipleMatch, collapse=""))
        leadingZeroes <- nchar(paste(profileResults$multipleMatch, collapse="")) - nchar(checkMultipleMatch)
        if(leadingZeroes == length(profileResults$multipleMatch)) {
            profileMorphology <- "ZeroFlat"
        }
        else if (sum(profileResults$multipleMatch[(leadingZeroes+1):length(profileResults$multipleMatch)]==0)==0) {
            if (secondMinimum==length(profileResults$multipleMatch)) {
                profileMorphology <- "NoSecondMinimum"
            }
            else {
                profileMorphology <- "NonZeroSecondMinimum"
            }
        }
        else {
           profileMorphology <- "ZeroSecondMinimum"
        }
        
        ## Guess the optimal parameter value using a different approach depending on morphology
        if (profileMorphology %in% c("ZeroFlat", "NoSecondMinimum")) {
            guessOptimum <- which.min(abs(sapply(2:(length(profileResults$unique)-1), function(x) (profileResults$unique[x-1]-profileResults$unique[x+1])/2)))+1
        }
        else {
            guessOptimum <- secondMinimum
        }
        
        profileResults$missingDataLoad <- round(sum(amDatasetFocal$multilocus==amDatasetFocal$missingCode)/length(amDatasetFocal$multilocus),3)
        profileResults$allelicDiversity <- round(mean(unlist(lapply(amUniqueResult$alleleFreq$loci, function(x) x$numAlleles))), 1)
        profileResults$guessOptimum <- rep(FALSE, nrow(profileResults))
        profileResults$guessOptimum[guessOptimum] <- TRUE
        profileResults$guessMorphology <- profileMorphology

        
        if (verbose) {
            cat("allelematch:  missing data load for input dataset is", profileResults$missingDataLoad[1], "\n")
            cat("allelematch:  allelic diversity for input dataset is", profileResults$allelicDiversity[1], "\n")
            cat("allelematch:  Best guess for optimal parameter at alleleMismatch=",
                profileResults$alleleMismatch[guessOptimum],
                " OR matchThreshold=",
                profileResults$matchThreshold[guessOptimum],
                " OR cutHeight=",
                profileResults$cutHeight[guessOptimum], "\n", sep="")
            cat("allelematch:  Best guess for unique profile morphology: ", profileMorphology, "\n", sep="")
            if (profileMorphology %in% c("NonZeroSecondMinimum", "NoSecondMinimum")) {
                cat("allelematch:  Use extra caution.  Detection of optimal parameter is more error prone with this morphology.\n")
            }
        }
    }
    
    if (doPlot) {
        #dev.new()
        layout(matrix(c(1,1,1,1,1,1,1,2,2,2), 1, 10))
        par(mar=c(5.1, 4.1, 2, 2))

        plot.default(c(min(profileResults[, profileType]), max(profileResults[, profileType])), c(0, max(profileResults[, c("unique", "unclassified", "multipleMatch")])),
                     type="n", axes=TRUE, xlab=profileType, ylab="Count", cex=2)
        points(profileResults[, profileType], profileResults[, "unique"], pch=19, col="red")
        lines(profileResults[, profileType], profileResults[, "unique"], lwd=2, lty="solid", col="red")
        points(profileResults[, profileType], profileResults[, "multipleMatch"], pch=22, col="black")
        lines(profileResults[, profileType], profileResults[, "multipleMatch"], lwd=1, lty="solid", col="black")
        points(profileResults[, profileType], profileResults[, "unclassified"], pch=24, col="black")
        lines(profileResults[, profileType], profileResults[, "unclassified"],  lwd=1, lty="dotted", col="black")   
        if (guessOptimum) {
            arrows(profileResults[, profileType][which(profileResults$guessOptimum)], max(profileResults[, c("unique", "unclassified", "multipleMatch")])*0.3,
                   profileResults[, profileType][which(profileResults$guessOptimum)], 0, length=0.15, angle=20, lwd=2)
        }
        
        par(mar=c(0,0,0,1))
        plot.default(c(0,100), c(0,100), type="n", axes=FALSE, ylab="", xlab="")
        text(0, 100, "allelematch", cex=2, adj=c(0,0.5))
        text(0, 97, "amUniqueProfile()", cex=1, adj=c(0,0.5))
        text(0, 88, paste("missingDataLoad=", profileResults$missingDataLoad[1], sep=""), cex=1, adj=c(0,0.5))
        text(0, 86, paste("allelicDiversity=", profileResults$allelicDiversity[1], sep=""), cex=1, adj=c(0,0.5))
        legend(x=0, y=50, legend=c("unique", "multipleMatch", "unclassified"), lwd=c(2, 1, 1), lty=c("solid", "solid", "dotted"), col=c("red", "black", "black"), pch=c(19, 22, 24))
        if (guessOptimum) {
            text(0, 35, "Best guess for optimum:", cex=1, adj=c(0,0.5))
            text(0, 33, paste(profileType, "=", signif(profileResults[, profileType][which(profileResults$guessOptimum)], 2), sep=""), cex=1, adj=c(0, 0.5))
            text(0, 25, "Profile morphology:", cex=1, adj=c(0,0.5))
            text(0, 23, profileMorphology, cex=1, adj=c(0, 0.5))
            if (profileMorphology %in% c("NonZeroSecondMinimum", "NoSecondMinimum")) {
                text(0, 21, "Caution with optimum", col="red", cex=1, adj=c(0, 0.5))
            }
        }  
    }    
        
    return(profileResults)
} 

##
## amCSSForHTML()
amCSSForHTML <- function() {
    return("
        html {
            height: 100%;
        }
        
        body {
            background-color: inherit; 
            color: inherit; 
            font-family: Verdana; 
            font-size: xx-small; 
            margin: 0; 
            height: 100%;
        
        }
        
        a:active {
            color: #CC0000; 
        }
        
        a:link {
            color: #CC0000; 
        }
        
        a:visited {
            color: #CC0000; 
        }
        
        .amMismatchAllele {
            background-color: #CC0000;
            color: white;
            font-weight: bold;
        }
        
        .amMissingAllele {
            background-color: #FFCCCC;
        }
        
        .amInterpolatedAllele {
            background-color: blue;
            color: white;
        }
        
        .amGrid {
            border-collapse: separate;
        }
        
        .amGridContent {
            padding: 0;	
            border: 1px solid #7EACB1; 			
        }
        

        .amGridUpperPanel, .amGridLowerPanel {
            padding: 3px;	
            border-left: 0;
            border-right: 0;	
            background-color: #F4FAFB; 
            color: #2A769D;	 
            font-family: Verdana; 
            font-size: xx-small; 	
        }
        
        .amGridUpperPanel {
            border-top: 0px;
            border-bottom: 1px solid;
            border-color: #7EACB1; 
        }
        
        .amGridMiddlePanel {
            border: 0;	
        }
        
        .amGridLowerPanel {
            border-top: 1px solid;
            border-bottom: 0px; 
            border-color: #C2D4DA; 
        }
        
        .amGridUpperPanel td, .amGridLowerPanel td {
            color: #2A769D;	 
            font-family: Verdana; 
            font-size: xx-small; 		
        }
        
        
        .amTable {
            border: 0;
            border-spacing: 0;
            border-collapse: collapse;
            empty-cells: show;
            width: 100%;
            font-family: Verdana; 
            font-size: xx-small; 			
        }
        
        .amTableSeparate {	
            border-collapse: separate;		
        }
        
        .amTable td {
            padding: 3px; 
            border-bottom: 1px solid; 
            border-top: 0px;
            border-left: 0px;
            border-right: 1px solid; 
            border-color: #C2D4DA;  
            white-space:nowrap;
        }
        
            
        .amTable .amTableHeader, .amTable .amTableHeader td {
            background-color: #B7D8DC;	
            color: #000000; 
            border-bottom: 1px solid; 
            border-right: 1px solid; 
            border-color: #7EACB1; 
            background-repeat: repeat-x;		
            vertical-align: top;
            white-space:nowrap;
        }
        
        .amPointer {
            cursor: pointer;
        }
        
        
        .amTableHeaderBtn {
            width: 100%;
            font-family: Verdana; 
            font-size: xx-small; 		
        }
        
        .amTableHeader .amTableHeaderBtn td {
            background: transparent;
            padding: 0;
            border: 0;
            white-space: nowrap;		
        }
        
        .amTableSelectRow {
            background-color: #FFFF66; 
            color: #000000;
        }\n")
}




##
## amHTML.amPairwise()
amHTML.amPairwise <- function(x, htmlFile=NULL, htmlCSS=amCSSForHTML()) {
    
    if (class(x) != "amPairwise") {
        stop("allelematch:  this function requires an \"amPairwise\" object", call.=FALSE)
    }
    
    
    if (is.null(htmlFile)) {   
        if (file.exists(Sys.getenv("TMP"))) {
            htmlFilePath <- Sys.getenv("TMP")
        }
        else if (file.exists(Sys.getenv("R_USER"))) {
            htmlFilePath <- Sys.getenv("R_USER")
        }
        else {
            htmlFilePath <- "."
        }
        htmlFile <- paste(htmlFilePath, "\\", class(x), "_", paste(sample(c(LETTERS[1:26], 0:9))[1:8], collapse=""), ".htm", sep="")
        
        usingTmpFile <- TRUE
    }
    else {
        usingTmpFile <- FALSE
    }
    
    fileID <- file(htmlFile, open="wt")
    
    ## Open page with correct HTML header
    cat("<!DOCTYPE html PUBLIC \"-//W3C//DTD XHTML 1.0 Transitional//EN\" \"http://www.w3.org/TR/xhtml1/DTD/xhtml1-transitional.dtd\">\n
         <html xmlns='http://www.w3.org/1999/xhtml' xml:lang='en' lang='en'><head>\n
        <meta http-equiv=\"Content-Type\" content=\"text/html; charset=utf-8\"></meta>\n", file=fileID)
    
    ## Title of page
    cat(paste(c("<title>allelematch: ", class(x), "() output</title>"), collapse=""), file=fileID, append=TRUE)
    
    ## Add CSS
    cat(paste(c("<style type=\"text/css\">", htmlCSS, "</style>"), collapse="\n"), file=fileID, append=TRUE)
    cat("</head><body>", file=fileID, append=TRUE)
    
    ## Page-wide DIV
    cat("<div style=\"margin-left:5%; margin-right:5%\">", file=fileID, append=TRUE)

    ## Page title area
    cat("<br><table cellspacing=\"0\" class=\"amGrid\"><tr><td class=\"amGridContent\">", file=fileID, append=TRUE)
    cat("<div class=\"amGridUpperPanel\">", file=fileID, append=TRUE)
    cat("<div style=\"font-size:x-small;\"><table>",file=fileID, append=TRUE)
    cat("<tr><td style=\"width:400px;\"><span style=\"font-size:20px;\">allelematch<br>pairwise analysis</span><br><br></td></tr>\n", file=fileID, append=TRUE)       
    
    headerHTML <-  matrix("", 6, 2)
    headerHTML[1,] <- c("\nfocal dataset N=", x$focalDatasetN)
        
    if (class(x) == "amPairwise") {
        if (x$focalIsComparison)
            headerHTML[2, ]<- c("focal dataset compared against itself", "")
        else
            headerHTML[2, ] <- c("comparison dataset N=", x$comparisonDatasetN)
        headerHTML[3, ] <- c("missing data represented by: ", x$missingCode)
        headerHTML[4, ] <- c("alleleMismatch (m-hat; maximum number of mismatching alleles): ", round(x$alleleMismatch, 2))
        headerHTML[5, ] <- c("matchThreshold (s-hat; lowest matching score returned): ", round(x$matchThreshold, 3))
        headerHTML[6, ] <- c("summary generated: ", date())
    }
    
    cat(apply(headerHTML, 1, function(x) paste("<tr><td><b>", x[1], "</b><em>", x[2], "</em></td></tr>\n",sep="")), file=fileID, append=TRUE)
    cat("</table></div></div></td></tr></table><br><br>\n", file=fileID, append=TRUE)

    y <- x$pairwise
    ## Do tables for each focal genotype
    for (i in 1:length(y)) {

        ## Table header
        cat("<table cellspacing=\"0\" class=\"amGrid\"><tr><td class=\"amGridContent\">", file=fileID, append=TRUE)
               
        ## Content above header names
        cat("<div class=\"amGridUpperPanel\">", file=fileID, append=TRUE)
        cat(paste(c("<span style=\"font-size:medium;\">(", i, " of ", length(y), ")</span></div>"), collapse=""), file=fileID, append=TRUE)
        
       
        ## Header names
        cat("<div class=\"amGridMiddlePanel\">\n", file=fileID, append=TRUE)
        cat("<table cellspacing=\"0\" class=\"amTable amTableSeparate\">\n", file=fileID, append=TRUE)
        cat("<tr class=\"amTableHeader\">\n", file=fileID, append=TRUE)
        if (!is.null(y[[i]]$focal$metaData))
            headerNames <- c("", "", dimnames(y[[i]]$focal$multilocus)[[2]], "Score")
        else
            headerNames <- c("", dimnames(y[[i]]$focal$multilocus)[[2]], "Score")
        sapply(headerNames, function(j)
                  cat(paste("<td class=\"amPointer\"><table cellspacing=\"0\" class=\"amTableHeaderBtn\"><tr><td>",
                         j,
                         "</td><td style=\"width:10px;\">&nbsp;</td></tr></table></td>\n", sep=""), file=fileID, append=TRUE))
        cat("</tr>\n", file=fileID, append=TRUE)
        
        ## Content of tables
        cat("<tr>\n", file=fileID, append=TRUE)
    
        ## Focal genotype row
        if (!is.null(y[[i]]$focal$metaData))
            focalRow <- c(y[[i]]$focal$index, y[[i]]$focal$metaData, y[[i]]$focal$multilocus, "FOCAL")
        else
            focalRow <- c(y[[i]]$focal$index, y[[i]]$focal$multilocus, "FOCAL")
        sapply(focalRow, function(row) cat(paste("<td class=\"amTableSelectRow\"><div>",
                                               row,
                                               "</div></td>\n", sep=""), file=fileID, append=TRUE))
        cat("</tr>\n", file=fileID, append=TRUE)
        
        ## Matching genotype rows
        for (k in 1:nrow(y[[i]]$match$multilocus)) {
            if (!is.null(y[[i]]$focal$metaData)) {
                matchRow <- c(y[[i]]$match$index[k], y[[i]]$match$metaData[k], y[[i]]$match$multilocus[k, ], y[[i]]$match$score[k])
                matchFlags <- c(1, 1, y[[i]]$match$flags[k, ], 1)
            }
            else {
                matchRow <- c(y[[i]]$match$index[k],  y[[i]]$match$multilocus[k, ], y[[i]]$match$score[k])
                matchFlags <- c(1, y[[i]]$match$flags[k, ], 1)
            }
            matchRowHTML <- "<tr>"
            for (m in 1:length(matchRow)) {
                if (matchFlags[m] == 1) {
                    matchRowHTML <- paste(matchRowHTML, "<td><div>", matchRow[m], "</div></td>", sep="")
                }
                else if (matchFlags[m] == 0) {
                    matchRowHTML <- paste(matchRowHTML, "<td><div class=\"amMismatchAllele\">", matchRow[m], "</div></td>", sep="")
                }
                else if(matchFlags[m] == 2) {
                    matchRowHTML <- paste(matchRowHTML, "<td><div class=\"amMissingAllele\">", matchRow[m], "</div></td>", sep="") 
                }
            }
            matchRowHTML <- paste(matchRowHTML, "</tr>\n")
            cat(matchRowHTML, file=fileID, append=TRUE) 
        }

 
        ## End of content tables
        cat("</table></div>", file=fileID, append=TRUE)
        
        ## Content below header names
        cat("<div class=\"amGridLowerPanel\">", file=fileID, append=TRUE)
        
        if ((y[[i]]$match$perfect==0) && (y[[i]]$match$partial==0)) {
            cat("<span>No matches found.</span>", file=fileID, append=TRUE)
        }
        else {
          cat(paste(c("<span>", (y[[i]]$match$perfect), " perfect matches found.  ", y[[i]]$match$partial, " partial matches found.</span>"), collapse=""), file=fileID, append=TRUE)
        }
        cat("</div>", file=fileID, append=TRUE)
        cat("</td></tr></table>", file=fileID, append=TRUE)
        cat("<br><br>", file=fileID, append=TRUE)
    }
    
    
    ## Page-wide
    cat("<span style=\"font-size:x-small;\">Generated by allelematch:  an R package<br></span>", file=fileID, append=TRUE)
    cat("<span style=\"font-size:x-small;\">To reference this analysis please use citation(\"allelematch\")<br><br></span>", file=fileID, append=TRUE)
    cat("</div>", file=fileID, append=TRUE)
    cat("</body></html>", file=fileID, append=TRUE)
    close(fileID)
    
    ## If htmlFile is not given, assume that this is being run
    ## for a quick view, therefore open produced html file in browser
    ## Also take this opportunity to clean up temp folder, if required
    if (usingTmpFile) {
        cat("Opening HTML file (", htmlFile, ") in default browser...\n", sep="")
        browseURL(htmlFile)
        if (htmlFilePath != ".") {
            oldTmpFiles <- Sys.glob(paste(htmlFilePath, "\\am*.htm", sep=""))
            if (length(oldTmpFiles) > 0) {
                deleteFiles <- file.remove(oldTmpFiles[sapply(file.info(Sys.glob(paste(htmlFilePath, "\\am*.htm", sep="")))$ctime, function(x) difftime(Sys.time(), x, units="hours")>24)])
                if (sum(deleteFiles) > 0) cat("Deleting allelematch HTML output older than 24 hours from temporary folder...\n")
            }
        }
    }
}

##
## amHTML.amCluster()
amHTML.amCluster <- function(x, htmlFile=NULL, htmlCSS=amCSSForHTML()) {
    
    if (class(x) != "amCluster") {
        stop("allelematch:  this function requires an \"amCluster\" object", call.=FALSE)
    }
    
    if (is.null(htmlFile)) {   
        if (file.exists(Sys.getenv("TMP"))) {
            htmlFilePath <- Sys.getenv("TMP")
        }
        else if (file.exists(Sys.getenv("R_USER"))) {
            htmlFilePath <- Sys.getenv("R_USER")
        }
        else {
            htmlFilePath <- "."
        }
        htmlFile <- paste(htmlFilePath, "\\", class(x), "_", paste(sample(c(LETTERS[1:26], 0:9))[1:8], collapse=""), ".htm", sep="")
        
        usingTmpFile <- TRUE
    }
    else {
        usingTmpFile <- FALSE
    }
    
    fileID <- file(htmlFile, open="wt")
    
    ## Open page with correct HTML header
    cat("<!DOCTYPE html PUBLIC \"-//W3C//DTD XHTML 1.0 Transitional//EN\" \"http://www.w3.org/TR/xhtml1/DTD/xhtml1-transitional.dtd\">\n
         <html xmlns='http://www.w3.org/1999/xhtml' xml:lang='en' lang='en'><head>\n
        <meta http-equiv=\"Content-Type\" content=\"text/html; charset=utf-8\"></meta>\n", file=fileID)
    
    ## Title of page
    cat(paste(c("<title>allelematch: ", class(x), "() output</title>"), collapse=""), file=fileID, append=TRUE)
    
    ## Add CSS
    cat(paste(c("<style type=\"text/css\">", htmlCSS, "</style>"), collapse="\n"), file=fileID, append=TRUE)
    cat("</head><body>", file=fileID, append=TRUE)
    
    ## Page-wide DIV
    cat("<div style=\"margin-left:5%; margin-right:5%\">", file=fileID, append=TRUE)

    ## Page title area
    cat("<br><table cellspacing=\"0\" class=\"amGrid\"><tr><td class=\"amGridContent\">", file=fileID, append=TRUE)
    cat("<div class=\"amGridUpperPanel\">", file=fileID, append=TRUE)
    cat("<div style=\"font-size:x-small;\"><table>",file=fileID, append=TRUE)
    cat("<tr><td style=\"width:400px;\"><span style=\"font-size:20px;\">allelematch<br>cluster analysis</span><br><br></td></tr>\n", file=fileID, append=TRUE)       
    
    headerHTML <-  matrix("", 12, 2)
    headerHTML[1,] <- c("\nFocal dataset N=", x$focalDatasetN)
    
    if (class(x) == "amCluster") {
        headerHTML[2, ] <- c("unique N=", nrow(x$unique$multilocus))
        headerHTML[3, ] <- c("unique (consensus) N=", length(x$cluster))
        headerHTML[4, ] <- c("unique (singletons) N=", length(x$singletons))
        headerHTML[5, ] <- c("missing data represented by: ", x$missingCode)
        headerHTML[6, ] <- c("missing data matching method: ", x$missingMethod)
        headerHTML[7, ] <- c("clustered genotypes consensus method: ", x$consensusMethod)
        headerHTML[8, ] <- c("hierarchical clustering method: ", x$clusterMethod)
        headerHTML[9, ] <- c("cutHeight (d-hat; dynamic tree cutting height): ", format(x$cutHeight, scientific=FALSE))
        headerHTML[10, ] <- c("run until only singletons: ", x$runUntilSingletons)
        headerHTML[11, ] <- c("runs: ", x$totalRuns)
        headerHTML[12, ] <- c("summary generated: ", date())
    }
    
    cat(apply(headerHTML, 1, function(x) paste("<tr><td><b>", x[1], "</b><em>", x[2], "</em></td></tr>\n",sep="")), file=fileID, append=TRUE)
    cat("</table></div></div></td></tr></table><br><br>\n", file=fileID, append=TRUE)

    ## Unique Genotypes
    y <- x$unique
    ## Table header
    cat("<table cellspacing=\"0\" class=\"amGrid\"><tr><td class=\"amGridContent\">", file=fileID, append=TRUE)
           
    ## Content above header names
    cat("<div class=\"amGridUpperPanel\">", file=fileID, append=TRUE)
    cat(paste(c("<span style=\"font-size:medium;\">Unique genotypes</span></div>"), collapse=""), file=fileID, append=TRUE)
    
    ## Header names
    cat("<div class=\"amGridMiddlePanel\">\n", file=fileID, append=TRUE)
    cat("<table cellspacing=\"0\" class=\"amTable amTableSeparate\">\n", file=fileID, append=TRUE)
    cat("<tr class=\"amTableHeader\">\n", file=fileID, append=TRUE)
    if (!is.null(y$metaData))
        headerNames <- c("", "", "", dimnames(y$multilocus)[[2]], "Type")
    else
        headerNames <- c("", "", dimnames(y$multilocus)[[2]], "Type")
    sapply(headerNames, function(j)
              cat(paste("<td class=\"amPointer\"><table cellspacing=\"0\" class=\"amTableHeaderBtn\"><tr><td>",
                     j,
                     "</td><td style=\"width:10px;\">&nbsp;</td></tr></table></td>\n", sep=""), file=fileID, append=TRUE))
    cat("</tr>\n", file=fileID, append=TRUE)
    
    ## Unique genotype rows
    for (k in 1:nrow(y$multilocus)) {
        cat("<tr onmouseout=\"this.style.cssText='background-color:none;';\" onmouseover=\"this.style.cssText='background-color:#E0FFFF';\" style=\"\">\n", file=fileID, append=TRUE)
        if (!is.null(y$metaData)) {
            matchRow <- c(paste("<a href=\"#", y$index[k], "\">Jump</a>", sep=""), paste("<a name=\"back_", y$index[k], "\">", y$index[k], "</a>", sep=""), y$metaData[k], y$multilocus[k, ], y$uniqueType[k])
        }
        else {
            matchRow <- c(paste("<a href=\"#", y$index[k], "\">Jump</a>", sep=""), paste("<a name=\"back_", y$index[k], "\">", y$index[k], "</a>", sep=""), y$multilocus[k, ], y$uniqueType[k])
        }
        sapply(matchRow, function(j)
              cat(paste("<td><div>",
                     j,
                     "</div></td>\n", sep=""), file=fileID, append=TRUE))
        cat("</tr>\n", file=fileID, append=TRUE)
    }

    ## End of content tables
    cat("</table></div>", file=fileID, append=TRUE)
    
    ## Content below header names
    cat("<div class=\"amGridLowerPanel\">", file=fileID, append=TRUE)
    cat(paste(c("<span>There were ", nrow(y$multilocus), " unique genotypes identified using the parameters supplied.</span>"), collapse=""), file=fileID, append=TRUE)
    cat("</div>", file=fileID, append=TRUE)
    cat("</td></tr></table>", file=fileID, append=TRUE)
    cat("<br><br>", file=fileID, append=TRUE)


    ## Clustered genotypes
    y <- x$cluster
    if (length(y) > 0) {
        ## Do tables for each cluster
        for (i in 1:length(y)) {
    
            ## Table header
            cat("<table cellspacing=\"0\" class=\"amGrid\"><tr><td class=\"amGridContent\">", file=fileID, append=TRUE)
                   
            ## Content above header names
            cat("<div class=\"amGridUpperPanel\">", file=fileID, append=TRUE)
            cat(paste(c("<a name=\"", y[[i]]$focal$index, "\"><span style=\"font-size:medium;\">(Cluster ", i, " of ", length(y), ")</span></a>",
                        "<br><a href=\"#back_", y[[i]]$focal$index, "\">Jump up</a><br></div>"), collapse=""), file=fileID, append=TRUE)
            
           
            ## Header names
            cat("<div class=\"amGridMiddlePanel\">\n", file=fileID, append=TRUE)
            cat("<table cellspacing=\"0\" class=\"amTable amTableSeparate\">\n", file=fileID, append=TRUE)
            cat("<tr class=\"amTableHeader\">\n", file=fileID, append=TRUE)
            if (!is.null(y[[i]]$focal$metaData))
                headerNames <- c("", "", "", dimnames(y[[i]]$focal$multilocus)[[2]], "Score")
            else
                headerNames <- c("", "", dimnames(y[[i]]$focal$multilocus)[[2]], "Score")
            sapply(headerNames, function(j)
                      cat(paste("<td class=\"amPointer\"><table cellspacing=\"0\" class=\"amTableHeaderBtn\"><tr><td>",
                             j,
                             "</td><td style=\"width:10px;\">&nbsp;</td></tr></table></td>\n", sep=""), file=fileID, append=TRUE))
            cat("</tr>\n", file=fileID, append=TRUE)
            
            ## Content of tables
            cat("<tr>\n", file=fileID, append=TRUE)
        
            ## Focal genotype row
            if (!is.null(y[[i]]$focal$metaData)) {
                focalRow <- c("CONSENSUS", y[[i]]$focal$index, y[[i]]$focal$metaData, y[[i]]$focal$multilocus, "")
                focalFlags <- c(1, 1, 1, y[[i]]$focal$flags[1, ], 1)
            }
            else {
                focalRow <- c("CONSENSUS", y[[i]]$focal$index, y[[i]]$focal$multilocus, "")
                focalFlags <- c(1, 1, y[[i]]$focal$flags[1, ], 1)
            }
            focalRowHTML <- "<tr>"
            for (m in 1:length(focalRow)) {
                    if (focalFlags[m] == 3) {
                        focalRowHTML <- paste(focalRowHTML, "<td class=\"amTableSelectRow\"><div class=\"amInterpolatedAllele\">", focalRow[m], "</div></td>", sep="")
                    }
                    else {
                        focalRowHTML <- paste(focalRowHTML, "<td class=\"amTableSelectRow\"><div>", focalRow[m], "</div></td>", sep="")
                    }
                }
            focalRowHTML <- paste(focalRowHTML, "</tr>\n")
            cat(focalRowHTML, file=fileID, append=TRUE) 
            
    
            ## Matching genotype rows
            for (k in 1:nrow(y[[i]]$match$multilocus)) {
                if (!is.null(y[[i]]$focal$metaData)) {
                    matchRow <- c("MEMBER", y[[i]]$match$index[k], y[[i]]$match$metaData[k], y[[i]]$match$multilocus[k, ], y[[i]]$match$score[k])
                    matchFlags <- c(1, 1, 1, y[[i]]$match$flags[k, ], 1)
                }
                else {
                    matchRow <- c("MEMBER", y[[i]]$match$index[k],  y[[i]]$match$multilocus[k, ], y[[i]]$match$score[k])
                    matchFlags <- c(1, 1, y[[i]]$match$flags[k, ], 1)
                }
                matchRowHTML <- "<tr>"
                for (m in 1:length(matchRow)) {
                    if (matchFlags[m] == 1) {
                        matchRowHTML <- paste(matchRowHTML, "<td><div>", matchRow[m], "</div></td>", sep="")
                    }
                    else if (matchFlags[m] == 0) {
                        matchRowHTML <- paste(matchRowHTML, "<td><div class=\"amMismatchAllele\">", matchRow[m], "</div></td>", sep="")
                    }
                    else if(matchFlags[m] == 2) {
                        matchRowHTML <- paste(matchRowHTML, "<td><div class=\"amMissingAllele\">", matchRow[m], "</div></td>", sep="") 
                    }
                }
                matchRowHTML <- paste(matchRowHTML, "</tr>\n")
                cat(matchRowHTML, file=fileID, append=TRUE)
    
            }
    
     
            ## End of content tables
            cat("</table></div>", file=fileID, append=TRUE)
            
            ## Content below header names
            cat("<div class=\"amGridLowerPanel\">", file=fileID, append=TRUE)
            
            if ((y[[i]]$match$perfect==0) && (y[[i]]$match$partial==0)) {
                cat("<span>No matches found.</span>", file=fileID, append=TRUE)
            }
            else {
              cat(paste(c("<span>", (y[[i]]$match$perfect), " perfect matches found.  ", y[[i]]$match$partial, " partial matches found.</span>"), collapse=""), file=fileID, append=TRUE)
            }
            cat("</div>", file=fileID, append=TRUE)
            cat("</td></tr></table>", file=fileID, append=TRUE)
            cat("<br><br>", file=fileID, append=TRUE)
        }
    }

    ## Singleton genotypes
    y <- x$singletons
    if (length(y) > 0) {
        ## Do tables for each cluster
        for (i in 1:length(y)) {
    
            ## Table header
            cat("<table cellspacing=\"0\" class=\"amGrid\"><tr><td class=\"amGridContent\">", file=fileID, append=TRUE)
                   
            ## Content above header names
            cat("<div class=\"amGridUpperPanel\">", file=fileID, append=TRUE)
            cat(paste(c("<a name=\"", y[[i]]$focal$index, "\"><span style=\"font-size:medium;\">(Singleton ", i, " of ", length(y), ")</span></a>",
                        "<br><a href=\"#back_", y[[i]]$focal$index, "\">Jump up</a><br></div>"), collapse=""), file=fileID, append=TRUE)
            
           
            ## Header names
            cat("<div class=\"amGridMiddlePanel\">\n", file=fileID, append=TRUE)
            cat("<table cellspacing=\"0\" class=\"amTable amTableSeparate\">\n", file=fileID, append=TRUE)
            cat("<tr class=\"amTableHeader\">\n", file=fileID, append=TRUE)
            if (!is.null(y[[i]]$focal$metaData))
                headerNames <- c("", "", "", dimnames(y[[i]]$focal$multilocus)[[2]], "Score")
            else
                headerNames <- c("", "", dimnames(y[[i]]$focal$multilocus)[[2]], "Score")
            sapply(headerNames, function(j)
                      cat(paste("<td class=\"amPointer\"><table cellspacing=\"0\" class=\"amTableHeaderBtn\"><tr><td>",
                             j,
                             "</td><td style=\"width:10px;\">&nbsp;</td></tr></table></td>\n", sep=""), file=fileID, append=TRUE))
            cat("</tr>\n", file=fileID, append=TRUE)
            
            ## Content of tables
            cat("<tr>\n", file=fileID, append=TRUE)
        
            ## Focal genotype row
            if (!is.null(y[[i]]$focal$metaData)) {
                focalRow <- c("SINGLETON", y[[i]]$focal$index, y[[i]]$focal$metaData, y[[i]]$focal$multilocus, "")
                focalFlags <- c(1, 1, 1, y[[i]]$focal$flags[1, ], 1)
            }
            else {
                focalRow <- c("SINGLETON", y[[i]]$focal$index, y[[i]]$focal$multilocus, "")
                focalFlags <- c(1, 1, y[[i]]$focal$flags[1, ], 1)
            }
            focalRowHTML <- "<tr>"
            for (m in 1:length(focalRow)) {
                    if (focalFlags[m] == 3) {
                        focalRowHTML <- paste(focalRowHTML, "<td class=\"amTableSelectRow\"><div class=\"amInterpolateAllele\">", focalRow[m], "</div></td>", sep="")
                    }
                    else {
                        focalRowHTML <- paste(focalRowHTML, "<td class=\"amTableSelectRow\"><div>", focalRow[m], "</div></td>", sep="")
                    }
                }
            focalRowHTML <- paste(focalRowHTML, "</tr>\n")
            cat(focalRowHTML, file=fileID, append=TRUE) 
            
            
            ## Matching genotype rows
            for (k in 1:nrow(y[[i]]$match$multilocus)) {
                if (!is.null(y[[i]]$focal$metaData)) {
                    matchRow <- c("CLOSEST", y[[i]]$match$index[k], y[[i]]$match$metaData[k], y[[i]]$match$multilocus[k, ], y[[i]]$match$score[k])
                    matchFlags <- c(1, 1, 1, y[[i]]$match$flags[k, ], 1)
                }
                else {
                    matchRow <- c("CLOSEST", y[[i]]$match$index[k],  y[[i]]$match$multilocus[k, ], y[[i]]$match$score[k])
                    matchFlags <- c(1, 1, y[[i]]$match$flags[k, ], 1)
                }
                matchRowHTML <- "<tr>"
                for (m in 1:length(matchRow)) {
                    if (matchFlags[m] == 1) {
                        matchRowHTML <- paste(matchRowHTML, "<td><div>", matchRow[m], "</div></td>", sep="")
                    }
                    else if (matchFlags[m] == 0) {
                        matchRowHTML <- paste(matchRowHTML, "<td><div class=\"amMismatchAllele\">", matchRow[m], "</div></td>", sep="")
                    }
                    else if(matchFlags[m] == 2) {
                        matchRowHTML <- paste(matchRowHTML, "<td><div class=\"amMissingAllele\">", matchRow[m], "</div></td>", sep="") 
                    }
                }
                matchRowHTML <- paste(matchRowHTML, "</tr>\n")
                cat(matchRowHTML, file=fileID, append=TRUE) 
            }
    
     
            ## End of content tables
            cat("</table></div>", file=fileID, append=TRUE)
            
            ## Content below header names
            cat("<div class=\"amGridLowerPanel\">", file=fileID, append=TRUE)
            
            cat(paste(c("<span>Closest match to singleton shown for diagnostic purposes</span>"), collapse=""), file=fileID, append=TRUE)
            
            cat("</div>", file=fileID, append=TRUE)
            cat("</td></tr></table>", file=fileID, append=TRUE)
            cat("<br><br>", file=fileID, append=TRUE)
        }
    }
    
    
    ## Page-wide
    cat("<span style=\"font-size:x-small;\">Generated by allelematch:  an R package<br></span>", file=fileID, append=TRUE)
    cat("<span style=\"font-size:x-small;\">To reference this analysis please use citation(\"allelematch\")<br><br></span>", file=fileID, append=TRUE)
    cat("</div>", file=fileID, append=TRUE)
    cat("</body></html>", file=fileID, append=TRUE)
    close(fileID)
    
    ## If htmlFile is not given, assume that this is being run
    ## for a quick view, therefore open produced html file in browser
    ## Also take this opportunity to clean up temp folder, if required
    if (usingTmpFile) {
        cat("Opening HTML file (", htmlFile, ") in default browser...\n", sep="")
        browseURL(htmlFile)
        if (htmlFilePath != ".") {
            oldTmpFiles <- Sys.glob(paste(htmlFilePath, "\\am*.htm", sep=""))
            if (length(oldTmpFiles) > 0) {
                deleteFiles <- file.remove(oldTmpFiles[sapply(file.info(Sys.glob(paste(htmlFilePath, "\\am*.htm", sep="")))$ctime, function(x) difftime(Sys.time(), x, units="hours")>24)])
                if (sum(deleteFiles) > 0) cat("Deleting allelematch HTML output older than 24 hours from temporary folder...\n")
            }
        }
    }
}    
    
##
## amHTML.amUnique()
amHTML.amUnique <- function(x, htmlFile=NULL, htmlCSS=amCSSForHTML()) {

  
    if (class(x) != "amUnique") {
        stop("allelematch:  this function requires an \"amUnique\" object", call.=FALSE)
    }
    
    if (is.null(htmlFile)) {   
        if (file.exists(Sys.getenv("TMP"))) {
            htmlFilePath <- Sys.getenv("TMP")
        }
        else if (file.exists(Sys.getenv("R_USER"))) {
            htmlFilePath <- Sys.getenv("R_USER")
        }
        else {
            htmlFilePath <- "."
        }
        htmlFile <- paste(htmlFilePath, "\\", class(x), "_", paste(sample(c(LETTERS[1:26], 0:9))[1:8], collapse=""), ".htm", sep="")
        
        usingTmpFile <- TRUE
    }
    else {
        usingTmpFile <- FALSE
    }
    
    fileID <- file(htmlFile, open="wt")
    
    ## Open page with correct HTML header
    cat("<!DOCTYPE html PUBLIC \"-//W3C//DTD XHTML 1.0 Transitional//EN\" \"http://www.w3.org/TR/xhtml1/DTD/xhtml1-transitional.dtd\">\n
         <html xmlns='http://www.w3.org/1999/xhtml' xml:lang='en' lang='en'><head>\n
        <meta http-equiv=\"Content-Type\" content=\"text/html; charset=utf-8\"></meta>\n", file=fileID)
    
    ## Title of page
    cat(paste(c("<title>allelematch: ", class(x), "() output</title>"), collapse=""), file=fileID, append=TRUE)
    
    ## Add CSS
    cat(paste(c("<style type=\"text/css\">", htmlCSS, "</style>"), collapse="\n"), file=fileID, append=TRUE)
    cat("</head><body>", file=fileID, append=TRUE)
    
    ## Page-wide DIV
    cat("<div style=\"margin-left:5%; margin-right:5%\">", file=fileID, append=TRUE)

    ## Page title area
    cat("<br><table cellspacing=\"0\" class=\"amGrid\"><tr><td class=\"amGridContent\">", file=fileID, append=TRUE)
    cat("<div class=\"amGridUpperPanel\">", file=fileID, append=TRUE)
    cat("<div style=\"font-size:x-small;\"><table>",file=fileID, append=TRUE)
    cat("<tr><td style=\"width:400px;\"><span style=\"font-size:20px;\">allelematch<br>unique analysis</span><br><br></td></tr>\n", file=fileID, append=TRUE)       
    

        
    if (class(x) == "amUnique") {
        headerHTML <-  matrix("", 14, 2)
        headerHTML[1, ] <- c("\nunique N=", x$focalDatasetN)
        headerHTML[2, ] <- c("samples N=", x$comparisonDatasetN)
        headerHTML[3, ] <- c("<br>loci N=", length(x$alleleFreq$loci))
        headerHTML[4, ] <- c("locus names: ", paste(sapply(x$alleleFreq$loci, function(x) x$name), collapse=", "))
        headerHTML[5, ] <- c("<br>missing data represented by: ", x$missingCode)
        headerHTML[6, ] <- c("clustered genotypes consensus method: ", x$consensusMethod)
        headerHTML[7, ] <- c("Psib calculated for: ", ifelse(x$doPsib!="all", "samples with no mismatches", "all samples (mismatches treated as missing data)"))      
        headerHTML[8, ] <- c("<br>alleleMismatch (m-hat; maximum number of mismatching alleles): ", round(x$alleleMismatch, 2))
        headerHTML[9, ] <- c("cutHeight (d-hat; dynamic tree cutting height): ", round(x$cutHeight, 3))
        headerHTML[10, ] <- c("matchThreshold (s-hat; lowest matching score returned): ", round(x$matchThreshold, 3))
        if (x$numUnclassified > 0) {
            headerHTML[11, ] <- c("<span style=\"color:red;\"><br>unclassified (samples that were not classified) N=</span>", paste("<span style=\"color:red;\">", x$numUnclassified, "</span>", sep=""))
        }
        else {
            headerHTML[11, ] <- c("<br>unclassified (samples that were not classified) N=", x$numUnclassified)
        }
        headerHTML[12, ] <- c("multipleMatch (samples that match more than one unique genotype) N=", x$numMultipleMatches)
        headerHTML[13, ] <- c("<br>Note: Unique genotypes are determined based on clustering of their scores followed by a dynamic tree cutting procedure (see supplementary documentation).",
                             "  Psib appears for reference purposes and is not used to determine unique genotypes.  It is calculated using allele frequencies in the unique genotype set.")
        headerHTML[14, ] <- c("<br>summary generated: ", date())
    }
    
    cat(apply(headerHTML, 1, function(x) paste("<tr><td><b>", x[1], "</b><em>", x[2], "</em></td></tr>\n",sep="")), file=fileID, append=TRUE)
    cat("</table></div></div></td></tr></table><br><br>\n", file=fileID, append=TRUE)


    ## Unique Genotypes
    y <- x$unique
    ## Table header
    cat("<table cellspacing=\"0\" class=\"amGrid\"><tr><td class=\"amGridContent\">", file=fileID, append=TRUE)
           
    ## Content above header names
    cat("<div class=\"amGridUpperPanel\">", file=fileID, append=TRUE)
    cat(paste(c("<span style=\"font-size:medium;\">Unique genotypes</span></div>"), collapse=""), file=fileID, append=TRUE)
    
    ## Header names
    cat("<div class=\"amGridMiddlePanel\">\n", file=fileID, append=TRUE)
    cat("<table cellspacing=\"0\" class=\"amTable amTableSeparate\">\n", file=fileID, append=TRUE)
    cat("<tr class=\"amTableHeader\">\n", file=fileID, append=TRUE)
    if (!is.null(y$metaData))
        headerNames <- c("", "", "", dimnames(y$multilocus)[[2]], "Psib", "Type")
    else
        headerNames <- c("", "", dimnames(y$multilocus)[[2]], "Psib", "Type")
    sapply(headerNames, function(j)
              cat(paste("<td class=\"amPointer\"><table cellspacing=\"0\" class=\"amTableHeaderBtn\"><tr><td>",
                     j,
                     "</td><td style=\"width:10px;\">&nbsp;</td></tr></table></td>\n", sep=""), file=fileID, append=TRUE))
    cat("</tr>\n", file=fileID, append=TRUE)
    
    ## Unique genotype rows
    for (k in 1:nrow(y$multilocus)) {
        
        psib <- signif(y$psib[k],2)
        if (is.na(psib)) {
            psib <- "&nbsp; &nbsp; ---"
        }
        else if (psib < 0.00001) {
            psib <- "<0.00001"
        }
        else {
            psib <- format(psib, scientific=FALSE)
        }
        cat("<tr onmouseout=\"this.style.cssText='background-color:none;';\" onmouseover=\"this.style.cssText='background-color:#E0FFFF';\" style=\"\">\n", file=fileID, append=TRUE)
        if (!is.null(y$metaData)) {
            matchRow <- c(paste("<a href=\"#", y$index[k], "\">Jump</a>", sep=""), y$index[k], y$metaData[k], y$multilocus[k, ], psib, y$rowFlag[k])
        }
        else {
            matchRow <- c(paste("<a href=\"#", y$index[k], "\">Jump</a>", sep=""), y$index[k], y$multilocus[k, ], psib, y$rowFlag[k])
        }
        sapply(matchRow, function(j)
              cat(paste("<td><div>",
                     j,
                     "</div></td>\n", sep=""), file=fileID, append=TRUE))
        cat("</tr>\n", file=fileID, append=TRUE)
    }

    ## End of content tables
    cat("</table></div>", file=fileID, append=TRUE)
    
    ## Content below header names
    cat("<div class=\"amGridLowerPanel\">", file=fileID, append=TRUE)
    cat(paste(c("<span>There were ", nrow(y$multilocus), " unique genotypes identified using the parameters supplied.</span>"), collapse=""), file=fileID, append=TRUE)
    cat("</div>", file=fileID, append=TRUE)
    cat("</td></tr></table>", file=fileID, append=TRUE)
    cat("<br><br>", file=fileID, append=TRUE)


    ## Unclassified genotypes
    if (!is.null(x$unclassified)) {
        y <- x$unclassified
        ## Table header
        cat("<table cellspacing=\"0\" class=\"amGrid\"><tr><td class=\"amGridContent\">", file=fileID, append=TRUE)
               
        ## Content above header names
        cat("<div class=\"amGridUpperPanel\">", file=fileID, append=TRUE)
        cat(paste(c("<span style=\"font-size:medium; color:red; font-weight:bold;\">Unclassified samples</span></div>"), collapse=""), file=fileID, append=TRUE)
        
        ## Header names
        cat("<div class=\"amGridMiddlePanel\">\n", file=fileID, append=TRUE)
        cat("<table cellspacing=\"0\" class=\"amTable amTableSeparate\">\n", file=fileID, append=TRUE)
        cat("<tr class=\"amTableHeader\">\n", file=fileID, append=TRUE)
        if (!is.null(y$metaData))
            headerNames <- c("", "", "", dimnames(y$multilocus)[[2]], "Type")
        else
            headerNames <- c("", "", dimnames(y$multilocus)[[2]], "Type")
        sapply(headerNames, function(j)
                  cat(paste("<td class=\"amPointer\"><table cellspacing=\"0\" class=\"amTableHeaderBtn\"><tr><td>",
                         j,
                         "</td><td style=\"width:10px;\">&nbsp;</td></tr></table></td>\n", sep=""), file=fileID, append=TRUE))
        cat("</tr>\n", file=fileID, append=TRUE)
        
        ## Not matched rows
        for (k in 1:nrow(y$multilocus)) {
            cat("<tr onmouseout=\"this.style.cssText='background-color:none;';\" onmouseover=\"this.style.cssText='background-color:#E0FFFF';\" style=\"\">\n", file=fileID, append=TRUE)
            if (!is.null(y$metaData)) {
                matchRow <- c(paste("<span style=\"color:red; font-weight:bold;\">!!!</span>", sep=""), y$index[k], y$metaData[k], y$multilocus[k, ], "UNCLASSIFIED")
            }
            else {
                matchRow <- c(paste("<span style=\"color:red; font-weight:bold;\">!!!</span>", sep=""), y$index[k], y$multilocus[k, ], "UNCLASSIFIED")
            }
            sapply(matchRow, function(j)
                  cat(paste("<td><div>",
                         j,
                         "</div></td>\n", sep=""), file=fileID, append=TRUE))
            cat("</tr>\n", file=fileID, append=TRUE)
        }
    
        ## End of content tables
        cat("</table></div>", file=fileID, append=TRUE)
        
        ## Content below header names
        cat("<div class=\"amGridLowerPanel\">", file=fileID, append=TRUE)
        cat(paste(c("<span>There were ", nrow(y$multilocus), " samples that were not classified using the parameters supplied.<br></span>
                    <span style=\"color:red;\">Unclassified samples are often not sufficiently similar to match unique genotypes, and not sufficiently different to be declared unique themselves.
                    <br>This is typically because of missing data.  Please see the supplementary documentation for more information.</span>"), collapse=""), file=fileID, append=TRUE)
        cat("</div>", file=fileID, append=TRUE)
        cat("</td></tr></table>", file=fileID, append=TRUE)
        cat("<br><br>", file=fileID, append=TRUE)
    }

    y <- x$pairwise
    ## Do tables for each focal genotype
    for (i in 1:length(y)) {

        ## Table header
        cat("<table cellspacing=\"0\" class=\"amGrid\"><tr><td class=\"amGridContent\">", file=fileID, append=TRUE)
               
        ## Content above header names
        cat("<div class=\"amGridUpperPanel\">", file=fileID, append=TRUE)
        cat(paste(c("<a name=\"", y[[i]]$focal$index, "\"><span style=\"font-size:medium;\">Unique genotype (", i, " of ", length(y), ")</span></a></div>"), collapse=""), file=fileID, append=TRUE)
        
       
        ## Header names
        cat("<div class=\"amGridMiddlePanel\">\n", file=fileID, append=TRUE)
        cat("<table cellspacing=\"0\" class=\"amTable amTableSeparate\">\n", file=fileID, append=TRUE)
        cat("<tr class=\"amTableHeader\">\n", file=fileID, append=TRUE)
        if (!is.null(y[[i]]$focal$metaData))
            headerNames <- c("", "", "", dimnames(y[[i]]$focal$multilocus)[[2]], "Psib", "Score", "Type")
        else
            headerNames <- c("", "", dimnames(y[[i]]$focal$multilocus)[[2]], "Psib", "Score", "Type")
        sapply(headerNames, function(j)
                  cat(paste("<td class=\"amPointer\"><table cellspacing=\"0\" class=\"amTableHeaderBtn\"><tr><td>",
                         j,
                         "</td><td style=\"width:10px;\">&nbsp;</td></tr></table></td>\n", sep=""), file=fileID, append=TRUE))
        cat("</tr>\n", file=fileID, append=TRUE)
        
        ## Content of tables
        cat("<tr>\n", file=fileID, append=TRUE)
    
        ## Focal genotype row
        if (!is.null(y[[i]]$focal$metaData))
            focalRow <- c("", y[[i]]$focal$index, y[[i]]$focal$metaData, y[[i]]$focal$multilocus, "", "", y[[i]]$focal$rowFlag)
        else
            focalRow <- c("", y[[i]]$focal$index, y[[i]]$focal$multilocus, "", "", y[[i]]$focal$rowFlag)
        sapply(focalRow, function(row) cat(paste("<td class=\"amTableSelectRow\"><div>",
                                               row,
                                               "</div></td>\n", sep=""), file=fileID, append=TRUE))
        cat("</tr>\n", file=fileID, append=TRUE)
    
        ## Matching genotype rows
        for (k in 1:nrow(y[[i]]$match$multilocus)) {
            psib <- signif(y[[i]]$match$psib[k],2)
            if (is.na(psib)) {
                psib <- "&nbsp; &nbsp; ---"
            }
            else if (psib < 0.00001) {
                psib <- "<0.00001"
            }
            else {
                psib <- format(psib, scientific=FALSE)
            }
           
            if (!is.null(y[[i]]$focal$metaData)) {
                if (y[[i]]$match$rowFlag[k]=="MULTIPLE_MATCH") {
                    matchRow <- c("<span style=\"color:red; font-weight:bold;\">!!!</span>", y[[i]]$match$index[k], y[[i]]$match$metaData[k], y[[i]]$match$multilocus[k, ], psib,  y[[i]]$match$score[k], "MULTIPLE_MATCH")
                    matchFlags <- c(1, 1, 1, y[[i]]$match$flags[k, ], 1, 1, 1)
                }
                else {
                    matchRow <- c("", y[[i]]$match$index[k], y[[i]]$match$metaData[k], y[[i]]$match$multilocus[k, ], psib,  y[[i]]$match$score[k], "MATCH")
                    matchFlags <- c(1, 1, 1, y[[i]]$match$flags[k, ], 1, 1, 1)
                }
            }
            else {
                if (y[[i]]$match$rowFlag[k]=="MULTIPLE_MATCH") {
                    matchRow <- c("<span style=\"color:red; font-weight:bold;\">!!!</span>", y[[i]]$match$index[k],  y[[i]]$match$multilocus[k, ], psib,  y[[i]]$match$score[k], "MULTIPLE_MATCH")
                    matchFlags <- c(1, 1, y[[i]]$match$flags[k, ], 1, 1, 1)
                }
                else {
                    matchRow <- c("", y[[i]]$match$index[k],  y[[i]]$match$multilocus[k, ], psib,  y[[i]]$match$score[k], "MATCH")
                    matchFlags <- c(1, 1, y[[i]]$match$flags[k, ], 1, 1, 1)
                }
            }
            matchRowHTML <- "<tr>"
            for (m in 1:length(matchRow)) {
                if (matchFlags[m] == 1) {
                    matchRowHTML <- paste(matchRowHTML, "<td><div>", matchRow[m], "</div></td>", sep="")
                }
                else if (matchFlags[m] == 0) {
                    matchRowHTML <- paste(matchRowHTML, "<td><div class=\"amMismatchAllele\">", matchRow[m], "</div></td>", sep="")
                }
                else if(matchFlags[m] == 2) {
                    matchRowHTML <- paste(matchRowHTML, "<td><div class=\"amMissingAllele\">", matchRow[m], "</div></td>", sep="") 
                }
            }
            matchRowHTML <- paste(matchRowHTML, "</tr>\n")
            cat(matchRowHTML, file=fileID, append=TRUE) 
        }

 
        ## End of content tables
        cat("</table></div>", file=fileID, append=TRUE)
        
        ## Content below header names
        cat("<div class=\"amGridLowerPanel\">", file=fileID, append=TRUE)
        
        if ((y[[i]]$match$perfect==0) && (y[[i]]$match$partial==0)) {
            cat("<span>No matches found.</span>", file=fileID, append=TRUE)
        }
        else {
            if (any(y[[i]]$match$rowFlag=="MULTIPLE_MATCH")) {
                multipleMatchMsg <- "<br><span style=\"color:red\">multipleMatch samples, which match two or more unique genotypes, should be manually reviewed.  Please see supplementary documentation for more information</span>"
            }
            else {
                multipleMatchMsg <- ""
            }
            if (x$doPsib!="all") {
                doPsibMsg <- "<br>Psib is calculated for samples that have no mismatches.  Differences among samples are due to loci with missing data."
            }
            else {
                doPsibMsg <- "<br>Psib is calculated for all samples.  Loci with one or both mismatching alleles are treated as missing data."
            }
            
            cat(paste(c("<span>Unique genotype compared against ", x$comparisonDatasetN, " samples, returning those with score>=", round(x$matchThreshold, 3),
                      multipleMatchMsg, doPsibMsg), collapse=""), file=fileID, append=TRUE)
        }
        cat("</div>", file=fileID, append=TRUE)
        cat("</td></tr></table>", file=fileID, append=TRUE)
        cat("<br><br>", file=fileID, append=TRUE)
    }
    
    
    ## Page-wide
    cat("<span style=\"font-size:x-small;\">Generated by allelematch:  an R package<br></span>", file=fileID, append=TRUE)
    cat("<span style=\"font-size:x-small;\">To reference this analysis please use citation(\"allelematch\")<br><br></span>", file=fileID, append=TRUE)
    cat("</div>", file=fileID, append=TRUE)
    cat("</body></html>", file=fileID, append=TRUE)
    close(fileID)
    
    ## If htmlFile is not given, assume that this is being run
    ## for a quick view, therefore open produced html file in browser
    ## Also take this opportunity to clean up temp folder, if required
    if (usingTmpFile) {
        cat("Opening HTML file (", htmlFile, ") in default browser...\n", sep="")
        browseURL(htmlFile)
        if (htmlFilePath != ".") {
            oldTmpFiles <- Sys.glob(paste(htmlFilePath, "\\am*.htm", sep=""))
            if (length(oldTmpFiles) > 0) {
                deleteFiles <- file.remove(oldTmpFiles[sapply(file.info(Sys.glob(paste(htmlFilePath, "\\am*.htm", sep="")))$ctime, function(x) difftime(Sys.time(), x, units="hours")>24)])
                if (sum(deleteFiles) > 0) cat("Deleting allelematch HTML output older than 24 hours from temporary folder...\n")
            }
        }
    }
}



##
## amCSV.amPairwise()
amCSV.amPairwise <- function(x, csvFile) {
    
    if (class(x) != "amPairwise") {
        stop("allelematch:  this function requires an \"amPairwise\" object", call.=FALSE)
    }
     
    y <- x$pairwise
    
    ## Determine the names of columns
    columnNames <- dimnames(y[[1]]$focal$multilocus)[[2]]
    if (!is.null(y[[1]]$match$metaData)) {
        columnNames <- c("comparisonMetaData", columnNames)
    }
    
    ## Add in other columns
    columnNames <- c("matchGroup", "nMatchGroup", "focalIndex", "comparisonIndex", "matchThreshold", "score", columnNames)
    
    ## Create empty matrix to hold data
    csvTable <- matrix(NA, 1, length(columnNames))
    dimnames(csvTable)[[2]] <- columnNames
    for (i in 1:length(y)) {
        
        if (!is.null(y[[1]]$match$metaData)) {
            thisMatch <- cbind(matchGroup=i,
                               nMatchGroup=nrow(y[[i]]$match$multilocus),
                               focalIndex=y[[i]]$focal$index,
                               comparisonIndex=y[[i]]$match$index,
                               matchThreshold=as.character(x$matchThreshold),
                               score=y[[i]]$match$score,
                               comparisonMetaData=y[[i]]$match$metaData,
                               y[[i]]$match$multilocus)
        }
        else {
            thisMatch <- cbind(matchGroup=i,
                               nMatchGroup=nrow(y[[i]]$match$multilocus),
                               focalIndex=y[[i]]$focal$index,
                               comparisonIndex=y[[i]]$match$index,
                               matchThreshold=as.character(x$matchThreshold),
                               score=y[[i]]$match$score,
                               y[[i]]$match$multilocus)
        }
        
        csvTable <- rbind(csvTable, thisMatch)                                   
            
    }
    csvTable <- csvTable[2:nrow(csvTable), ]

    write.csv(csvTable, file=csvFile, row.names=FALSE)
}

##
## amCSV.amCluster()
amCSV.amCluster <- function(x, csvFile) {
     
    if (class(x) != "amCluster") {
        stop("allelematch:  this function requires an \"amCluster\" object", call.=FALSE)
    }
    
    y <- x$unique
    if (!is.null(y$metaData)) {
        csvTable <- cbind(index=y$index, metaData=y$metaData, y$multilocus)
    }
    else {
        csvTable <- cbind(index=y$index, y$multilocus)
    }
    
    write.csv(csvTable, file=csvFile, row.names=FALSE)
}
#
## amCSV.amUnique()
amCSV.amUnique <- function(x, csvFile, uniqueOnly=FALSE) {
        
    if (class(x) != "amUnique") {
        stop("allelematch:  this function requires an \"amUnique\" object", call.=FALSE)
    }

    if (uniqueOnly) {
        ## Pretend it's an amCluster object
        class(x) <- "amCluster"
        amCSV.amCluster(x, csvFile)
    }
    else {
        
        ## Determine the names of columns
        columnNames <- dimnames(x$pairwise[[1]]$focal$multilocus)[[2]]
        if (!is.null(x$pairwise[[1]]$match$metaData)) {
            columnNames <- c("metaData", columnNames)
        }
        
        ## Add in other columns
        columnNames <- c("uniqueGroup", "rowType", "uniqueIndex", "matchIndex", "nUniqueGroup", "alleleMismatch", "matchThreshold", "cutHeight", "Psib", "score", columnNames)
        
        ## Create empty matrix to hold data
        csvTable <- matrix(NA, 1, length(columnNames))
        dimnames(csvTable)[[2]] <- columnNames
        
        ## Unclassifieds...if there are any
        ## If two or more unclassifieds then use cbind()
        if (x$numUnclassified >= 2) {
            y <- x$unclassified
            unclassified <- cbind(uniqueGroup="",
                                  rowType="UNCLASSIFIED",
                                  uniqueIndex=y$index,
                                  matchIndex="",
                                  nUniqueGroup="",
                                  alleleMismatch=as.character(x$alleleMismatch),
                                  matchThreshold=as.character(x$matchThreshold),
                                  cutHeight=as.character(x$cutHeight),
                                  Psib="",
                                  score="")
            if (!is.null(y$metaData)) {
                unclassified <- cbind(unclassified,
                                   metaData=y$metaData,
                                   y$multilocus)
            }
            else {
                unclassified <- cbind(unclassified,
                                   y$multilocus)
            }
            csvTable <- rbind(csvTable, unclassified)
        }
        ## If only one unclassified use c()
        else if (x$numUnclassified == 1) {
            y <- x$unclassified
            unclassified <- c(uniqueGroup="",
                                  rowType="UNCLASSIFIED",
                                  uniqueIndex=y$index,
                                  matchIndex="",
                                  nUniqueGroup="",
                                  alleleMismatch=as.character(x$alleleMismatch),
                                  matchThreshold=as.character(x$matchThreshold),
                                  cutHeight=as.character(x$cutHeight),
                                  Psib="",
                                  score="")
            if (!is.null(y$metaData)) {
                unclassified <- c(unclassified,
                                   metaData=y$metaData,
                                   y$multilocus)
            }
            else {
                unclassified <- c(unclassified,
                                   y$multilocus)
            }
            csvTable <- rbind(csvTable, unclassified)
        }
        
        
        ## Unique genotype groups
        y <- x$pairwise
        for (i in 1:length(y)) {
            
            ## Clean up some variables      
            rowFlag <- y[[i]]$match$rowFlag
            rowFlag[rowFlag==""] <- "MATCH"
            Psib <- y[[i]]$match$psib
            Psib[is.na(Psib)] <- ""
            
            ## Focal row first
            thisMatch <- cbind(uniqueGroup=i,
                                   rowType=y[[i]]$focal$rowFlag,
                                   uniqueIndex=y[[i]]$focal$index,
                                   matchIndex=y[[i]]$focal$index,
                                   nUniqueGroup=nrow(y[[i]]$match$multilocus),
                                   alleleMismatch=as.character(x$alleleMismatch),
                                   matchThreshold=as.character(x$matchThreshold),
                                   cutHeight=as.character(x$cutHeight),
                                   Psib=y[[i]]$match$psib[1],
                                   score="")
            
            if (!is.null(y[[1]]$match$metaData)) {
                thisMatch <- cbind(thisMatch,
                                   metaData=y[[i]]$focal$metaData,
                                   y[[i]]$focal$multilocus)
            }
            else {
                thisMatch <- cbind(thisMatch,
                                   y[[i]]$focal$multilocus)
            }
            csvTable <- rbind(csvTable, thisMatch)

      
            ## Match rows next..if there are any
            ## If two or more matches then use cbind()
            if (sum(y[[i]]$match$index != y[[i]]$focal$index) >= 2) {

                thisMatch <- cbind(uniqueGroup=i,
                                       rowType=rowFlag[y[[i]]$match$index != y[[i]]$focal$index],
                                       uniqueIndex=y[[i]]$focal$index,
                                       matchIndex=y[[i]]$match$index[y[[i]]$match$index != y[[i]]$focal$index],
                                       nUniqueGroup=nrow(y[[i]]$match$multilocus),                                   
                                       alleleMismatch=as.character(x$alleleMismatch),
                                       matchThreshold=as.character(x$matchThreshold),
                                       cutHeight=as.character(x$cutHeight),
                                       Psib=Psib[y[[i]]$match$index != y[[i]]$focal$index],
                                       score=y[[i]]$match$score[y[[i]]$match$index != y[[i]]$focal$index])
                if (!is.null(y[[1]]$match$metaData)) {
                    thisMatch <- cbind(thisMatch,
                                       metaData=y[[i]]$match$metaData[y[[i]]$match$index != y[[i]]$focal$index],
                                       y[[i]]$match$multilocus[y[[i]]$match$index != y[[i]]$focal$index, ])
                }
                else {
                    thisMatch <- cbind(thisMatch,
                                       y[[i]]$match$multilocus[y[[i]]$match$index != y[[i]]$focal$index, ])
                }
                csvTable <- rbind(csvTable, thisMatch)
            }
            ## If only one match use c()
            else if (sum(y[[i]]$match$index != y[[i]]$focal$index) == 1) {
                
                thisMatch <- c(uniqueGroup=i,
                                       rowType=rowFlag[y[[i]]$match$index != y[[i]]$focal$index],
                                       uniqueIndex=y[[i]]$focal$index,
                                       matchIndex=y[[i]]$match$index[y[[i]]$match$index != y[[i]]$focal$index],
                                       nUniqueGroup=nrow(y[[i]]$match$multilocus),                                   
                                       alleleMismatch=as.character(x$alleleMismatch),
                                       matchThreshold=as.character(x$matchThreshold),
                                       cutHeight=as.character(x$cutHeight),
                                       Psib=Psib[y[[i]]$match$index != y[[i]]$focal$index],
                                       score=y[[i]]$match$score[y[[i]]$match$index != y[[i]]$focal$index])
                if (!is.null(y[[1]]$match$metaData)) {
                    thisMatch <- c(thisMatch,
                                       metaData=y[[i]]$match$metaData[y[[i]]$match$index != y[[i]]$focal$index],
                                       y[[i]]$match$multilocus[y[[i]]$match$index != y[[i]]$focal$index, ])
                }
                else {
                    thisMatch <- c(thisMatch,
                                       y[[i]]$match$multilocus[y[[i]]$match$index != y[[i]]$focal$index, ])
                }
                csvTable <- rbind(csvTable, thisMatch)
            }
                
        }
        csvTable <- csvTable[2:nrow(csvTable), ]
    
        write.csv(csvTable, file=csvFile, row.names=FALSE)
    }
}




