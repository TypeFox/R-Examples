### estimateWeights uses an Rcpp module for faster computation of co-occurrence counts.
### See R documentation for more info.

estimateWeights <- function(cuesOutcomes, removeDuplicates=TRUE, saveCounts=FALSE, verbose=FALSE, trueCondProb=TRUE, addBackground=FALSE, hasUnicode=FALSE, ...) {
    ## Internal functions.
    ## Convert UTF8 code points to integers
    toInt = function(char) {
        if (char != '_') { 
            return(paste(utf8ToInt(char),"*",sep='')) 
        } else { 
            return(char)
        }
    }

    ## Convert integers to UTF8 code points
    toUTF8 = function(token) {
        token = sub('_','',token)
        tok = type.convert(token,as.is=T)
            if (is.numeric(tok)) {
                return(intToUtf8(tok)) 
            } else { 
                return(token)
            }
    }

    ## Operates on lists of strings
    convertCues = function(cuelist) { 
        return(paste(unlist(lapply(unlist(strsplit(cuelist,'')), FUN=toInt)),collapse=''))
    }

    ## Operates on lists of strings
    convertBack = function(cuelist) { 
        return(paste(unlist(lapply(unlist(strsplit(cuelist,'*',fixed=T)), FUN=toUTF8)),collapse=''))
    } 

    if (!is.data.frame(cuesOutcomes)) {
        stop("Error: The argument 'cuesOutcomes' must be a dataframe.")
    }

    basename <- NULL
    basename = paste(substitute(cuesOutcomes))
    coocFile = paste(basename,".coocCues.rds",sep='')
    coocOutFile = paste(basename,".coocCuesOutcomes.rds",sep='')
    # Check if pre-computed counts exist. If so, load them
    if (file.exists(coocFile) && file.exists(coocOutFile) ) {
        message(paste(c("NOTE: Loading pre-computed coocurrence matrices.\nIgnoring DataFrame '", basename, "' Provided.\nPlease remove the files ",coocFile," and ",coocOutFile, " if this behavior is not desired.")),sep="")
        if (removeDuplicates && verbose) {
            warning("Did not remove duplicates because there were pre-computed cooccurrence matrices availabe. Please Manually remove the files and run again to make sure that duplicates are removed.")
        }
        flush.console()
        coocCues = readRDS(coocFile)
        coocCuesOutcomes = readRDS(coocOutFile)
    } else {
        ## check for valid column names
        if (!("Cues" %in% colnames(cuesOutcomes))) {
            stop("The 'Cues' column is missing from your dataframe. Please correct the column name and try again. ")
        }
        if (!("Outcomes" %in% colnames(cuesOutcomes))) {
            stop("The 'Outcomes' column is missing from your dataframe. Please correct the column name and try again. ")
        }
        if (!("Frequency" %in% colnames(cuesOutcomes))) {
            warning("The 'Frequency' column is missing from your dataframe. A column of constant frequencies (1) was added.")
            cuesOutcomes$Frequency=1
        }
        NA.cues <- which(is.na(cuesOutcomes$Cues))
        NA.outcomes <- which(is.na(cuesOutcomes$Outcomes))
        if(length(NA.cues)>0)
            stop(paste("NA's in 'Cues': ",length(NA.cues)," cases.",sep=""))
        if(length(NA.outcomes)>0)
            stop(paste("NA's in 'Outcomes': ",length(NA.outcomes)," cases.",sep=""))

        ## Fixing unicode sorting errors.
        if (hasUnicode) {
            cuesOutcomes$Cues = unlist(lapply(cuesOutcomes$Cues,FUN=convertCues))
        }

        ## Call Rcpp function to process all events.
        coocCues = matrix()
        # learnLegacy is not exported (private)
        CuAndCo = learnLegacy(DFin=cuesOutcomes, RemoveDuplicates=removeDuplicates, verbose=verbose)
        coocCues = CuAndCo[[1]]
        coocCuesOutcomes = CuAndCo[[2]]
        ## Recommended for removal by Brian Ripley
        #        rm(CuAndCo)
        # gc()
    }
    # At this point we should have cooc counts for Cue-Cue and Cue-Outcome

    if (verbose) message("Starting to process matrices.")
    ## Check sanity of arguments
    if ((addBackground) & (!trueCondProb)) {
        message("*WARNING: Can't add background rates without true conditional probabilities. \n*ACTION: Proceeding without background rates.")
        addBackground = FALSE
    }
    ## If requested, add background rates for Cue-Cue cooccurrence
    if (addBackground & trueCondProb) {
        ## Check first if Environ has already been computed
        if (sum(coocCues["Environ",]) == 0 & sum(coocCues[,"Environ"] == 0)) {
            cueTotals = diag(coocCues)
#            grandTotal = sum(cueTotals)
            coocCues["Environ",] = cueTotals
            coocCues[,"Environ"] = cueTotals
            coocCues["Environ","Environ"] = sum(cuesOutcomes$Frequency)
        }
    }
    else {
        ## remove rows and columns reserved for background rates
        coocCues=coocCues[!rownames(coocCues) %in% "Environ", !colnames(coocCues) %in% "Environ" ,drop=FALSE]
        coocCuesOutcomes=coocCuesOutcomes[!rownames(coocCuesOutcomes) %in% "Environ",,drop=FALSE]
    }

    ## Check for single cue and outcome.
    if ((nrow(coocCuesOutcomes) <2) & (ncol(coocCuesOutcomes) <2)) {
        stop("Your data had only one unique cue and one unique outcome, making it impossible to estimate the 'weight'. Please make sure that your training data has more than one and cue or more than one outcome.")
    }

    ## Save the cooc matrices for later reuse (after doing Background rates and normalization.
    if (saveCounts) {
        if (verbose) message("Completed Event Counts. Saving Cooc Data for future calculations.")
        flush.console()
        saveRDS(coocCues, file=coocFile)
        if (verbose) message(paste("Saved ",coocFile))
        flush.console()
        saveRDS(coocCuesOutcomes, file=coocOutFile)
        if (verbose) message(paste("Saved ",coocOutFile))
        flush.console()
    }


    # At this point we begin to normalize the counts.

    if (trueCondProb) {
        ##Convert Cue-Outcome counts to Cue-Outcome Probabilities using diagonal
        cueTotals = diag(coocCues) 
        cueTotals[cueTotals == 0] = 1
        condProbsCues = coocCues/cueTotals
        probsOutcomesGivenCues = coocCuesOutcomes/cueTotals
    } else {
        ## use the original algorithm for normalization
        rowsums = rowSums(coocCuesOutcomes)
        rowsums[rowsums == 0] = 1
        condProbsCues = coocCues/rowsums
        probsOutcomesGivenCues = coocCuesOutcomes/rowsums
    }


    if (verbose) message("Starting to calculate pseudoinverse.")
    flush.console()
    n = dim(condProbsCues)[1]
    if (n < 20000) {
        pseudoinverse = ginv(condProbsCues)
        ## Could be faster!!! Do some tests!
        ## pseudoinverse = (t(solve(crossprod(condProbsCues),condProbsCues)))
    } else {
        ## Use an approximation of the pseudoinverse here to make this feasible
        ## average hardware.
        if (verbose) message("Number of cues was too large for standard pseudoinverse. Switching to lower-rank approximation.")
        pseudoinverse = random.pseudoinverse(condProbsCues,verbose=verbose)
    }
    ## Calculate the weights by multiplying the pseudoinver of the c-c
    ## counts by the probabilites of the outcomes given the cues.
    weightMatrix = pseudoinverse %*% probsOutcomesGivenCues
    ## Deal with Unicode issue.
    if (hasUnicode) {
        rownames(weightMatrix) = unlist(lapply(rownames(coocCues),FUN=convertBack))
    } else {
        rownames(weightMatrix) = rownames(coocCues)
    }
    colnames(weightMatrix) = colnames(coocCuesOutcomes)
    if (verbose) message("Completed calculations. Returning weight matrix.")
    flush.console()
    return(weightMatrix)
  }



