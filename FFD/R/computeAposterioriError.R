## Freedom from disease nach Ziller et al.
##
## Functions to analyze the achieved statistical error 
## (alpha-error) for different sampling strategies.
##
## Ian Kopacka
## 2010-07-06


###############################################################################
###############################################################################
###############################################################################
###############################################################################
## Berechne den Alpha-Fehler der gezogenen Stichprobe.
## Eingabeparameter:
##   alphaErrorVector...Vektor (numeric). Alpha-Fehler (= 1-Herdensensitivitaet)
##                      der gezogenen Betriebe.
##   nPopulation........Numeric. Populationsgroesse (= Gesamtanzahl der Betriebe)
##   nDiseased..........Numeric. Anzahl der erkrankten Betriebe lt. 
##                      Designpraevalenz
##   method............."exact" for exact error, "approx" for approximation
##                      recommended for nDiseased > 7
## Ausgabe: Numeric. Alpha-Fehler bei gezogener Stichprobe (exakt).
##
## Version 10: 
computeAposterioriError <- function(alphaErrorVector, nPopulation, 
    nDiseased, method = "default"){
    
    if (nDiseased < 1) return(1)

    ## Determine sample method:
    if (!(method %in% c("exact", "approx", "approximation", "default"))){
        stop("Undefined value for argument 'method'; must be either 'exact' or 'approx'.")
    } 
    if (method == "default"){
        if (nDiseased > 6){
            method <- "approx"
        } else {
            method <- "exact"
        }
    }
    
    ## Error check:
    nSample <- length(alphaErrorVector)  # sample size
    if(nPopulation < nSample) stop("nPopulation must not be smaller than nSample")
    if(nPopulation < nDiseased) stop("nPopulation must not be smaller than nDiseased")
        
    if (method == "exact"){  
        ## EXACT EVALUATION OF A-POSTERIORI ALPHA-ERROR:
        ################################################
              
        ## Range of possible number of diseased in sample:
        nDiseasedSampleVector <- max(1, nSample - nPopulation + nDiseased) : 
            min(nDiseased, nSample)
        ## Vector of probabilites that there are y infected herds in the 
        ## sample for all y in nDiseasedSampleVector:
        probabilitiesDiseasedVector <- choose(nPopulation - nSample,
            nDiseased - c(0, nDiseasedSampleVector)) / choose(nPopulation, 
            nDiseased)
        ## Error check:
		if (any(is.nan(probabilitiesDiseasedVector))){
		probabilitiesDiseasedVector <- exp(lfactorial(nPopulation - nSample) - 
			lfactorial(nDiseased - c(0, nDiseasedSampleVector)) - 
		    lfactorial(nPopulation - nSample - nDiseased + 
			    c(0, nDiseasedSampleVector)) - 
            lfactorial(nPopulation) + lfactorial(nDiseased) +
			lfactorial(nPopulation-nDiseased))
        }		
        if (any(is.nan(probabilitiesDiseasedVector))){
			stop("Binoial Coefficient could not be evaluated.")
		}
        ## Vector of the conditional probabilities of finding 0 test positive
        ## herds in the sample, given that there are y infected herds in the 
        ## sample for all y in nDiseasedSampleVector (costly computation due
        ## to combinatorical aspects):
    
        ## Grouping of alpha-values:
        alphaDf <- data.frame(alpha = as.numeric(as.character(names(table(alphaErrorVector)))), 
            freq = as.vector(table(alphaErrorVector)))
        ## Compute all elements of the form choose(ni, ji)*alphai^ji:
        matPot <- sapply(1:min(nDiseased, nSample), function(x) 
            choose(alphaDf$freq,x)*alphaDf$alpha^x)
        if(dim(alphaDf)[1] == 1) matPot <- t(as.matrix(matPot)) 
    
        ## Evaluate conditional probabilities (will later be done in C-function):
        potList <- NULL    
        
        conditionalProbabilitesVector <- 0*nDiseasedSampleVector
        for (ii in seq(along = nDiseasedSampleVector)){
            if(!is.null(potList)){
                potList <- iterationExpList10(nDiseasedSampleVector[ii],potList) 
            } else {
                potList <- sumDecomp10(nDiseasedSampleVector[ii])
            }
                conditionalProbabilitesVector[ii] <- prodCombN10(matPot, potList)
        }
    
        ## Evaluate alpha error:
        alphaError <- sum(probabilitiesDiseasedVector * 
            c(1,conditionalProbabilitesVector))
    } else {
        ## APPROXIMATE EVALUATION OF A-POSTERIORI ALPHA-ERROR:
        ######################################################
        
        ## Compute mean herd sensitivity:
        herdSensMean <- mean(1-alphaErrorVector)
        alphaError <- computePValue(nPopulation = nPopulation, 
            nSample = nSample, nDiseased = nDiseased, 
            sensitivity = herdSensMean, specificity = 1)    
    }
        
    return(alphaError)    
} 

###################################################################
###################################################################
## Helper functions:

prodCombN10 <- function(matPot,potList){
    ## Initialisiere Summe:
    sumOut <- 0
    prodCombNinner <- function(rv,potList,matPot){
        potVec <- potList[[rv]]
            if(length(potVec) == 1){
                out <- sum(matPot[,potVec])
            } else {    
                tablePot <- as.vector(table(potVec))
                if(length(tablePot) == 1){                    
                    ## C-version:
                    #dyn.load("combNsymmWrapper.dll")
                    outList <- .C("combNsymmWrapper", vector = matPot[,potVec[1]], 
                        nVector = as.integer(length(matPot[,potVec[1]])), 
                        nGroups = as.integer(tablePot), result = 0.0, PACKAGE = "FFD")
                    out <- outList$result                    
                    
                } else {           
                    ### C-version:
                    potVec1 <- potVec
                    tablePot <- table(potVec1)

                    alphaVec <- as.vector(matPot[,as.numeric(names(tablePot))])
                    potVec <- as.vector(tablePot)                    
                    nPot <- length(potVec)
                    nAlpha <- dim(matPot)[1]
                    indVec <- rep(1,nAlpha)                    
                    ## C:
                    #dyn.load("combNWrapper.dll")
                    outList <- .C("combNWrapper", alphaVec = alphaVec,
                        nAlpha = as.integer(nAlpha), nPot = as.integer(nPot),
                        potVec = as.integer(potVec), indVec = as.integer(indVec),
                        result = 0.0, PACKAGE = "FFD")
                    out <- outList$result  
                    
                }
            }            
    }
    outVec <- sapply(seq(along = potList), function(rv) prodCombNinner(rv,potList,matPot))
    return(sum(outVec))    
}

###################################################################

sumDecomp10 <- function(nExp){
    if(nExp > 1){
        expList <- sumDecomp10(nExp-1)
        return(iterationExpList10(nExp, expList))    
    } else {
        return(list(1))
    }      
}

###################################################################

iterationExpList10 <- function(nExp, expList){
    ## For all with a length greater than 1 and where the last element
    ## is smaller than the previous one add 1 to the last element:    
    l1 <- c(nExp,lapply(expList, function(x){
        nVec <- length(x)
        if((nVec > 1) && (x[nVec] < x[nVec-1])){
            x[nVec] <- x[nVec] + 1
            return(x)
        } else {
            return(NULL)
        }        
    }))
    l1 <- l1[!sapply(l1, is.null)]
 
    ## Append 1:
    l2 <- lapply(expList, function(x) c(x,1))
    return(c(l1,l2))
}

###################################################################
###################################################################
