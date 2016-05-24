"drm" <- function(
formula, curveid, pmodels, weights, data = NULL, subset, fct, 
type = c("continuous", "binomial", "Poisson", "quantal", "event"), bcVal = NULL, bcAdd = 0, 
start, na.action = na.fail, robust = "mean", logDose = NULL, 
control = drmc(), lowerl = NULL, upperl = NULL, separate = FALSE,
pshifts = NULL)
{
#    ## Matching 'adjust' argument
#    adjust <- match.arg(adjust)

    ## Matching argument values
    type <- match.arg(type)
    
    ## Loading MASS
#    require(MASS, quietly = TRUE)  # used for boxcox and ginv

    ## Setting na.action option
    options(na.action = deparse(substitute(na.action)))

    ## Setting control parameters
    useD <- control$"useD"
    constrained <- control$"constr"
#    maxDose <- control$"maxDose"
    maxIt <- control$"maxIt"
    optMethod <- control$"method"
    relTol <- control$"relTol"
    warnVal <- control$"warnVal"
#    zeroTol <- control$"zeroTol"
#    bcConstant <- bcAdd  
    rmNA <- control$"rmNA"  # in drmEM...
    errorMessage <- control$"errorm"  # in drmOpt
    noMessage <- control$"noMessage"  # reporting finding control measurements? 
#    trace <- control$"trace"
#    otrace <- control$"otrace"
    dscaleThres <- control$"dscaleThres"
    rscaleThres <- control$"rscaleThres"
        
    ## Setting warnings policy
    options(warn = warnVal)

#    ## Setting adjustment
#    if (adjust == "none") {boxcox <- FALSE; varPower <- FALSE}
#    if (adjust == "bc1") {boxcox <- TRUE; varPower <- FALSE}
#    if (adjust == "vp") {boxcox <- FALSE; varPower <- TRUE}

#    if ( (!is.null(bcVal)) && (is.numeric(bcVal))) {boxcox <- bc}
        
#    if (!(robust == "mean"))
#    {
#        boxcox <- FALSE
#        varPower <- FALSE
#    }

    ## Handling 'start' argument
    if (missing(start)) {selfStart <- TRUE} else {selfStart <- FALSE}

    ## Handling 'fct' argument
    if ( (!is.list(fct)) && (!is.function(fct)) ) {stop("No function or list given in argument 'fct'")}
    if (is.function(fct)) 
    {
        fct <- fct2list(fct, 2)
    }
    
    ## Converting a user specified list
    if (is.null(names(fct))) {fct$"fct" <- fct[[1]]; fct$"ssfct" <- fct[[2]]; fct$"names" <- fct[[3]]}
    
    if (!is.function(fct$"fct")) 
    {
        stop("First entry in list to 'fct' NOT a function")
    } else {
        drcFct <- fct$"fct"
    }
    
    if (is.null(fct$"ssfct")) {noSSfct <- TRUE} else {noSSfct <- FALSE}
    if ((!is.function(fct$"ssfct")) && selfStart)
    {
        stop("Neither self starter function nor starting values provided")
    } else {
        ssfct <- fct$"ssfct"
    }
    
    if (is.null(fct$"names") || (!is.character(fct$"names"))) 
    {
        stop("Parameter names (as vector a strings) are NOT supplied")
    } else {
        parNames <- fct$"names" 
        numNames <- length(parNames)
    }
    
#    ## Coercing two arguments in 'ssfct' into one argument
#    lenASS <- length(formals(ssfct))
#    if (lenASS > 1)
#    {
#        stop("Self starter function should only have one argument, which takes a data frame")
##        ssTemp <- ssfct
##        ssfct <- function(dataset) {ssTemp(dataset[, head(1:lenASS, -1)], dataset[, lenASS])}
#    }

    ## Checking whether or not first derivates are supplied    
    isDF <- is.function(fct$"deriv1")
    if ( (useD) && (isDF) )
    {
        dfct1 <- fct$"deriv1"  # deriv1  # [[4]]
#        drcDer2 <- fct$deriv2  # [[5]]
    } else {
        dfct1 <- NULL
    }

    ## Checking whether or not second derivates are supplied    
    if ( (useD) && (is.function(fct$"deriv2")) )
    {
        dfct2 <- fct$"deriv2"
    } else {
        dfct2 <- NULL
    }    
#    fct$"anovaYes"$"bin" <- NULL
#    fct$"anovaYes"$"cont" <- TRUE

    ## Storing call details
    callDetail <- match.call()

    ## Handling the 'formula', 'curveid' and 'data' arguments
    anName <- deparse(substitute(curveid))  # storing name for later use
    if (length(anName) > 1) {anName <- anName[1]}  # to circumvent the behaviour of 'substitute' in do.call("multdrc", ...)
    if (nchar(anName) < 1) {anName <- "1"}  # in case only one curve is analysed


    mf <- match.call(expand.dots = FALSE)   
    nmf <- names(mf) 
    mnmf <- match(c("formula", "curveid", "data", "subset", "na.action", "weights"), nmf, 0) 

    mf[[1]] <- as.name("model.frame")
    mf <- eval(mf[c(1,mnmf)], parent.frame())  #, globalenv())
    mt <- attr(mf, "terms")
        
    dose <- model.matrix(mt, mf)[,-c(1)]  # with no intercept
    resp <- model.response(mf, "numeric")
    origDose <- dose
    origResp <- resp  # in case of transformation of the response    
    lenData <- length(resp)
    numObs <- length(resp)

    xDim <- ncol(as.matrix(dose))
#    if (xDim > 1)
#    {
#        stop("drm() is only designed for 1-dim. dose vectors")
#    }

#    dimData <- xDim + 1  # dimension of dose plus 1 dimensional response
    
#    varNames <- names(mf)
#    varNames <- varNames[c(2:dimData,1)]

#    print(names(mf))
#    print(model.extract(mf, "weights"))
#    print(model.weights(mf))
    varNames <- names(mf)[c(2, 1)]  
    varNames0 <- names(mf) 
    # only used once, but mf is overwritten later on

    ## Retrieving weights
    wVec <- model.weights(mf)
    if (is.null(wVec))
    {
        wVec <- rep(1, numObs)
    }

#    ## Extracting variable for heterogeneous variances
#    vvar <- model.extract(mf, "hetvar")    
    
    ## Finding indices for missing values
    missingIndices <- attr(mf, "na.action")
    if (is.null(missingIndices)) {removeMI <- function(x){x}} else {removeMI <- function(x){x[-missingIndices,]}}

    ## Handling "curveid" argument
    assayNo <- model.extract(mf, "curveid") 
    if (is.null(assayNo))  # in case not supplied
    {
        assayNo <- rep(1, numObs)
    }
    uniqueNames <- unique(assayNo)
    colOrder <- order(uniqueNames)
#    print(colOrder)
    uniqueNames <- as.character(uniqueNames)

    ## Re-enumerating the levels in 'assayNo' and 'pmodels'
    assayNoOld <- assayNo
#    ciOrigIndex <- uniqueNames  # unique(assayNo)
#    ciOrigLength <- length(unique(assayNoOld))

    ## Detecting control measurements 
    
    ## Defining helper function     
    colConvert <- function(vec)
    {
        len <- length(vec)
        assLev <- unique(vec)

        retVec <- rep(0,len)
        j <- 1
        for (i in 1:length(assLev)) {retVec[vec == assLev[i]] <- j; j <- j + 1}

        return(retVec)
    } 
    assayNo <- colConvert(assayNoOld)
    assayNames <- as.character(unique(assayNoOld))
    numAss <- length(assayNames)
    
    
#    lenDose <- unlist(lapply(tapply(dose, assayNoOld, unique), length))
#    conDose <- names(lenDose)[lenDose == 1]
#    nconDose <- names(lenDose)[lenDose > 1]
#    if (length(conDose) > 0) 
#    {        
#        if (!noMessage) 
#        {
#            cat(paste("Control measurements detected for level: ", conDose, "\n", sep = ""))
#        }
#            
#        assayNo[assayNoOld %in% conDose] <- nconDose[1]
#        ciOrigIndex <- unique(assayNo)
##        ciOrigLength <- length(unique(assayNoOld))  # numAss
#      
#      
#        ## Updating names, number of curves and the enumeration (starting from 1)
#        assayNames <- nconDose
##        numAss <- length(assayNames) 
#        assayNo <- colConvert(assayNo)
#                       
#        cm <- NULL
##        }
#         
#        uniqueDose <- lapply(tapply(dose, assayNoOld, unique), length)
#        udNames <- names(uniqueDose[uniqueDose == 1])
#        if (length(udNames) > 0) 
#        {
#            cm <- udNames
#            if (!noMessage) {cat(paste("Control measurements detected for level: ", udNames, "\n", sep = ""))}
#            ## add a check to see if at least one component in pmodels results in a single column
            
##            conInd <- assayNoOld%in%udNames
##            assayNo[conInd] <- (assayNo[!conInd])[1]
##            cm <- NULL
##assayNew <- assayNo
##assayNew[conInd] <- (assayNo[!conInd])[1]
##print(assayNew)
##
#            conInd <- assayNoOld%in%udNames
#            assayNo[conInd] <- (assayNo[!conInd])[1]
#            ciOrigIndex <- unique(assayNo)
#            ciOrigLength <- numAss
#            
#            ## Updating names, number of curves and the enumeration (starting from 1)
#            assayNames <- as.character(unique(assayNoOld[!conInd]))
#            numAss <- length(assayNames) 
#            assayNo <- colConvert(assayNo)
#                       
#            cm <- NULL
#

## New -commented out
#    } else {
#        cm <- NULL
#        ciOrigIndex <- unique(assayNo)
##        ciOrigLength <- numAss 
#        
#        assayNames <- as.character(unique(assayNoOld))          
#        assayNo <- colConvert(assayNoOld)  # re-enumerating from 1 to numAss       
#    }
#    numAss <- length(assayNames)   
#    print(ciOrigIndex)
#    print(ciOrigLength)


    if (xDim > 1) {tempDoseVec <- dose[, 1]} else {tempDoseVec <- dose} 
#    uniqueDose <- lapply(tapply(dose, assayNoOld, unique), length)
    uniqueDose <- lapply(tapply(tempDoseVec, assayNoOld, unique), length)
    udNames <- names(uniqueDose[uniqueDose == 1])
    if (length(udNames) > 0) 
    {
        cm <- udNames
        if (!noMessage) 
        {
            cat(paste("Control measurements detected for level: ", udNames, "\n", sep = ""))
            
            if (separate)
            {
                stop("Having a common control when fitting separate models does not make sense!\n")
            }
        }
        conInd <- assayNoOld %in% udNames
        assayNo[conInd] <- (assayNo[!conInd])[1]
        ciOrigIndex <- unique(assayNo)
        ciOrigLength <- numAss
        
        ## Updating names, number of curves and the enumeration (starting from 1)
        assayNames <- as.character(unique(assayNoOld[!conInd]))
        numAss <- length(assayNames)
        assayNo <- colConvert(assayNo)
        cm <- NULL
    } else {
        cm <- NULL
        ciOrigIndex <- unique(assayNo)
        ciOrigLength <- numAss
    }
#    print(assayNo)
    
    ## Pooling data from different curves
    if ((separate) && (numAss < 2))
    {
#        warning("Nothing to pool", call. = FALSE)
        warning("Only one level: separate = TRUE has no effect", call. = FALSE)
        separate <- FALSE 
    }    
    if ((separate) && (!missing(pmodels)))
    {
        warning("Separate fitting switched off", call. = FALSE)
        separate <- FALSE
    }
    if (separate)
    {
#        return(idrm(dose, resp, assayNo, wVec, fct, type))
        return(idrm(dose, resp, assayNoOld, wVec, fct, type, control))
    }    
    
    ## Handling "pmodels" argument
    pmodelsList <- list()
    if (missing(pmodels)) 
    {
#        pmodels <- as.data.frame(matrix(assayNo, numObs, numNames))
#        
        if (length(unique(assayNo)) == 1) 
        {
            for (i in 1:numNames) 
            {
                pmodelsList[[i]] <- matrix(1, numObs, 1)
            }
        } else {
            modelMat <- model.matrix(~ factor(assayNo) - 1, level = unique(assayNo))  # no intercept term
            colnames(modelMat) <- assayNames
            for (i in 1:numNames) 
            {
                pmodelsList[[i]] <- modelMat
#                print(head(modelMat))
#                pmodelsList[[i]] <- pmodelsList[[i]][, colOrder]
            }
        }   
    } else {
        ## Handling a list or data.frame argument of "pmodels"
        if (is.null(data)) 
        {
            pmodels <- eval(substitute(pmodels), envir = .GlobalEnv)
        } else {
            pmodels <- eval(substitute(pmodels), envir = data, enclos = parent.frame())
        }
   
        if (is.data.frame(pmodels))
        {
            lenCol <- ncol(pmodels)
            pmodelsMat <- matrix(0, numObs, lenCol)    
    
            for (i in 1:lenCol) 
            {
                if (length(unique(pmodels[,i])) == 1) 
                {
                    pmodelsList[[i]] <- matrix(1, numObs, 1)
                    pmodelsMat[,i] <- rep(1, numObs)    
                } else {
                    mf <- eval(model.frame(~factor(pmodels[,i]) - 1), parent.frame())  # converting to factors
                    mt <- attr(mf, "terms")    
    
                    mf2 <- model.matrix(mt, mf)
                    ncmf2 <- ncol(mf2)

                    mf3 <- removeMI(mf2)
                    pmodelsList[[i]] <- mf3
                    pmodelsMat[, i] <- mf3 %*% c(1:ncmf2)
                }
            }
        } else {

            if (is.list(pmodels))
            {   
                lenCol <- length(pmodels)
                pmodelsMat <- matrix(0, length(resp), lenCol)
    
                for (i in 1:lenCol) 
                {
                    if (paste(as.character(pmodels[[i]]), collapse = "") == "~1") 
                    {
                        pmodelsList[[i]] <- matrix(1, numObs, 1)
                        pmodelsMat[,i] <- rep(1, numObs)
                    } else {
                        mf <- eval(model.frame(pmodels[[i]], data=data), parent.frame())   
                        mt <- attr(mf, "terms")    
                        
                        mf2 <- model.matrix(mt, mf)
                        ncmf2 <- ncol(mf2)

                        mf3 <- removeMI(mf2)                    
                        pmodelsList[[i]] <- mf3  
                    
                        pmodelsMat[,i] <- mf3%*%c(1:ncmf2)                    
                    }
                }
            }
        }     
#        pmodelsOld <- pmodels
#        pmodels <- as.data.frame(pmodelsMat)  # pmodelsMat not used any more
    }
#    for (i in 1:numNames) {pmodels[, i] <- colConvert(pmodels[, i])}

    
    ## Re-setting na.action
    options(na.action = "na.omit")  # the default

    ## Transforming dose value if they are provided as log dose
    if ( !is.null(logDose) && is.numeric(logDose) ) 
    {
       origDose <- dose
       dose <- logDose^dose
    }

#    ## Handling one-dimensional x     
#    if (xDim == 1)
#    {
        ## Defining ANOVA model
#        bcc <- rep(bcAdd, numObs)    
#        if (numAss > 1) 
#        {
#            anovaFormula <- (resp + bcc) ~ offset(bcc) + doseFactor*factor(assayNo)
#            alternative <- 2
#        } else {
#            anovaFormula <- (resp + bcc) ~ offset(bcc) + doseFactor
#            alternative <- 1
#        }
         

#        ## Checking whether there is enough df to perform Box-Cox transformation
#        if ( boxcox && ( (numObs - numAss*length(unique(dose))) < numObs/10) )
#        {
#            if (boxcox) {warning("Box-Cox transformation based on clustering of dose values", call. = FALSE)}
#            doseFactor <- factor(cutree(hclust(dist(dose), method = "average"), numObs/3))  
#            # constructing groups containing roughly 3 observations 
#
#            ## Re-defining ANOVA model
#            if (numAss > 1) 
#            {
#                anovaFormula <- (resp + bcc) ~ offset(bcc) + doseFactor*factor(assayNo)
#                dset <- data.frame(doseFactor, resp, assayNo, bcc)
#                alternative <- 2
#            } else {
#                anovaFormula <- (resp + bcc) ~ offset(bcc) + doseFactor
#                dset <- data.frame(doseFactor, resp, bcc)
#                alternative <- 1
#            }
#        } else {
#            doseFactor <- factor(dose)
#        }
         
#        ## Fitting ANOVA model
#        if (type == "continuous")
#        {
#            testList <- drmLOFls()
##            if (varPower) {testList <- drmLOFvp()}
##            if (!is.null(vvar)) {testList <- drmLOFhv()}
#        }
#        if (type == "binomial")
#        {
#            testList <- drmLOFbinomial()
#        }
#        if (type == "Poisson")
#        {
#            testList <- drmLOFPoisson()
#        }
#        
##        if (varPower) {testList <- mdrcVp(anovaYes = TRUE)} else {testList <- drmEMls(anovaYes = TRUE)}
##        if (!is.null(vvar)) {testList <- mdrcHetVar(anovaYes = TRUE)}
#
#        gofTest <- testList$"gofTest"            
#        lofTest <- testList$"anovaTest"
#        if (!is.null(lofTest))
#        {
#            dset <- data.frame(dose, factor(dose), resp, assayNo, bcc)
#            anovaModel0 <- lofTest(anovaFormula, dset)
#        } else {
#            anovaModel0 <- NULL
#            alternative <- 0
#        }
        
#        ## Fitting ANOVA model
#        testList <- switch(type, 
#        "continuous" = drmLOFls(), 
#        "binomial" = drmLOFbinomial(), 
#        "Poisson" = drmLOFPoisson())
#
#        gofTest <- testList$"gofTest"            
#        lofTest <- testList$"anovaTest"
#        if (!is.null(lofTest))
#        {
#            afList <- anovaFormula(dose, resp, assayNo, bcAdd)
#            anovaForm <- afList$"anovaFormula"
#            anovaData <- afList$"anovaData"        
#
#            anovaModel0 <- lofTest(anovaForm, anovaData)
#        } else {
#            anovaModel0 <- NULL
#        }


#        ## Applying the Box-Cox transformation (lambda is defined here!)
#        bcResult <- drmBoxcox(boxcox, anovaFormula, dset)  
#        lambda <- bcResult[[1]]
#        boxcoxci <- bcResult[[2]] 
#        boxcox <- bcResult[[3]]      

#        lambda <- 0
#        isNumeric <- is.numeric(boxcox)
#        if ( (isNumeric) || (is.logical(boxcox) && boxcox)  ) 
#        {
#            if (!isNumeric)
#            {
#                profLik <- boxcox(anovaFormula, lambda = seq(-2.6, 2.6, 1/10), plotit = FALSE, data = dset)  
#                # boxcox in MASS
#                
#                maxIndex <- which.max(profLik$y)
#                lambda <- (profLik$x)[maxIndex]
#                boxcoxci <- drmBoxcoxCI(profLik)
#            }
#            if (isNumeric)
#            {
#                lambda <- boxcox
#                boxcoxci <- c(NA, NA)                
#            }
#        } else {
#            lambda <- NA
#            boxcoxci <- c(NA, NA)
#        }
 

#        ## Using self starter 
#        if (!noSSfct)
#        {
#            ## Calculating initial estimates for the parameters using the self starter
#            startMat <- matrix(0, numAss, numNames)
#            doseresp <- data.frame(dose, origResp)
#    
#            for (i in 1:numAss)
#            {
#                indexT1 <- (assayNo == i)
#                if (any(indexT1)) 
#                {
#                    isfi <- is.finite(dose)  # removing infinite dose values
#                    logVec <- indexT1 & isfi
#                    startMat[i, ] <- ssfct(doseresp[logVec, ])  # ssfct(dose[logVec], origResp[logVec] )
#                } else {
#                    startMat[i, ] <- rep(NA, numNames)
#                }
#
#                ## Identifying a dose response curve only consisting of control measurements
#                if (sum(!is.na(startMat[i, ])) == 1) {upperPos <- (1:numNames)[!is.na(startMat[i, ])]}
#            }
##            colMat <- matrix(0, numNames, numAss)
##            maxParm <- rep(0, numNames)  # storing the max number of parameters
#        }

        
#    ## Handling multi-dimensional x   
#    } else { 
#        stop("Currently multi-dimensional dose values are not supported") 
#        alternative <- NULL
#        anovaModel0 <- NULL
#        anovaModel <- NULL
#        gofTest <- NULL
        
#        if (!is.null(bcVal))
#        {
#            lambda <- boxcox                  
#            boxcoxci <- c(NA, NA)
#        } else {
#            lambda <- NA                  
#            boxcoxci <- NULL
#        }       

#        ## Using self starter
#        if (!noSSfct)
#        {
#            ## Calculating initial estimates for the parameters using the self starter
#            startMat <- matrix(0, numAss, numNames)
#            doseresp <- data.frame(dose, origResp)
#    
#            for (i in 1:numAss)
#            {
#                indexT1 <- (assayNo == i)
#                if (any(indexT1)) 
#                {
#                    startMat[i, ] <- ssfct(doseresp[indexT1, ])  # ssfct(dose[indexT1], origResp[indexT1])
#                } else {
#                    startMat[i, ] <- rep(NA, numNames)
#                }
#                
#                ## Identifying a dose response curve only consisting of control measurements
#                if (sum(!is.na(startMat[i,]))==1) {upperPos <- (1:numNames)[!is.na(startMat[i,])]}
#            }
#            colMat <- matrix(0, numNames, numAss)
#            maxParm <- rep(0, numNames)  # storing the max number of parameters
#        }
#    }

    ## Finding parameters for the control measurements which will not be estimated
    pmodelsList2 <- list()
    for (i in 1:numNames)
    { 
        colNames <- colnames(pmodelsList[[i]])

        if ( (!is.null(cm)) && (!is.null(colNames)) ) 
        {
            accm <- as.character(cm)
            pos <- grep(accm, colNames)
            if (length(pos) == 0) 
            {
                candCol <- pmodelsList[[i]][, 1]
                if ( !(length(assayNoOld[candCol==1])==0) && (all(assayNoOld[candCol==1] == accm)) )
                {
                    pos <- 1  # the control measurements correspond to the "Intercept" term
                }
            }  
        } else {pos <- numeric(0)}


        ## Defining 'pmodelsList2' from 'pmodelsList'
        if ((length(pos) > 0) && !(upperPos == i) )
        {
            pmodelsList2[[i]] <- as.matrix(pmodelsList[[i]][, -pos])  # column is removed
        } else {
            pmodelsList2[[i]] <- as.matrix(pmodelsList[[i]])  # column is kept
        } 
    }
    
    for (i in 1:numNames)
    {
        if (ncol(pmodelsList[[i]]) > numAss) 
        {
            pmodelsList2[[i]] <- model.matrix(~factor(assayNo) - 1)
            colnames(pmodelsList2[[i]]) <- assayNames 
        } else {
            pmodelsList2[[i]] <- as.matrix(pmodelsList[[i]])  # columns are kept
        }
    }
    
    ## Constructing vectors 'ncclVec' and 'parmPos' used below
    ncclVec <- rep(0, numNames)
    for (i in 1:numNames)
    {
        ncclVec[i] <- ncol(pmodelsList2[[i]])  # ncol(as.matrix(pmodelsList2[[i]]))
    }
    parmPos <- c(0, cumsum(ncclVec)[-numNames])

    ## Constructing parameter names                       
    pnList <- drmParNames(numNames, parNames, pmodelsList2)
    parmVec <- pnList[[1]]
    parmVecA <- pnList[[2]]
    parmVecB <- pnList[[3]]   

    ## Defining with indices for the individual parameters in the model
    parmIndex <- list()
    for (i in 1:numNames)
    {
        parmIndex[[i]] <- parmPos[i] + 1:ncclVec[i]
    }

    ## Scaling of dose and response values 
    scaleFct <- fct$"scaleFct"
    if (!is.null(scaleFct))  # && (is.null(lowerl)) && (is.null(upperl)) )
    # currently the scaling interferes with constraining optimization
    {
        # Defining scaling for dose and response values 
        doseScaling <- 10^(floor(log10(median(dose))))   
#        if ( (is.na(doseScaling)) || (doseScaling < 1e-10) )  # changed May 16 2012
        if ( (is.na(doseScaling)) || (doseScaling < dscaleThres) )
        {
            doseScaling <- 1
        }

        respScaling <- 10^(floor(log10(median(resp)))) 
#        if ( (is.na(respScaling)) || (respScaling < 1e-10) || (!identical(type, "continuous")) || (!is.null(bcVal)) )  # changed May 16 2012
        if ( (is.na(respScaling)) || (respScaling < rscaleThres) || (!identical(type, "continuous")) || (!is.null(bcVal)) )
        {
            respScaling <- 1
        }   
#        print(resp)
#        print(median(resp))

#        doseScaling <- 1
#        respScaling <- 1

## Starting values need to be calculated after BC transformation!!!        

        # Retrieving scaling vector
        longScaleVec <- rep(scaleFct(doseScaling, respScaling), as.vector(unlist(lapply(parmIndex, length))))
        
    } else {
        doseScaling <- 1
        respScaling <- 1
        
        longScaleVec <- 1
#        
#        startVecSc <- startVec
    }      
#    print(c(doseScaling, respScaling, longScaleVec))    

    ## Constructing vector of initial parameter values
    startVecList <- list()
    
    ## Calculating initial estimates for the parameters using the self starter
    if(!noSSfct)
    {
        startMat <- matrix(0, numAss, numNames)
        lenASS <- length(formals(ssfct))
        if (lenASS > 1)  
        # in case doseScaling and respScaling arguments are available
        # scaling is done inside ssfct()
        {
            doseresp <- data.frame(x = dose, y = origResp)
            ssFct <- function(dframe){ssfct(dframe, doseScaling, respScaling)}
        } else {
        # scaling is explicitly applied to the dose and response values
            doseresp <- data.frame(x = dose / doseScaling, y = origResp / respScaling)
            ssFct <- ssfct
        }
#        doseresp <- data.frame(x = dose / doseScaling, y = origResp / respScaling)
#        doseresp <- data.frame(dose, origResp)

# Not sure this indicator is needed?! Only used once below!
# Note is.finite() only works with vectors!
# Commented out 2010-12-13   
        isfi <- is.finite(dose)  # removing infinite dose values
    
        if (identical(type, "event"))
        {
            dr2 <- doseresp[, 3]
#            print(doseresp[, 2:3])
            isFinite <- is.finite(doseresp[, 2])
            respVec <- rep(NA, length(dr2))
            respVec[isFinite] <- cumsum(dr2[isFinite]) / sum(dr2)
#            doseresp[, 3] <- cumsum(dr2[isFinite]) / sum(dr2)
##            doseresp[!is.finite(doseresp[, 2]), 1] <- NA              
#            doseresp <- doseresp[isFinite, c(1, 3)]
#            names(doseresp) <- c("x", "y") 
            doseresp <- (data.frame(x = doseresp[, 1], y = respVec))[isFinite, ]
#            print(doseresp)           
        } else {
            isFinite <- is.finite(doseresp[, 2])
        }
    
        ## Finding starting values for each curve
        for (i in 1:numAss)
        {
            indexT1 <- (assayNo[isFinite] == i)
            if (any(indexT1)) 
            {
# Commented out 2010-12-13            
#                logVec <- indexT1 & isfi
                logVec <- indexT1
                                
#                startMat[i, ] <- ssfct(doseresp[logVec, ])  # ssfct(dose[logVec], origResp[logVec] )
#                startMat[i, ] <- ssfct(doseresp[logVec, ], doseScaling, respScaling)
                startMat[i, ] <- ssFct(doseresp[logVec, ])
            } else {
                 startMat[i, ] <- rep(NA, numNames)
            }
    
            ## Identifying a dose response curve only consisting of control measurements
            if (sum(!is.na(startMat[i, ])) == 1) 
            {
                upperPos <- (1:numNames)[!is.na(startMat[i, ])]
#                print(upperPos)
            }
        }
#        print(startMat)
          
#        startMat2 <- matrix(unlist(lapply(split(doseresp, assayNo[isFinite]), ssFct)), nrow = numAss, byrow = TRUE)
#        upperPos2 <- c(rep(1:numNames, numAss))[t(is.na(startMat))]
#        print(upperPos2)
#        print(startMat2)
  

#       New approach?
#        timeUsed <- 0 
#        ssFctWrapper <- function(dframeSubset)
#        {
#            ssFct(dframeSubset[is.finite(dframeSubset[1, ]), ])
#            timeUsed <- timeUsed + system.time(ssFct(dframeSubset[is.finite(dframeSubset[1, ]), ]))[3]
#        }  
#        startMat2 <- matrix(as.vector(unlist(lapply(split(doseresp, assayNo), ssFctWrapper))), 
#        numAss, numNames, byrow = TRUE)
#        print(startMat2)  # for comparison  
  
    
        ## Transforming matrix of starting values into a vector
        nrsm <- nrow(startMat)
        for (i in 1:numNames)
        {
            sv <- rep(0, max(nrsm, ncclVec[i]))
            indVec <- 1:ncclVec[i]
            sv[1:nrsm] <- startMat[, i]
            sv <- sv[!is.na(sv)]
            
            isZero <- (sv == 0)
            sv[isZero] <- mean(sv)
            
            startVecList[[i]] <- sv[indVec]
#            print(startVecList[[i]])
        }
        startVec <- unlist(startVecList)        
    } else {
        startVec <- start  # no checking if no self starter function is provided!!!
    }
    
    ## Checking the number of start values provided
    if (!selfStart && !noSSfct) 
    {
        lenReq <- length(startVec)  # generated from self starter
        if (length(start) == lenReq) 
        {
            startVec <- start / longScaleVec
        } else {
            stop(paste("Wrong number of initial parameter values. ", lenReq, " values should be supplied", sep = ""))
        }
    }

    ## Converting parameters
    if (selfStart)
    {
        startVec <- drmConvertParm(startVec, startMat, assayNo, pmodelsList2) 
    }

    # Scaling starting values (currently not done in drmEMls)
#    startVecSc <- startVec / longScaleVec 
    startVecSc <- startVec
#    print(startVecSc)    

    ## Defining function which converts parameter vector to parameter matrix            
    parmMatrix <- matrix(0, numObs, numNames)
    parm2mat <- function(parm)
    {
#        parmMatrix <- matrix(0, lenData, numNames)
        for (i in 1:numNames)
        {
#           print(as.matrix(pmodelsList2[[i]]))
#           print(parmPos[i] + 1:ncclVec[i])
#           print(parm[parmPos[i] + 1:ncclVec[i]])          
#           parmMatrix[, i] <- pmodelsList2[[i]] %*% parm[parmPos[i] + 1:ncclVec[i]]
           parmMatrix[, i] <- pmodelsList2[[i]] %*% parm[parmIndex[[i]]]
        }
        return(parmMatrix)
    }        

    ## Defining non-linear function
#    if (!is.null(fctList))
#    {
#        ivList <- list()
#        ivList2 <- list()        
#        matList <- list()
#        svList <- list()
#        for (i in 1:numAss)
#        {
#            indexT1 <- (assayNo == i)
#            isfi <- is.finite(dose)  # removing infinite dose values
#
#            ivList[[i]] <- indexT1
##            svList[[i]] <- fctList[[i]]$"ssfct"( doseresp[(indexT1 & isfi), ] )
#            logVec <- indexT1 & isfi
#            svList[[i]] <- fctList[[i]]$"ssfct"(doseresp[logVec, ])  # dose[logVec], origResp[logVec])
#            matList[[i]] <- c( sum(indexT1), length(svList[[i]]) )
#            
#            ivList2[[i]] <- match(fctList[[i]]$names, fct$names)             
#        }
#
#    
#        posVec <- rep(0, numAss)
#        for (i in 1:numAss)
#        {
#            posVec[i] <- matList[[i]][2]
#        }
#        posVec <- cumsum(posVec)
#        posVec <- c(0, posVec)
##        print(posVec)
#    
#        drcFct1 <- function(dose, parm)
#        {
#            retVec <- rep(0, numObs)
#            for (i in 1:numAss)
#            {
#                iVec <- ivList[[i]]
#                pMat <- matrix(parm[(posVec[i]+1):posVec[i+1]], matList[[i]][1], matList[[i]][2], byrow = TRUE) 
#                retVec[iVec] <- fctList[[i]]$"fct"( dose[iVec], pMat )
#            }
#            return(retVec)
#        }
#        
#        startVec <- as.vector(unlist(svList))
#    } else {

    ## Defining model function 
    multCurves <- modelFunction(dose, parm2mat, drcFct, cm, assayNoOld, upperPos, fct$"retFct", 
                                doseScaling, respScaling, isFinite = rep(TRUE, lenData), pshifts)

#    drcFct1 <- function(dose, parm)
#    {
#        drcFct(dose, parm2mat(parm))
#    }
##    }
#
#
#    ## Defining model function        
#    if (!is.null(fct$"retFct"))
#    {
#        drcFct <- fct$"retFct"(doseScaling, respScaling)  #, numObs)
#        drcFct1 <- function(dose, parm)
#        {
#            drcFct(dose, parm2mat(parm))
#        }
#    }
#
#    if (is.null(cm)) 
#    {
#        multCurves <- function(dose, parm)
#        {          
#           drcFct1(dose, parm)  # fctList           
#        }
#    } else {
#        iv <- assayNoOld == cm
#        niv <- !iv
#        fctEval <- rep(0, numObs)
#
#        multCurves <- function(dose, parm) 
#        {
#            parmVal <- parm2mat(parm)           
#            fctEval[iv] <- parmVal[iv, upperPos, drop = FALSE]
#            fctEval[niv] <- drcFct(dose[niv], parmVal[niv, , drop = FALSE])
#
#            fctEval
#        }
#    }
##    print(startVec)
##    print(multCurves(dose, startVec))

    
    ## Defining first derivative (if available) ... used once in drmEMls()
    if (!is.null(dfct1))
    {
        dmatfct <- function(dose, parm)
        {
            dfct1(dose, parm2mat(parm))
        }
    } else {
        dmatfct <-NULL
    } 

    ## Box-Cox transformation is applied
    if (!is.null(bcVal))  # (boxcox)
    {
#        varPower <- FALSE  # not both boxcox and varPower at the same time

        ## Defining Box-Cox transformation function
        bcfct <- function(x, lambda, bctol, add = bcAdd)
        {
            if (abs(lambda) > bctol)
            {
                return(((x + add)^lambda - 1)/lambda)
            } else {
                return(log(x + add))    
            }
        }
        
        ## Setting the tolerance for Box-Cox transformation being the logarithm transformation 
        ##  (same as in boxcox.default in MASS package)
        bcTol <- 0.02 
        
#        resp <- bcfct(resp, lambda, bcTol)
        resp <- bcfct(resp, bcVal, bcTol)

        multCurves2 <- function(dose, parm)
        {
            bcfct(multCurves(dose, parm), bcVal, bcTol)
        } 
    } else {multCurves2 <- multCurves}
#    print(startVec)
#    print(multCurves2(dose, startVec))    

    ## Defining estimation method -- perhaps working for continuous data
#    robustFct <- drmRobust(robust, match.call(), numObs, length(startVec))  
    robustFct <- drmRobust(robust, callDetail, numObs, length(startVec))  

    if (type == "continuous")
    {
        ## Ordinary least squares estimation
        estMethod <- drmEMls(dose, resp, multCurves2, startVecSc, robustFct, wVec, rmNA, dmf = dmatfct, 
        doseScaling = doseScaling, respScaling = respScaling)
        
#        if (adjust == "vp")  #(varPower)
#        {        
#            estMethod <- drmEMvp(dose, resp, multCurves2)  # mdrcVp(dose, resp, multCurves2)
#            lenStartVec <- length(startVec)
#            
#            start2ss <- estMethod$"ssfct"(cbind(dose, resp))
#            if (missing(start2))
#            {
#                startVec <- c(startVec, start2ss)            
#            } else {
#                if (length(start2) == 2)  # canonical 2?
#                {
#                    startVec <- c(startVec, start2)            
#                }
#            }
##            startVec <- c(startVec, estMethod$"ssfct"(cbind(dose, resp)))
#            parmVec <- c(parmVec, "Sigma", "Power")
#            
#            startVecSc <- startVec
#        }
                  
#        if (!is.null(vvar))
#        {
#            estMethod <- mdrcHetVar(dose, resp, multCurves2, vvar)
#            lenStartVec <- length(startVec)            
#            startVec <- c(startVec, estMethod$"ssfct"(cbind(dose, resp)))
#            parmVec <- c(parmVec, as.character(unique(vvar)))        
#        }
    }
    if (identical(type, "binomial"))
    {
        estMethod <- drmEMbinomial(dose, resp, multCurves2, startVecSc, robustFct, wVec, rmNA, 
        doseScaling = doseScaling)        
    } 
    if (identical(type, "Poisson"))
    {
        estMethod <- drmEMPoisson(dose, resp, multCurves2, startVecSc, weightsVec = wVec, 
                                  doseScaling = doseScaling)
    }
    if (identical(type, "event"))
    {
        estMethod <- drmEMeventtime(dose, resp, multCurves2, doseScaling = doseScaling)
    }  
#    if (identical(type, "standard"))
#    {
#        estMethod <- drmEMstandard(dose, resp, multCurves2, doseScaling = doseScaling)
#    }      
#    if (identical(type, "Wadley"))
#    {
#        estMethod <- drmEMWadley(dose, resp, multCurves2, doseScaling = doseScaling)
#        startVecSc <- c(startVecSc, max(resp) * 1.3)
#    }                  
    opfct <- estMethod$opfct            


    ## Re-fitting the ANOVA model to incorporate Box-Cox transformation (if necessary)
#    if (type == "continuous")
#    {    
#        if (!is.na(lambda))
#        {
#            dset <- data.frame(dose, doseFactor, resp, assayNo, bcc)  # dataset with new resp values        
#            anovaModel0 <- (testList$"anovaTest")(anovaFormula, dset)            
##            anovaModel <- anovaModel0$"anovaFit"
#        }
#    }

    ## Defining lower and upper limits of parameters
#    if (constrained)
#    {
    if (!is.null(lowerl)) 
    {
        if (!is.numeric(lowerl) || !((length(lowerl) == sum(ncclVec)) || (length(lowerl) == numNames)))
        {
            stop("Not correct 'lowerl' argument")
        } else {
            if (length(lowerl) == numNames) 
            {
                lowerLimits <- rep(lowerl, ncclVec)
            } else {
                lowerLimits <- lowerl
            }
        }
        constrained <- TRUE 
        
    } else {  ## In case lower limits are not specified
        lowerLimits <- rep(-Inf, length(startVec))
    }

    if (!is.null(upperl)) 
    {
        if (!is.numeric(upperl) || !((length(upperl) == sum(ncclVec)) || (length(upperl) == numNames)))
        {
            stop("Not correct 'upperl' argument")
        } else {
            if (length(upperl) == numNames) 
            {
                upperLimits <- rep(upperl, ncclVec)
            } else {
                upperLimits <- upperl
            }
        } 
        constrained <- TRUE
                
    } else {  ## In case upper limits are not specified
        upperLimits <- rep(Inf, length(startVec))
    }
    
    lowerLimits <- lowerLimits  / longScaleVec
    upperLimits <- upperLimits  / longScaleVec
        
#    if (all(!is.finite(lowerLimits)) && all(!is.finite(upperLimits))) 
#    {
#        stop("No constraints are imposed via 'lowerl' and 'upperl' arguments")
#    }
#    }

    ## Optimising
    
    ## Setting derivatives
    opdfctTemp <- estMethod$"opdfct1"
    appFct <- function(x, y){tapply(x, y, sum)}   
    
    if (!is.null(opdfctTemp))
    {
        opdfct1 <- function(parm)
        {
#            print(as.vector(apply(opdfctTemp(parm), 2, appFct, assayNo)))
            as.vector(apply(opdfctTemp(parm), 2, appFct, assayNo))
        }
    } else {
        opdfct1 <- NULL
    }  

    ## Manipulating before optimisation
        
#    ## Scaling x values
#if (FALSE)
#{
#    sxInd <- fct$"sxInd"
#    sxYN <- !is.null(sxInd) && ((max(dose)<1e-2) || (min(dose)>1e2) || (diff(range(dose))>1e2) )
#    if ( sxYN && (is.null(fctList)) )
#    {
##        if (!is.null(fctList))
##        {
##            parmIndX <- rep(0, numAss)
##            for (i in 1:numAss)
##            {
##                parmIndX[i] <- fctList[[i]]$"sxInd"
##            }
##            parmIndX <- cumsum(parmIndX)
##        } else {
#            parmIndX <- parmPos[sxInd] + 1:ncclVec[sxInd]        
##        }
#    
#        scaleXConstant <- median(dose)
#        sxFct <- scaleX(scaleXConstant)  # , scaleX(dose, maxDose)        
#        if (adjust == "vp")
#        {
#            dose <- sxFct(dose)    
#            opfct <- drmEMvp(dose, resp, multCurves2)$"opfct"
#        }
#        
#        startVec[parmIndX] <- sxFct(startVec[parmIndX])
#    }
##    print(startVec)  # 2
#}

#    ## Scaling y values
#    ##  based on the original response value
#    ##  not the transformed values
#    syInd <- fct$"syInd"
#    lensy <- length(syInd)
#    parmIndY <- list()
#    
#    lyLim <- 1e-2
#    uyLim <- 1e2
#    syYN <- !is.null(syInd) && ((max(origResp)<lyLim) || (min(origResp)>uyLim) || (diff(range(origResp))>uyLim)) 
#    if ( syYN && (is.null(fctList)) )
#    {
##        if (!is.null(fctList))
##        {
##            parmIndY <- rep(0, numAss)
##            for (i in 1:numAss)
##            {
##                parmIndY[[i]] <- fctList[[i]]$"syInd"
##            }
##            parmIndY <- cumsum(as.vector(unlist(parmIndY)))
##        } else {
#            for (i in 1:lensy)
#            {
#                parmIndY[[i]] <- parmPos[syInd[i]] + c(1:ncclVec[syInd[i]])
#            }
#            tempPIY <- as.vector(unlist(parmIndY))
#            parmIndY <- tempPIY
##        }
#        if (adjust == "bc1")
#        { 
#            scaleYConstant <- bcfct(median(origResp), lambda, bcTol)  # median(origResp)
#        } else {
#            scaleYConstant <- median(origResp)  
#        }        
#        syFct <- scaleY(median(origResp))  # scaleY(scaleYConstant)
#        startVec[parmIndY] <- syFct(startVec[parmIndY])
#    }
#    # scaling of y values through 'opfct' definition
##    print(startVec)  # 3


    ## Testing nonlinear function
#    print(startVecSc)
#    print(multCurves2(dose, startVecSc))
#    print(opfct(startVecSc))
##    print(dose)
##    print(resp)
     
    ## Scaling objective function
#    if (type == "continuous")
#    {
#        ofVal <- opfct(startVec)
#        if ( !is.nan(ofVal) && ( (ofVal < 1e-2) || (ofVal >1e2) ) )
#        {
#            opfct2 <- function(c){opfct(c)/opfct(startVec)}
#        } else {
#            opfct2 <- opfct
#        }
#    } else {
#        opfct2 <- opfct
#    }
#    opfct2 <- opfct  # only used once below

    ## Optimising the objective function previously defined
    startVecSc <- as.vector(startVecSc)  # removing names
    nlsFit <- drmOpt(opfct, opdfct1, startVecSc, optMethod, constrained, warnVal, 
    upperLimits, lowerLimits, errorMessage, maxIt, relTol, parmVec = parmVec, traceVal = control$"trace",
    matchCall = callDetail, silentVal = control$"otrace") 
#    matchCall = match.call()) 
        
    if (!nlsFit$convergence) {return(nlsFit)}
    
    if (identical(type, "event"))
    {
#        dose <- dose[isFinite, 2]
#        resp <- (as.vector(unlist(tapply(resp, assayNo, function(x){cumsum(x) / sum(x)}))))[isFinite]

#        orderDose <- order(dose0)
#        dose1 <- dose0[orderDose]

        assayNo0 <- assayNo[isFinite]
        dose0 <- dose[, 2]
        dose1 <- dose0[isFinite]
        dose <- as.vector(unlist(tapply(dose1, assayNo0, function(x){unique(sort(x))})))
        
        ## Rescaling per curve id
        idList <- split(data.frame(dose0, resp), assayNo)
#        print(idList)
        
        respFct <- function(idListElt)
        {
            doseVec <- idListElt[, 1]
            dose2 <- unique(sort(doseVec))
            orderDose <- order(doseVec)
            resp1 <- tapply(idListElt[orderDose, 2], doseVec[orderDose], sum)  # obtaining one count per time interval
            resp2 <- cumsum(resp1) / sum(resp1)
            
            cbind(dose2, resp2)[is.finite(dose2), , drop = FALSE]
        }
        drList <- lapply(idList, respFct)
        lapList <- lapply(drList, function(x){x[, 1]})
        dose <- as.vector(unlist(lapList))
        resp <- as.vector(unlist(lapply(drList, function(x){x[, 2]})))
        
#        listCI <- split(assayNoOld[isFinite], assayNoOld[isFinite])
#        splitFactor <- factor(assayNoOld[isFinite], exclude = NULL)
        splitFactor <- factor(assayNo, exclude = NULL)        
        listCI <- split(splitFactor, splitFactor)
        lenVec <- as.vector(unlist(lapply(lapList, length)))
#        print(listCI)
#        print(lenVec)
        plotid <- as.factor(as.vector(unlist(mapply(function(x,y){x[1:y]}, listCI, lenVec))))
#        plotid <- plotid[complete.cases(plotid)] 
        levels(plotid) <- unique(assayNoOld)
    } else {
        plotid <- NULL
    }

#    if (identical(type, "Wadley"))
#    {
#        longScaleVec <- c(longScaleVec, 1)
#    
#    }        
         
#    print(nlsFit) 
     
    ## Adjusting for pre-fit scaling 
    if (!is.null(scaleFct))
    {
        # Scaling the sums of squares value back
        nlsFit$value <- nlsFit$value * (respScaling^2)
   
        # Scaling estimates and Hessian back
        nlsFit$par <- nlsFit$par * longScaleVec
        nlsFit$hessian <- nlsFit$hessian * (1/outer(longScaleVec/respScaling, longScaleVec/respScaling))
    }
    
    if (!is.null(fct$"retFct"))
    {
        drcFct <- fct$"retFct"(1, 1)  #, numObs)  # resetting the scaling
        drcFct1 <- function(dose, parm)
        {
            drcFct(dose, parm2mat(parm)[isFinite, , drop = FALSE])
        }
    }
    
    
#    print(nlsFit$par)
#    nlsFit$value <- opfct(nlsFit$par)  # used in the residual variance

    ## Manipulating after optimisation
    
#    ## Adjusting for scaling of y values
#    if ( syYN && (is.null(fctList)) )
#    {
#        nlsFit$value <- syFct(syFct(nlsFit$value, down = FALSE), down = FALSE)
#        startVec[parmIndY] <- syFct(startVec[parmIndY], down = FALSE)
#        nlsFit$par[parmIndY] <- syFct(nlsFit$par[parmIndY], down = FALSE)
#
#        scaleFct1 <- function(hessian) 
#                     {
#                         newHessian <- hessian
#                         newHessian[, parmIndY] <- syFct(newHessian[, parmIndY], down = FALSE)
#                         newHessian[parmIndY, ] <- syFct(newHessian[parmIndY, ], down = FALSE)
#                         return(newHessian)
#                     }                
#    } else {
#        scaleFct1 <- function(x) {x}    
#    }

    
#    ## Adjusting for scaling of x values
#if (FALSE)
#{    
#    if ( sxYN && (is.null(fctList)) )  # (!is.null(sxInd))
#    {
#        if (adjust == "vp")
#        {
#            dose <- sxFct(dose, down = FALSE)
#        }
#        startVec[parmIndX] <- sxFct(startVec[parmIndX], down = FALSE)
#        nlsFit$par[parmIndX] <- sxFct(nlsFit$par[parmIndX], down = FALSE)
#
#        scaleFct2 <- function(hessian) 
#                     {
#                         newHessian <- scaleFct1(hessian)
#                         newHessian[, parmIndX] <- sxFct(newHessian[, parmIndX], down = FALSE)
#                         newHessian[parmIndX, ] <- sxFct(newHessian[parmIndX, ], down = FALSE)
#                         return(newHessian)
#                     }                
#    } else {
#        scaleFct2 <- function(hessian) 
#        {
#            scaleFct1(hessian)
#        }
#    }
#}

#    ## Handling variance parameters
#    varParm <- NULL
#        
#    if (varPower)
#    {
#        varParm <- list(type = "varPower", index = 1:lenStartVec)        
#    }
#    if (!is.null(vvar))
#    {
#        varParm <- list(type = "hetvar", index = 1:lenStartVec)
#    }

    # Testing against the ANOVA (F-test)
    nlsSS <- nlsFit$value
    nlsDF <- numObs - length(startVec)

    ## Constructing a plot function
        
    ## Picking parameter estimates for each curve. Does only work for factors not changing within a curve!
    if (!is.null(cm)) {iVec <- (1:numAss)[!(uniqueNames==cm)]} else {iVec <- 1:numAss}
    
    pickCurve <- rep(0, length(iVec))
    for (i in iVec)
    {
       pickCurve[i] <- (1:numObs)[assayNo == i][1]
    }
    parmMat <- matrix(NA, numAss, numNames)

    fixedParm <- (estMethod$"parmfct")(nlsFit)
#    print(nlsFit$par)
#    print(fixedParm)
    parmMat[iVec, ] <- (parm2mat(fixedParm))[pickCurve, ]

    indexMat2 <- parm2mat(1:length(fixedParm))
    indexMat2 <- indexMat2[!duplicated(indexMat2), ]

#    if(!is.null(fctList))
#    {
#         parmMat <- matrix(NA, numAss, numNames)
#         for (i in 1:numAss)
#         {
#             parmMat[i, ivList2[[i]]] <- fixedParm[(posVec[i]+1):posVec[i+1]]
#         }
#    }    
        
    if (!is.null(cm))
    {
#        conPos <- upperPos
#        print(conPos)
        parmMat[-iVec, upperPos] <- (parm2mat(fixedParm))[assayNoOld == cm, , drop = FALSE][1, upperPos]  
        # 1: simply picking the first row
    }
    rownames(parmMat) <- assayNames


    pmFct <- function(fixedParm)
    {
        if (!is.null(cm)) {iVec <- (1:numAss)[!(uniqueNames == cm)]} else {iVec <- 1:numAss}
    
        if (!is.null(cm))
        {
#            conPos <- conList$"pos"
            parmMat[-iVec, upperPos] <- (parm2mat(fixedParm))[assayNoOld == cm, , drop = FALSE][1, upperPos]  
            # 1: simply picking the first row
        }
        rownames(parmMat) <- assayNames
    
        return(parmMat)
    }
    parmMat <- pmFct(fixedParm)  # (estMethod$"parmfct")(nlsFit) )
#    print(pmFct(1:length(fixedParm)))

#    ## Scaling parameters
#    if (!is.null(fct$scaleFct))
#    {
#        scaleFct <- function(parm)
#        {
#            fct$scaleFct(parm, xScaling, yScaling)
#        }
#    
#        parmMat <- apply(parmMat, 1, scaleFct)
#    }
#

    ## Constructing design matrix allowing calculations for each curve
#    colPos <- 1
#    rowPos <- 1
#    Xmat <- matrix(0, numAss*numNames, length(nlsFit$par))
#    Xmat <- matrix(0, numAss*numNames, length(fixedParm))


#    if (!is.null(fctList)) {omitList <- list()}
#    for (i in 1:numNames)
#    {
#        indVec <- iVec
#        lenIV <- length(indVec)
#
#        nccl <- ncol(pmodelsList2[[i]])  # min(maxParm[i], ncol(pmodelsList2[[i]]))        
#
#        XmatPart <- matrix(0, lenIV, nccl)
#        k <- 1
#        if (!is.null(fctList)) {omitVec <- rep(TRUE, lenIV)}
#        for (j in indVec)
#        {
#            if (!is.null(fctList))
#            {
#                parPresent <- !is.na(match(i, ivList2[[j]]))
#                omitVec[k] <- parPresent
#            }
#            
#            XmatPart[k, ] <- (pmodelsList2[[i]])[(1:lenData)[assayNo == j][1], 1:nccl]
#            k <- k + 1
#        }
#        if (!is.null(fctList))
#        {
#            XmatPart <- XmatPart[omitVec, , drop = FALSE]
#            nccl <- nccl - sum(!omitVec)
#            omitList[[i]] <- omitVec
#        }
#
#        Xmat[rowPos:(rowPos+lenIV-1), colPos:(colPos+nccl-1)] <- XmatPart
#        colPos <- colPos + nccl
#        rowPos <- rowPos + lenIV
#    }
#    Xmat <- Xmat[1:(rowPos-1), 1:(colPos-1)]


    ## Defining the plot function    
    pfFct <- function(parmMat)
    {
        plotFct <- function(dose)
        {
#            if (xDim == 1) {lenPts <- length(dose)} else {lenPts <- nrow(dose)}
            if (is.vector(dose)) 
            {
                lenPts <- length(dose)
            } else {
                lenPts <- nrow(dose)
            }
#            print(lenPts)
#            print(ciOrigLength)

            curvePts <- matrix(NA, lenPts, ciOrigLength)  # numAss)
            for (i in 1:numAss)
            {
#                if (!is.null(fctList)) 
#                {
#                    drcFct <- fctList[[i]]$"fct"
#                    numNames <- matList[[i]][2]    
#                }
            
                if (i %in% iVec)
                {
#                    parmChosen <- parmMat[i, ]
                    parmChosen <- parmMat[i, complete.cases(parmMat[i, ])]  # removing NAs 
#                    print(parmChosen)                   
                    
                    parmMat2 <- matrix(parmChosen, lenPts, numNames, byrow = TRUE)
#                    print(parmMat2)
                    curvePts[, ciOrigIndex[i]] <- drcFct(dose, parmMat2)
                } else { curvePts[, i] <- rep(NA, lenPts)}
            }
            return(curvePts)
        }
    
        return(plotFct)    
    }
#    print(parmMat)
    plotFct <- pfFct(parmMat)
#    plotFct(0:10)


    ## Computation of fitted values and residuals
    if (identical(type, "event"))
    {
        multCurves2 <- modelFunction(dose, parm2mat, drcFct, cm, assayNoOld, upperPos, fct$"retFct", doseScaling, respScaling, isFinite)
    }
    predVec <- multCurves2(dose, fixedParm)        
    resVec <- resp - predVec
    resVec[is.nan(predVec)] <- 0    

    diagMat <- matrix(c(predVec, resVec), length(dose), 2)
    colnames(diagMat) <- c("Predicted values", "Residuals")


    ## Adjusting for robust estimation: MAD based on residuals, centered at 0, is used as scale estimate    
    if (robust%in%c("median", "trimmed", "tukey", "winsor"))
    {
        nlsFit$value <- (mad(resVec, 0)^2)*nlsDF 
    }
#    if (robust=="winsor")
#    {
#        K <- 1 + length(startVec)*var(psi.huber(resVec/s, deriv=1))
#    }
    if (robust%in%c("lms", "lts"))  # p. 202 i Rousseeuw and Leroy: Robust Regression and Outlier Detection
    {  
        scaleEst <- 1.4826*(1+5/(numObs-length(nlsFit$par)))*sqrt(median(resVec^2))                                 
        w <- (resVec/scaleEst < 2.5)
        nlsFit$value <- sum(w*resVec^2)/(sum(w)-length(nlsFit$par))    
    }

    
    ## Adding meaningful names for robust methods
    robust <- switch(robust, median="median", trimmed="metric trimming", tukey="Tukey's biweight", 
                             winsor="metric Winsorizing", lms="least median of squares",
                             lts="least trimmed squares")


    ## Collecting summary output
    sumVec <- c(NA, NA, NA, nlsSS, nlsDF, numObs)  # , alternative)
    sumList <- list(lenData = numObs, 
    alternative = NULL,  # alternative, 
    df.residual = numObs - length(startVec))


    ## The function call
#    callDetail <- match.call()
#    if (is.null(callDetail$fct)) {callDetail$fct <- substitute(l4())}


    ## The data set
    if (!is.null(logDose)) 
    {
        dose <- origDose
    }
    dataSet <- data.frame(origDose, origResp, assayNo, assayNoOld, wVec)
    
#    print(varNames0)
    if (identical(type, "event"))
    {
        names(dataSet) <- c(varNames0[c(2, 3, 1)], anName, anName, "weights")
    } else {
        names(dataSet) <- c(varNames0[c(2, 1)], anName, anName, "weights")
    }


#    ## Box-Cox information
#    bcVec <- c(lambda, boxcoxci)
#    if (all(is.na(bcVec))) {bcVec <- NULL}
#    if (!is.null(bcVec)) {bcVec <- c(bcVec, bcAdd)}


    ## Evaluating goodness-of-fit test
#    if (!is.null(gofTest)) {gofTest <- gofTest(resp, weights, predVec, sumList$"df.residual")}


#    ## Adjusting in case 'fctList' is specified
#    if (!is.null(fctList))
#    {
#        omitAllVec <- as.vector(unlist(omitList))
#
#        parmVec <- parmVec[omitAllVec] 
#        parmVecA <- parmVecA[omitAllVec] 
#        parmVecB <- parmVecB[omitAllVec] 
#        
#        orderVec <- match(as.vector(parmMat), nlsFit$par)
#        orderVec <- orderVec[complete.cases(orderVec)]       
#        
#        nlsFit$par <- nlsFit$par[orderVec]
#        nlsFit$hessian <- nlsFit$hessian[orderVec, orderVec]
#    }


    ## Constructing an index matrix for use in ED and SI
# (commented out Dec 7 2011, replaced by definition below of the index matrix)
#    hfct1 <- function(x)  # helper function
#    {
#        uniVec <- unique(x[!is.na(x)])
#        rv <- rep(NA, length(x))
#        for (i in 1:length(uniVec))
#        {
#            rv[abs(x-uniVec[i]) < 1e-12] <- i
#        }
#        rv
#    }
#    hfct2 <- function(x)
#    {
#        length(unique(x))
#    }
##    parmMat <- t(parmMat)
#    mat1 <- t(apply(t(parmMat), 1, hfct1))  # , 1:ncol(parmMat)))
#    cnccl <- head(cumsum(ncclVec), -1)
##    mat2 <- mat1
#    if (nrow(mat1) == 1) {mat1 <- t(mat1)}  # in case of only one curve
#    mat1[-1, ] <- mat1[-1, ] + cnccl

    ## Matrix of first derivatives evaluated at the parameter estimates
    if (isDF)
    {
#        print((parmMat[assayNo, , drop = FALSE])[isFinite, , drop = FALSE])
        deriv1Mat <- fct$"deriv1"(dose, (parmMat[assayNo, , drop = FALSE])[isFinite, , drop = FALSE])
    } else {
        deriv1Mat <- NULL
    }
#    deriv1Mat <- NULL

    ## Box-Cox information
    if (!is.null(bcVal))
    {
        bcVec <- list(lambda = bcVal, ci = c(NA, NA), bcAdd = bcAdd)
    } else {
        bcVec <- NULL
    }

    ## Parameter estimates
    coefVec <- nlsFit$par
    names(coefVec) <- parmVec
    
    ## Constructing the index matrix
#    parmMat <- t(parmMat)
    indexMat <- apply(t(parmMat), 2, function(x){match(x, coefVec)})

    ## Constructing data list ... where is it used?
    wName <- callDetail[["weights"]]
    if (is.null(wName)) 
    {
        wName <- "weights"
    } else {
        wName <- deparse(wName)
    }
#    dataList <- list(dose = as.vector(origDose), origResp = as.vector(origResp), weights = wVec, 
    dataList <- list(dose = origDose, origResp = as.vector(origResp), weights = wVec, 
    curveid = assayNoOld, resp = as.vector(resp),
    names = list(dName = varNames[1], orName = varNames[2], wName = wName, cNames = anName, rName = ""))
    if (identical(type, "event"))
    {
        dataList <- list(dose = dose, origResp = resp, weights = wVec[isFinite], 
        curveid = assayNoOld[isFinite], plotid = plotid, resp = resp,
        names = list(dName = varNames[1], orName = varNames[2], wName = wName, cNames = anName, rName = ""))
    }


    ## What about naming the vector of weights?

    ## Returning the fit
#    returnList <- list(varParm, nlsFit, list(plotFct, logDose), sumVec, startVec, list(parmVec, parmVecA, parmVecB), 
    returnList <- list(NULL, nlsFit, list(plotFct, logDose), sumVec, startVecSc * longScaleVec, 
#    returnList <- list(nlsFit, list(plotFct, logDose), sumVec, startVecSc * longScaleVec, 
    list(parmVec, parmVecA, parmVecB), 
    diagMat, callDetail, dataSet, t(parmMat), fct, robust, estMethod, numObs - length(startVec), 
#    anovaModel0, gofTest, 
#    sumList, NULL, pmFct, pfFct, type, mat1, logDose, cm, deriv1Mat, 
    sumList, NULL, pmFct, pfFct, type, indexMat, logDose, cm, deriv1Mat, 
    anName, data, wVec, 
    dataList,
    coefVec, bcVec,
    indexMat2)
    
    names(returnList) <- c("varParm", "fit", "curve", "summary", "start", "parNames", "predres", "call", "data", 
#    names(returnList) <- c("fit", "curve", "summary", "start", "parNames", "predres", "call", "data", 
    "parmMat", "fct", "robust", "estMethod", "df.residual", 
#    "anova", "gofTest", 
    "sumList", "scaleFct", "pmFct", "pfFct", "type", "indexMat", "logDose", "cm", "deriv1",
    "curveVarNam", "origData", "weights",
    "dataList", "coefficients", "boxcox", "indexMat2")
    ## Argument "scaleFct" not used anymore
    class(returnList) <- c("drc")  # , class(fct))

    return(returnList)
}
