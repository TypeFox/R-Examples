if (!isGeneric("show"))
  setGeneric("show", function(object)
             standardGeneric("show"))
if (!isGeneric("summary"))
  setGeneric("summary", function(object, ...)
             standardGeneric("summary"))
if (!isGeneric("plot"))
        setGeneric("plot", function(x, y, ...)
                standardGeneric("plot"))

## ################## ##
## INTERNAL CONSTANTS ##
## ################## ##

BETA2_MIN_CONST <- 0.01
TAU2_MIN_CONST <- 0.01
DENSITY_CONST <- 1e-300

## ######################## ##
## MAIN FUNCTIONS (EXPORTED)##
## ######################## ##

runScampi <- function(peptides, proteins, edgespp, rescaling=TRUE,
                      method="all", quantifyPeptides=TRUE, numIter=10,
                      verbose=FALSE)
{
  ## Purpose: Estimate a protein abundance score for each protein in the
  ##          dataset, based on the input peptide abundance scores and the
  ##          connectivity information between peptides and proteins.
  ##          Optionally, the peptide abundances can be estimated as well
  ##          to compare the predicted values with the input measurements.
  ## ----------------------------------------------------------------------
  ## Arguments: peptides:         Data frame with peptide information. The
  ##                              following columns are required: 'pepId'
  ##                              (unique identification number for each
  ##                              distinct peptide sequence, numbering from
  ##                              1:n where n=number of distinct peptide
  ##                              sequences), 'pepSeq' (peptide sequence,
  ##                              optionally including modifications and
  ##                              charge states), and 'pepQty' (peptide
  ##                              abundance score). An additional column
  ##                              'pepObs' (peptide observability or
  ##                              identification score) is used if provided.
  ##                              Each row in the data frame describes one
  ##                              observed distinct peptide sequence.
  ##            proteins:         Data frame with the protein information.
  ##                              The following columns are required:
  ##                              'protId' (unique identification number for
  ##                              each distinct protein sequence, numbering
  ##                              from (n+1):(n+m) where m=number of distinct
  ##                              protein sequences), 'protName' (protein
  ##                              identifier or protein sequence). Each row
  ##                              describes a distinct protein sequence to
  ##                              which at least one of the observed peptides
  ##                              is matching.
  ##            edgespp:          Data frame with two mandatory columns:
  ##                              'pepId' and 'protId'. Each row defines an
  ##                              edge of the bipartite graph.
  ##            rescaling:        If TRUE, the peptide abundance scores are
  ##                              logarithmized (log10). If this
  ##                              transformation has not yet been done during
  ##                              preprocessing, it is strongly recommended
  ##                              to stick to the default: 'rescaling=TRUE'.
  ##            method:           Describes which method should be used for
  ##                              the parameter estimation. Available:
  ##                              'method="LSE"', 'method="MLE"' and
  ##                              'method="all"' (default).
  ##            quantifyPeptides: If 'TRUE' (default) do also re-quantify
  ##                              the peptides and assess the peptide
  ##                              abundance scores.
  ##            numIter:          Only used with 'method="MLE"', number of
  ##                              successful MLE optimizytion runs to
  ##                              perform.
  ##            verbose:          If 'TRUE', detailed output is provided.
  ## ----------------------------------------------------------------------
  ## Author: Sarah Gerster, Date: 15 Aug 2012, 08:30

  cl <- match.call()

  ## generate object of class 'scampi'
  dataIn <- scampi(peptides=peptides, proteins=proteins, edgespp=edgespp)

  ## check input data
  dataChecked <- checkInputData(dataIn, rescaling=rescaling,
                                verbose=verbose)
    
  ## preprocess input data
  tmpPrepro <- preprocessInputData(dataChecked, verbose=verbose)
  dataPrepro <- tmpPrepro[["dataPrepro"]]
  ppGraph <- tmpPrepro[["ppGraph"]]
  ccList <- tmpPrepro[["ccList"]]
  rm(tmpPrepro)

  ## estimate model parameters
  scampiParam <- estimateModelParameters(method=method, ccList,
                                         dataPrepro@peptides,
                                         numIter=numIter,
                                         verbose=verbose)

  ## compute protein and peptide abundance scores
  scampiRes <- quantifyProteins(dataPrepro, ccList, scampiParam,
                                quantifyPeptides=quantifyPeptides,
                                verbose=verbose)

  ## prepare output object
  scampiOut <- scampiVal(call = cl,
                         peptides = scampiRes@peptides,
                         proteins = scampiRes@proteins,
                         edgespp = scampiRes@edgespp,
                         parameters = scampiParam,
                         ppGraph = ppGraph,
                         ccList = ccList)
  
  return(scampiOut)
}


iterateScampi <- function(peptides, proteins, edgespp, rescaling=TRUE,
                          method="LSE", numIter=2, numMLEIter=10,
                          thresh=2, verbose=FALSE)
{
  ## Purpose: Estimate a protein abundance score for each protein in the
  ##          dataset, based on the input peptide abundance scores and the
  ##          connectivity information between peptides and proteins. The
  ##          expected values for the peptide abundances are computed as
  ##          well. Comparing these values with the initial measurements
  ##          allows to detect outliers in the input data. Several
  ##          iterations of abundance estimation and outlier removal can
  ##          then be performed.
  ## ----------------------------------------------------------------------
  ## Arguments: peptides:         Data frame with peptide information. The
  ##                              following columns are required: 'pepId'
  ##                              (unique identification number for each
  ##                              distinct peptide sequence, numbering from
  ##                              1:n where n=number of distinct peptide
  ##                              sequences), 'pepSeq' (peptide sequence,
  ##                              optionally including modifications and
  ##                              charge states), and 'pepQty' (peptide
  ##                              abundance score). An additional column
  ##                              'pepObs' (peptide observability or
  ##                              identification score) is used if provided.
  ##                              Each row in the data frame describes one
  ##                              observed distinct peptide sequence.
  ##            proteins:         Data frame with the protein information.
  ##                              The following columns are required:
  ##                              'protId' (unique identification number for
  ##                              each distinct protein sequence, numbering
  ##                              from (n+1):(n+m) where m=number of distinct
  ##                              protein sequences), 'protName' (protein
  ##                              identifier or protein sequence). Each row
  ##                              describes a distinct protein sequence to
  ##                              which at least one of the observed peptides
  ##                              is matching.
  ##            edgespp:          Data frame with two mandatory columns:
  ##                              'pepId' and 'protId'. Each row defines an
  ##                              edge of the bipartite graph.
  ##            rescaling:        If TRUE, the peptide abundance scores are
  ##                              logarithmized (log10). If this
  ##                              transformation has not yet been done during
  ##                              preprocessing, it is strongly recommended
  ##                              to stick to the default: 'rescaling=TRUE'.
  ##            method:           Describes which method should be used for
  ##                              the parameter estimation. Available:
  ##                              'method="LSE"', 'method="MLE"' and
  ##                              'method="all"' (default).
  ##            numIter:          Number of estimation/outlier-removal
  ##                              iterations to be performed.
  ##            numMLEIter:       Only used with 'method="MLE"', number of
  ##                              successful MLE optimizytion runs to
  ##                              perform.
  ##            thresh:           Constant to tune the outlier selection
  ##                              process (IQR criterion).
  ##            verbose:          If 'TRUE', detailed output is provided.
  ## ----------------------------------------------------------------------
  ## Author: Sarah Gerster, Date: 16 Aug 2012, 14:10

  scampiOut <-
    list(iteration1 = list(scampiRes=runScampi(peptides, proteins, edgespp,
                             rescaling=rescaling, method=method,
                             quantifyPeptides=TRUE, numIter=numMLEIter,
                             verbose=verbose)))
  iterCount <- 2
  
  while (iterCount <= numIter) {
    peptides <- scampiOut[[paste("iteration", iterCount-1,
                                 sep="")]][["scampiRes"]]@peptides
    proteins <- scampiOut[[paste("iteration", iterCount-1,
                                 sep="")]][["scampiRes"]]@proteins
    edgespp <- scampiOut[[paste("iteration", iterCount-1,
                                sep="")]][["scampiRes"]]@edgespp
    
    ## select outlying peptides
    peptides[,"resid"] <- 
      peptides[,"pepQty"] - peptides[,paste(method,"Score",sep="")]
    
    qt1 <- quantile(peptides[,"resid"], 0.25) 
    qt3 <- quantile(peptides[,"resid"], 0.75) 
    iqr <- qt3 - qt1
    ind.remove.top <-  which(sort(peptides[,"resid"], decreasing=TRUE) > 
                             (qt3 + thresh*iqr))
    ind.remove.bot <-  which(sort(peptides[,"resid"]) < (qt1 - thresh*iqr))

    if ((length(ind.remove.bot)==0) && (length(ind.remove.top)==0)) {
      message("No outliers to remove after ", iterCount-1,
              " iterations. Stopping now.")
      break
    }
    
    peptides[,"outlier"] <- FALSE
    if (!(length(ind.remove.bot)==0)) {
      peptides[order(peptides[,"resid"])[ind.remove.bot],"outlier"] <- TRUE
    }
    if (!(length(ind.remove.top)==0)) {
      peptides[order(peptides[,"resid"], 
                     decreasing=TRUE)[ind.remove.top],"outlier"] <- TRUE
    }

    ## remove outliers from peptide dataframe
    peptides.outlier <- peptides[peptides[,"outlier"],
                                 c("pepSeq", "pepQty", "pepObs")]
    peptides.new <- peptides[!peptides[,"outlier"],]
    nPeptides <- nrow(peptides.new)
    
    ## get rid of surplus edges 
    edgespp.new <- edgespp[edgespp[,"pepId"] %in% peptides.new[,"pepId"],]
    
    ## get rid of surplus proteins
    proteins.new <- proteins[proteins[,"protId"] %in% edgespp.new[,"protId"],]
    nProteins <- nrow(proteins.new)
    
    ## renumber
    peptides.new[,"newId"] <- 1:nPeptides
    proteins.new[,"newId"] <- (nPeptides+1):(nPeptides+nProteins)
    edgespp.new[,"pepId"] <- sapply(edgespp.new[,"pepId"], replacePepNumbers,
                                    peptides.new)
    edgespp.new[,"protId"] <- sapply(edgespp.new[,"protId"], replaceProtNumbers,
                                     proteins.new)
    peptides.new[,"pepId"] <- peptides.new[,"newId"]
    proteins.new[,"protId"] <- proteins.new[,"newId"]
    peptides <- peptides.new[,-ncol(peptides.new)]
    proteins <- proteins.new[,-ncol(proteins.new)]
    edgespp <- edgespp.new

    ## re-run scampi
    scampiOut[[paste("iteration", iterCount, sep="")]] <-
              list(scampiRes=runScampi(peptides, proteins, edgespp,
                     rescaling=FALSE, method=method,
                     quantifyPeptides=TRUE, numIter=numMLEIter,
                     verbose=verbose),
                   peptideOutliers=peptides.outlier)
    
    ## increase counter to keep track of iterations
    iterCount <- iterCount + 1 
  }

  return(scampiOut)
}

## ################## ##
## EXPORTED FUNCTIONS ##
## ################## ##

checkInputData <- function(scampiData, ...) UseMethod("checkInputData")

checkInputData.scampi <- function(scampiData, rescaling=TRUE,
                                  verbose=FALSE, ...) 
{
  ## Purpose: Perform some consistency checks on the input data to ensure
  ##          that the protein quantification and, optionally, the peptide
  ##          reassessment will work smoothly.
  ## ----------------------------------------------------------------------
  ## Arguments: scampiData: object of class 'scampi'
  ##            rescaling:  if TRUE, transform peptide abundance scores by
  ##                        log10()
  ##            verbose:    if not FALSE, provide minimal progress info  
  ## ----------------------------------------------------------------------
  ## Author: Sarah Gerster, Date: 10 Aug 2012, 10:23

  if (verbose)
    message("*** 'checkInputData.scampi' checks input data frames for ",
            "consistency ***")

  peptides <- scampiData@peptides
  proteins <- scampiData@proteins
  edgespp <- scampiData@edgespp

  ## check if all required variables are present
  stopifnot(all(c("pepId", "pepSeq", "pepQty") %in%
                colnames(peptides)))
  stopifnot(all(c("protId", "protName") %in% colnames(proteins)))
  stopifnot(all(c("pepId", "protId") %in% colnames(edgespp)))
  
  ## peptides: additional specific checks
  if(length(colnames(peptides)) >
     (ifelse("pepObs" %in% colnames(peptides), 4, 3))){
    message("  ** only the following variables of the peptides data ",
            "frame are used for the computations: 'pepId', 'pepSeq', ",
            "'pepObs' and 'pepQty'")
  }
  stopifnot(length(unique(peptides[,"pepSeq"])) ==
            length(peptides[,"pepSeq"]))
  stopifnot(length(unique(peptides[,"pepId"])) ==
            length(peptides[,"pepId"]))
  if (any(!is.finite(peptides[,"pepQty"])))
    stop("check your peptide abundances, they hold NA, NaN, Inf or -Inf")
  if ("pepObs" %in% colnames(peptides)) {
    if (any(!is.finite(peptides[,"pepObs"])))
      stop("check your peptide observabilities, ",
           "they hold NA, NaN, Inf or -Inf")
  } else {
    peptides[,"pepObs"] <- 1
  }
  if (rescaling) {
    if (any(peptides[,"pepQty"] <= 0))
      stop("check your peptide abundances, they contain zero or values < 0")
    peptides[,"pepQty"] <- log10(peptides[,"pepQty"])
  }
  nPeptides <- nrow(peptides)
  
  ## proteins: additional specific checks
  if(length(colnames(proteins)) > 2 && verbose) 
    message("  ** only the following variables of the proteins data ",
            "frame are used for the computations: 'protId' and 'protName'")
  stopifnot(length(unique(proteins[,"protName"])) ==
            length(proteins[,"protName"]))
  stopifnot(length(unique(proteins[,"protId"])) ==
            length(proteins[,"protId"]))
  nProteins <- nrow(proteins)

  ## edgespp: additional specific checks 
  if(length(colnames(edgespp)) > 2 && verbose) 
    message("  ** only the following variables of the edgespp data ",
            "frame are used for the computations: 'pepId' and 'protId'")
  stopifnot(nrow(unique(edgespp)) == nrow(edgespp))

  ## check numbering of the peptides and proteins
  if (sum(unique(as.numeric(edgespp[,"pepId"])) <= nPeptides) !=
      nPeptides)
    stop("error in 'pepId'")
  if (sum(unique(as.numeric(edgespp[,"protId"])) > (nPeptides) &
          unique(as.numeric(edgespp[,"protId"])) <=
          (nPeptides + nProteins))
      != nProteins)
    stop("error in 'protId'")

  ## check that each peptide and each protein occurs in at least one edge
  stopifnot(all(peptides[,"pepId"] %in% edgespp[,"pepId"]))
  stopifnot(all(proteins[,"protId"] %in% edgespp[,"protId"]))

  ## check that all edge indicies are included in the peptides and
  ## proteins data frames, respectively
  stopifnot(all(edgespp[,"pepId"] %in% peptides[,"pepId"]))
  stopifnot(all(edgespp[,"protId"] %in% proteins[,"protId"]))

  if (verbose)
    message("--> consistency check completed successfully <--")
  
  return(scampi(peptides=peptides, proteins=proteins, edgespp=edgespp))
}
          
preprocessInputData <- function(scampiData, ...) UseMethod("preprocessInputData")

preprocessInputData.scampi <- function(scampiData, verbose=FALSE, ...) 
{
  ## Purpose: Preprocess the graph structure and connected components
  ##          as much as possible to speed up computations for parameter
  ##          estimation and peptide/protein abundance predictions.
  ## ----------------------------------------------------------------------
  ## Arguments: scampiData: object of class 'scampi'
  ##            verbose:    if not FALSE, provide minimal progress info  
  ## ----------------------------------------------------------------------
  ## Dependencies: package  'RBGL'
  ##               function 'buildBipartiteGraph()'
  ##               function 'preprocessCCs()'
  ##               function 'getMyCCNr()'
  ## ----------------------------------------------------------------------
  ## Author: Sarah Gerster, Date: 10 Aug 2012, 10:23

  if (verbose)
    message("*** 'preprocessInputData.scampi' performs some pre-",
            "calculations on the (sub-)graphs to speed up further ",
            "computations ***")
  
  peptides <- scampiData@peptides
  proteins <- scampiData@proteins
  edgespp <- scampiData@edgespp
  nPeptides <- nrow(peptides)             # nb peps in data set
  nProteins <- nrow(proteins)             # nb prots in data set

  if (verbose)
    message("  ** preprocess graph")
  ## build the bipartite graph
  ppGraph <- buildBipartiteGraph(nPeptides, nProteins, edgespp, verbose)
  
  ## get a list of all the connected components
  ccList <- connectedComp(ppGraph)
  
  ## preprocess some parameter independent values for each connected
  ## component
  ccList <- lapply(ccList, preprocessCC, peptides,
                          edgespp, nPeptides, verbose)
  if (verbose)
    message("  ** preprocess proteins")
  ## preprocess some parameter independent values for each protein
  proteins[,"ccInd"] <-                   
    sapply(proteins[,"protId"], getMyCCNr,
           group="proteins", ccList, verbose)
  
  if (verbose)
    message("  ** preprocess peptides")
  ## preprocess some parameter independent values for each peptide  
  peptides[,"ccInd"] <-                   # prepro peptides
    sapply(peptides[,"pepId"], getMyCCNr, group="peptides", 
           ccList, verbose)
  peptides[,"nrNeiProt"] <- 
    sapply(peptides[,"pepId"], getMyNeighborhoodSize, 
           group="peptides", edgespp, verbose)

  if (verbose)
    message("--> data preprocessing completed successfully <--")
  
  return(list(dataPrepro = scampi(peptides=peptides, proteins=proteins,
                edgespp=edgespp),
              ppGraph = ppGraph, ccList = ccList))
}

estimateModelParameters <- function(method="all", ccList,
                                    peptides=NULL, numIter=10,
                                    verbose=FALSE) 
{
  ## Purpose: Estimate the model paramteres alpha, beta, mu and tau from
  ##          the data
  ## ----------------------------------------------------------------------
  ## Arguments: method:       method to be used for the parameter esti-
  ##                          mation; can be 'all' (MLE and LSE), 'LSE' or
  ##                          'MLE'
  ##            ccList: list of pre-processed connected components
  ##            peptides:     data frame with pre-processed peptide info
  ##                          (only used for LSE)
  ##            numIter:      number of successful munerical optimizations
  ##                          to perform (only used for MLE)
  ##            verbose:      if TRUE: print minimal progress info, if >1:
  ##                          print (much) more info for sub-processes
  ## ----------------------------------------------------------------------
  ## Dependencies: function sampleInitialVector()
  ## ----------------------------------------------------------------------
  ## Author: Sarah Gerster, Date: 14 Aug 2012, 15:00

  if (verbose)
    message("*** 'estimateModelParameters.scampi' estimates the model ",
            "parameter with the method(s) of choice (provided by the ",
            " argument 'method') ***")

  res <- list()
  
  if (method %in% c("all", "LSE")) {
    lseParam <- estimateLSEparam(peptides, ccList, verbose)
    res[["LSE"]] <- lseParam
  }

  if (method %in% c("all", "MLE")) {
    mleParamMat <- matrix(NA, ncol=5, nrow=numIter)
    colnames(mleParamMat) <- c("alphaH","betaH","muH","tauH", "fnval")
    runIter <- 1

    if (method == "all") {
      run <- try(estimateMLEparam(c(lseParam[1], lseParam[2], lseParam[3],
                                    lseParam[4]), ccList, verbose),
                 TRUE)
      ## if run was successful, store it
      if (length(run) == 5) {
        mleParamMat[runIter,] <- run
        runIter <- runIter + 1
      }
    }
    
    while (runIter <= numIter) {
      if (verbose)
        message("--> next trial with new random starting values")
      ## run optimization process with a random start values
      run <- try(estimateMLEparam(sampleInitialVector(), ccList,
                                  verbose), TRUE)
      ## if run was successful, store it
      if (length(run) == 5) {
        mleParamMat[runIter,] <- run
        if (verbose)
          message("      (finished ", runIter, "successful optimization ",
                  "iterations)")
        runIter <- runIter + 1
      }
    }

    ## select run with minimal functrion value
    mleParam <- mleParamMat[which.min(mleParamMat[,5]),]
    res[["MLE"]] <- mleParam
  }    
  
  if (verbose)
    message("--> parameter estimated successfully <--")
  
  return(res)
}

quantifyProtein <- function(protInfo, ccList, param, verbose=FALSE)
{
  ## Purpose: compute the protein abundance score for the specified
  ##          parameter values
  ## ----------------------------------------------------------------------
  ## Arguments: protInfo:     vector with 2 elements: protId and index of 
  ##                          connected component which contains protId
  ##            ccList: list of pre-processed connected components
  ##            param:        vector with 4 elements: alpha, beta, mu and
  ##                          tau (model parameter estimates)
  ## ----------------------------------------------------------------------
  ## Author: Sarah Gerster, Date: 15 Aug 2012, 09:45

  protID <- as.numeric(protInfo[1])
  cc <- ccList[[as.numeric(protInfo[2])]]
  
  if (verbose > 1)
    message("  ** 'quantifyProtein' is computing expected abundance for ",
            "protein ", protID)
  
  ## which peptides have a common edge with protein protID?
  edgeToC <- as.numeric(cc$peptides) %in%
    cc$edges[cc$edges[,"protId"]==protID,"pepId"]
  ## compute covariance between abundance of protein 'protID' and
  ## peptide quantity: 'Gcu'
  Gcu <- cc$pepObs * param["betaH"] * edgeToC

  ## return expected value for the abundance of protein 'protID'
  if (is.matrix(cc$distance)){
    res <- crossprod(cc$pepQty - param["alphaH"] -
                     cc$pepObs*param["betaH"]*param["muH"]*diag(cc$distance),
                     solve(cc$covU)) %*% Gcu + param["muH"]
    protVar <- 1 - Gcu %*% solve(cc$covU) %*% Gcu
  } else {
    res <- crossprod(cc$pepQty - param["alphaH"] -
                     cc$pepObs*param["betaH"]*param["muH"]*cc$distance,
                     solve(cc$covU)) %*% Gcu + param["muH"]
    protVar <- 1 - Gcu %*% solve(cc$covU) %*% Gcu
  }
  
  if (verbose > 1)
    message("  --> protein quantified successfully")
  
  ## return(res)
  return(c(res, protVar))
}

quantifyProteins <- function(scampiData, ccList, paramList,
                             quantifyPeptides=FALSE, verbose=FALSE)
{
  ## Purpose: compute the protein (and optionally the peptide) abundance
  ##          score for the specified parameter values
  ## ----------------------------------------------------------------------
  ## Arguments: scampiData:       object of class scampi 
  ##            ccList:     list of pre-processed connected
  ##                              components
  ##            paramList         named list of vectors with 4 elements:
  ##                              alpha, beta, mu and tau (model parameter
  ##                              estimates);
  ##                              name of element = parameter estimation
  ##                              method
  ##            quantifyPeptides: if TRUE, also compute peptide abundance
  ##                              scores
  ## ----------------------------------------------------------------------
  ## Dependencies: function getCovU()
  ##               function quantifyProtein()
  ##               function quantifyPeptide()
  ## ----------------------------------------------------------------------
  ## Author: Sarah Gerster, Date: 15 Aug 2012, 09:45
  
  if (verbose)
    message("  ** 'quantifyProteins' is computing expected abundance for ",
            "all proteins in 'proteins'")

  proteins <- scampiData@proteins
  peptides <- scampiData@peptides
  
  for (methIter in names(paramList)) {
    ccList <-                                  # precompute covU
      lapply(ccList, getCovU, beta=paramList[[methIter]]["betaH"],
             tau=paramList[[methIter]]["tauH"])
  
  
    ## proteins[,paste(methIter, "Score", sep="")] <- # compute prot ab.
    ##   apply(proteins[, c("protId", "ccInd")], 1, quantifyProtein,
    ##         ccList, param=paramList[[methIter]], verbose=verbose)
    
    proteins[,c(paste(methIter, "Score", sep=""),    # compute prot ab.
                paste(methIter, "Var", sep=""))] <-   
                  t(apply(proteins[, c("protId", "ccInd")], 1,
                          quantifyProtein, ccList,
                          param=paramList[[methIter]], verbose=verbose))

    if (quantifyPeptides) {
      peptides[,paste(methIter, "Score", sep="")] <- # compute pep ab.
        apply(peptides[,c("pepId","ccInd")], 1, quantifyPeptide,
              ccList, param=paramList[[methIter]], verbose=verbose)
      peptides[,paste(methIter, "Resid", sep="")] <- # compute pep resid.
        peptides[,"pepQty"] - peptides[,paste(methIter, "Score", sep="")]
    }
  }
  
  if (verbose) {
    if (quantifyPeptides) {
      message("--> proteins and peptides quantified successfully <--")
    } else {
      message("--> proteins quantified successfully <--")
    }
  }
  
  return(scampi(peptides = peptides,
                proteins = proteins,
                edgespp = scampiData@edgespp))
}


quantifyPeptide <- function(pepInfo, ccList, param, verbose=FALSE)
{
  ## Purpose: compute the peptide abundance score for the specified
  ##          parameter values
  ## ----------------------------------------------------------------------
  ## Arguments: protInfo:     vector with 2 elements: 'pepId' and index of 
  ##                          connected component which contains 'pepId'
  ##            ccList: list of pre-processed connected components
  ##            param:        vector with 4 elements: 'alphaH', 'betaH', 
  ##                          'muH' and 'tauH' (model parameter estimates)
  ## ----------------------------------------------------------------------
  ## Author: Sarah Gerster, Date: 15 Aug 2012, 09:45
  pepID <- as.numeric(pepInfo[1])
  cc <- ccList[[as.numeric(pepInfo[2])]]
  
  if (verbose > 1)
    message("  ** 'quantifyPeptide' is computing expected abundance for ",
            "peptide ", pepID)
  
  if (length(cc$peptides) > 1) {
    pepInd <- which(as.numeric(cc$peptides)==pepID)
    
    ## estimate ProtConc without using pepID
    protSum <- 0
    covU <- cc$covU[-pepInd,-pepInd]
    U <- cc$pepQty[-pepInd]
    s <- cc$pepObs[-pepInd]
    d <- cc$distance[-pepInd,-pepInd]
    peps <- cc$peptides[-pepInd]
    edgeToP <- as.numeric(cc$proteins) %in%
      cc$edges[cc$edges[,"pepId"]==pepID,"protId"]
    
    ## if the protein is a neighbor of the peptide pepID, add its
    ## contribution to the peptide intensity
    for (protID in cc$proteins[edgeToP]) {
      edgeToC <- peps %in%
        cc$edges[cc$edges[,"protId"]==as.numeric(protID),"pepId"]
        Gcu <- s * param["betaH"] * edgeToC
        if (is.matrix(d)){
          protSum <- protSum + crossprod(U - param["alphaH"] -
                                         s*param["betaH"]*param["muH"]*diag(d),
                                         solve(covU)) %*% Gcu
        } else {
          protSum <- protSum + (U - param["alphaH"] -
                                s*param["betaH"]*param["muH"]*d)*Gcu/covU
        }
    }
    res <- param["alphaH"] + cc$pepObs[pepInd]*param["betaH"]*param["muH"]*
      cc$distance[pepInd,pepInd] + cc$pepObs[pepInd]*param["betaH"]*protSum
  } else {
    res <- param["alphaH"] + cc$pepObs*param["betaH"]*param["muH"]*cc$distance
  }

  names(res) <- NULL
  
  if (verbose > 1)
    message("  --> peptide quantified successfully")
  
  return(res)
}


## ################## ##
## INTERNAL FUNCTIONS ##
## ################## ##

buildBipartiteGraph <- function(nPeptides, nProteins, edgesPP,
                                verbose=FALSE)
{
  ## Purpose: Build a bipartite graph (with peptides and proteins). The
  ##          graph holds 'nPeptides' peptides (labeled from 1:nPeptides) 
  ##          and 'nProteins' (labeled from (nPeptides+1):(nPeptides+
  ##          nProteins)) connected according to the information in
  ##          'edgesPP'.
  ## ----------------------------------------------------------------------
  ## Arguments: nPeptides: Number of peptides in the dataset
  ##            nProteins: Number of proteins in the dataset
  ##            edgesPP:   Data frame with 2 columns: 'pepId' and 'protId'.
  ##                       Each row defines an edge of the bipartite graph.
  ## ----------------------------------------------------------------------
  ## Dependencies: package  'graph'
  ## ----------------------------------------------------------------------
  ## Author: Sarah Gerster, Date: 31 Aug 2011, 13:15

  if (verbose)
    message("  ** 'buildBipartiteGraph' sets up the input bipartite graph")
  
  ## prepare the edge list needed as input to 'graphNEL' based on the
  ## information in edgesPP
  Vpp <- as.character(1:(nPeptides+nProteins))
  edLpp <- vector("list", length=length(Vpp))
  names(edLpp) <- Vpp
  ## first set edges from peptides to proteins ...
  for (i in 1:nPeptides) {
    ind.i <- edgesPP[,"pepId"] == i
    edLpp[[i]] <- edgesPP[ind.i,"protId"]
  }
  ## ... then from proteins to peptides (undirected graph)
  for (j in (nPeptides+1):(nPeptides+nProteins)) {
    ind.j <- edgesPP[,"protId"] == j
    edLpp[[j]] <- edgesPP[ind.j,"pepId"]
  }

  return(new("graphNEL", nodes=Vpp, edgeL=edLpp, edgemode="undirected"))
}


preprocessCC <- function(cc, peptides, edgesPP, nPeptides, verbose=FALSE)
{
  ## Purpose: Prepare (parameter independent) data for the connected
  ##          component 'cc' in order to save time during subsequent
  ##          computations (parameter estimation and abundance prediction).
  ##          [cc := connected component]   
  ##          -> in the output, each cc has:
  ##             * elements (IDs of peptides & proteins)
  ##             * peptides (IDs of peptides)
  ##             * proteins (IDs of proteins)
  ##             * corresponding peptide intensities
  ##             * corresponding peptide scores
  ##             * "distance matrix" between peptides and proteins
  ##             * edges
  ## ----------------------------------------------------------------------
  ## Arguments: cc:        connected component
  ##            peptides:  data frame with the peptide information:
  ##                       -> requires columns: 'pepId', 'pepSeq',
  ##                          'pepObs' & 'pepQty'
  ##                       -> one row per observed peptide
  ##            edgesPP:   data frame with edge information
  ##                       -> requires columns: 'pepId' and 'protId'.
  ##                       Each row defines an edge of the bipartite graph.
  ##            nPeptides: number of peptides in the dataset
  ## ----------------------------------------------------------------------
  ## Dependencies: function 'getDmat()'
  ## ----------------------------------------------------------------------
  ## Author: Sarah Gerster, Date: 13 Aug 2012, 13:34

  if (verbose > 1)
    message("  ** 'preprocessCC' precomputes several values for 'cc'")
  
  res <- list()

  ## IDs of all elements in the cc
  res$elements <- as.numeric(cc)
  
  ## IDs of all proteins in the cc
  ind.prot <- as.integer(cc) > nPeptides
  res$proteins <- as.numeric(cc[ind.prot])
  
  ## IDs all the peptides in the cc
  res$peptides <- as.numeric(cc[!ind.prot])

  ## peptide intensities of all peptides in the cc
  ## (same order as in res$peptides)
  res$pepQty <- peptides[peptides[,"pepId"] %in% res$peptides,"pepQty"]
  
  ## peptide scores of all peptides in the cc
  ## (same order as in res$peptides)
  res$pepObs <- peptides[peptides[,"pepId"] %in% res$peptides,"pepObs"]
  
  ## "distance matrix": number of connected proteins for each peptide pair
  ## -> d[i,i] = |Ne(i)|
  ## -> d[i,j] = |Ne(i) AND Ne(j)|
  res$distance <- getDmat(res$peptides, edgesPP, verbose)

  ## prepare an element with the list of edges in 'cc'
  res$edges <- edgesPP[edgesPP[,"protId"] %in% as.numeric(res$proteins),]
  
  if (verbose > 1)
    message("  --> connected component successfully preprocessed")
  
  return(res)
}

getDmat <- function(pepArr, edgesPP, verbose=FALSE)
{
  ## Purpose: compute the "distance/connectivity matrix" for all peptide
  ##          pairs
  ##          -> d[i,i] = |Ne(i)|
  ##          -> d[i,j] = |Ne(i) AND Ne(j)|
  ##          => symmetric matrix
  ## ----------------------------------------------------------------------
  ## Arguments: pepArr:  IDs of the peptides
  ##            edgesPP: data frame with edge information
  ##                     -> requires columns: 'pepId' and 'protId'.
  ##                     Each row defines an edge of the bipartite graph.
  ## ----------------------------------------------------------------------
  ## Dependencies: function 'getMyProteins()'
  ## ----------------------------------------------------------------------
  ## Author: Sarah Gerster, Date: 13 Aug 2012, 14:02

  if (verbose > 1)
    message("    ** 'getDmat' sets up the connectivity/distance matrix")
  
  if (length(pepArr) == 1) {
    ## if there is only one peptide in the CC, return its number of
    ## neighboring proteins
    return(length(getMyProteins(pepArr, edgesPP)))
  } else {
    ## if there are more peptides, set up a matrix
    d <- matrix(0, nrow=length(pepArr), ncol=length(pepArr))
    ## compute entries for the upper triangular part
    ## -> d[i,j] = |Ne(i) AND Ne(j)|
    for (i in 1:(length(pepArr)-1)) {
      proti <- getMyProteins(pepArr[i],edgesPP)
      for (j in (i+1):length(pepArr)) {
        protj <- getMyProteins(pepArr[j],edgesPP)
        d[i,j] <- sum(proti %in% protj)
      }
    }
    ## put together the whole matrix (transpose triangular part and
    ## add diagonal elements (d[i,i] = |Ne(i)|))
    d <- d + t(d) +
      diag(unlist(lapply(lapply(pepArr, getMyProteins, edgesPP), length)))
  
  if (verbose > 1)
    message("    --> connectivity matrix successfully set up")
    
    return(d)
  }
}

getMyProteins <- function(pepID, edgesPP, verbose=FALSE)
{
  ## Purpose: Return the list of all proteins with an edge to peptide
  ##          'pepID'
  ## ----------------------------------------------------------------------
  ## Arguments: pepID:   ID of the peptide for which the list of matching
  ##                     proteins is needed
  ##            edgesPP: data frame with 2 columns: 'pepId' and 'protId'
  ##                     -> each row defines an edge of the bipartite graph
  ##            verbose: ???
  ## ----------------------------------------------------------------------
  ## Author: Sarah Gerster, Date: 13 Aug 2012, 14:18

  if (verbose > 1)
    message("    ** 'getMyProteins' finds the IDs of the proteins having a ",
            "common edge to peptide 'pepID'")
  
  return(edgesPP[edgesPP[,"pepId"]==pepID, "protId"])
}


getMyCCNr <- function(elID, group, ccList, verbose=FALSE)
{
  ## Purpose: Find connected component with 'ID' in its element 'group'
  ## ----------------------------------------------------------------------
  ## Arguments: elId:         index to look for (typically protein or
  ##                          peptide ID)
  ##            group:        element of 'connectedComp' in which we
  ##                          should look for 'elID'
  ##            ccList: list of pre-processed connected components
  ##            verbose:      ???
  ## ----------------------------------------------------------------------
  ## Dependencies: function isInCC()
  ## ----------------------------------------------------------------------
  ## Author: Sarah Gerster, Date:  13 Aout 2011, 17:40

  if (verbose > 1)
    message("    ** 'getMyCCNr' looks up the index of the connected ",
            "component holding 'elID' in its element 'group'")
  
  return(which(sapply(1:length(ccList), isInCC,
                      elID=elID, group=group, ccList)))
}


isInCC <- function(ccID, elID, group, ccList, verbose=FALSE)
{
  ## Purpose: Return TRUE if 'elID' belongs to element 'group' of
  ##          'ccList[[ccID]]'
  ## ----------------------------------------------------------------------
  ## Arguments: ccID:          connected component index 
  ##            elID:          index to look for (typically protein or
  ##                           peptide ID)
  ##            group:         element of 'ccList' in which we should
  ##                           look for 'elId'
  ##            ccList:  list of pre-processed connected component
  ##            verbose:       ???
  ## ----------------------------------------------------------------------
  ## Author: Sarah Gerster, Date:  13 Aug 2012, 17:43

  if (verbose > 1)
    message("    ** 'isInCC' returns TRUE/FALSE depending if the connected ",
            "component holds 'elID' in its element 'group' or not.")
  
  return(elID %in% unlist(ccList[[ccID]][group]))
}


getMyNeighborhoodSize <- function(elID, group="peptides", edgesPP,
                                  verbose=FALSE)
{
  ## Purpose: Return nighborhood size of elID
  ## ----------------------------------------------------------------------
  ## Arguments: elID:    index to look for (typically protein or peptide ID)
  ##            group:   type of element ("peptides" or "proteins") to
  ##                     which 'elID' belongs
  ##            edgesPP: data frame with 2 columns: 'pepId' and 'protId'
  ##                     -> each row defines an edge of the bipartite graph
  ##            verbose: if true, print detailed information
  ## ----------------------------------------------------------------------
  ## Author: Sarah Gerster, Date: 14 Aug 2012, 14:18

  if (verbose > 1)
    message("    ** 'getMyNeighborhoodSize' returns the number of ",
            "neighboring peptides or proteins to 'elID'")
  
  if (group == "peptides") {
    myId <- "pepId"
    neighId <- "protId"
  } else if (group == "proteins") {
    myId <- "protId"
    neighId <- "pepId"
  } else {
    stop("unknown 'group'")
  }
  
  ## return number of neighbors of element elID
  return(length(edgesPP[edgesPP[,myId]==elID, neighId]))
}

getProtId <- function(protName, proteins)
{
  ## Purpose: Return identification number assigend to protein 'protName'
  ## ----------------------------------------------------------------------
  ## Arguments: protName: protein accession
  ##            proteins: data frame holding all proteins in the bipartite
  ##                      graph.
  ##                      -> columns: 'protId', 'protName', ...
  ##                      -> one row per observed protein
  ## ----------------------------------------------------------------------
  ## Author: Sarah Gerster, Date:  5 Sep 2011, 09:13
  if (protName %in% proteins[,"protName"]) {
    return(proteins[proteins[,"protName"]==protName,"protId"])
  } else {
    return(NA)
  }
}

plotCC <- function(cc, peptides, proteins) {
  eds <- cc$edges
  pepseqs <- peptides[as.numeric(cc$peptides),"pepSeq"]
  protnames <- as.character(proteins[as.numeric(cc$proteins)-nrow(peptides),
                                     "protName"])
  nPeptides <- length(pepseqs)
  nProteins <- length(protnames)

  ## renumber peptides and proteins
  peps <- data.frame(pepId=as.numeric(cc$peptides), newId=1:nPeptides)
  prots <- data.frame(protId=as.numeric(cc$proteins),
                      newId=(1+nPeptides):(nPeptides+nProteins))
  eds[,"pepId"] <- sapply(eds[,"pepId"], replacePepNumbers, peps)
  eds[,"protId"] <- sapply(eds[,"protId"], replaceProtNumbers, prots)
  ccplot <- buildBipartiteGraph(nPeptides, nProteins, eds)
  nAttrs <- list()
  #nAttrs$label <- c(pepseqs,protnames)
  #nAttrs$label <- c(paste(pepseqs, cc$pepQty,sep="\n"), protnames)
  nAttrs$label <- c(round(cc$pepQty,2), protnames)
  names(nAttrs$label) <- 1:(nPeptides+nProteins)
  plot(ccplot, nodeAttrs=nAttrs, lwd=3)  
}

getTopNScore <- function(protID, peptides, edgesPP, n=3,
                         method="strict")
{
  ## Purpose: Compute the topN score (average over the n highest peptide
  ##          matching intensities) for the protein 'protID'.
  ## ----------------------------------------------------------------------
  ## Arguments: protID:   ID of protein for which the concentration should
  ##                      be estimated.
  ##            peptides: Data frame with the peptide information:
  ##                      -> columns: 'pepId', 'pepSeq', 'pepObs' &
  ##                         'pepQty'
  ##                      -> one row per observed peptide
  ##            edgesPP:  Data frame with 2 columns: 'pepId' and 'protId'.
  ##                      Each row defines an edge of the bipartite graph.
  ##            n:        Number of peptides used for the average.
  ##            method:   If strict, return NA when less than n peptides
  ##                      are available for a given protein. Else, return
  ##                      mean of whatever is available.
  ## ----------------------------------------------------------------------
  ## Dependencies: function getMyPeptides()
  ## ----------------------------------------------------------------------
  ## Author: Sarah Gerster, Date: 15 Mar 2011, 12:09

  ## get peptides sharing an edges with 'protID'
  pepIDs <- getMyPeptides(protID=protID, edgesPP=edgesPP)

  ## get intensities of the peptides
  pepQty <- peptides[peptides[,"pepId"] %in% pepIDs, "pepQty"]

  if (length(pepQty) >= n) {
    return(mean(sort(pepQty, decreasing=TRUE)[1:n]))
  } else {
    if (method=="strict") {
      ## not able to quantify proteins with less than n matching peptides
      return(NA)
    } else {
      ## make exceptions for proteins identified by less peptides and
      ## quantify them anyway by averaging whatever we have
      return(mean(pepQty))
    }
  }
}

getMyPeptides <- function(protID, edgesPP, group="proteins", ccList=NULL)
{
  ## Purpose: Return the list of all peptides with an edge to protein
  ##          'protID'
  ## ----------------------------------------------------------------------
  ## Arguments: protID:  ID of the protein for which the list of matching
  ##                     peptides is needed
  ##            edgesPP: Data frame with 2 columns: 'pepId' and 'protId'.
  ##                     Each row defines an edge of the bipartite graph.
  ## ----------------------------------------------------------------------
  ## Author: Sarah Gerster, Date: 31 Aug 2011, 14:18

  if (group=="proteins") {
    ## return list of peptides with an edge to protID
    return(edgesPP[edgesPP[,"protId"]==protID,"pepId"])
  } else if (group=="ccs") {
    return(length(ccList[[protID]]$peptides))
  } else {
    stop("Unknown 'group'.")
  }
}

estimateLSEparam <- function(peptides, ccList, verbose=FALSE)
{
  ## Purpose: Estimate the model paramteres alpha, beta, mu and tau from 
  ##          the input data with a least squared error approach on the
  ##          elements of the covariance matrix
  ## ----------------------------------------------------------------------
  ## Arguments: peptides:      dataframe with the peptide information:
  ##                           -> required columns: 'pepId', 'pepSeq',
  ##                              'pepObs' & 'pepQty'
  ##                           -> one row per observed peptide
  ##            ccList:  list of pre-processed connected components
  ##            verbose:       If true, print information during function
  ##                           evaluation
  ## ----------------------------------------------------------------------
  ## Dependencies: function getSampleCov()
  ##               function getBetaContribLSE()
  ##               function getTauContribLSE()
  ##               global constant TAU2_MIN_CONST
  ##               global constant BETA22_MIN_CONST
  ## ----------------------------------------------------------------------
  ## Author: Sarah Gerster, Date: 14 Aug 2011, 15:00

  if (verbose)
    message("  ** compute parameter estimates with LSE approach")

  if((all(peptides[,"nrNeiProt"] == 1) &&
      all(peptides[,"pepObs"] == mean(peptides[,"pepObs"])))) {
    muHat <- 0
    
    ## mean of all peptide intensities
    alphaHat <- mean(peptides[,"pepQty"])
    meanVect <- rep(alphaHat, nrow(peptides))
  } else {
    ## mean of all peptide intensities
    if (verbose) message("    * Mean vector")
    y.pepqty <- peptides[,"pepQty"]
    x.pepobsnrnei <- peptides[,"pepObs"]*peptides[,"nrNeiProt"]
    fit <- lm(y.pepqty ~ x.pepobsnrnei)
    alphaHat <- fit$coefficients[1]     # alphaHat is the intercept
                                        # of the fit
    betaMuHat <- fit$coefficients[2]    # slope of the fit is the
                                        # product of betaHat and
                                        # muHat
    meanVect <- alphaHat +
      peptides[,"pepObs"]*betaMuHat*peptides[,"nrNeiProt"]
  }
  
  ## estimate the covariance matrix for each CC from the data
  if (verbose) message("    * Estimate sample covariance matrices")
  ccList <- lapply(ccList, getSampleCov, myMean=meanVect)
  
  ## estimate beta from the off-diagonal elements of cov(U)
  if (verbose) message("    * Compute beta2")
  t.res <- as.data.frame(lapply(ccList, getBetaContribLSE))
  betaHat2 <- max(sum(t.res[1,])/sum(t.res[2,]), BETA2_MIN_CONST)
    
  ## tau out of diagonal elements of cov(U)
  if (verbose) message("    * Compute tau2")
  t.res <- sum(unlist(lapply(ccList, getTauContribLSE, betaHat2)))
    
  ## tau2 is maximum between estimated value (might be a negative number...)
  ## and a preset minimal (positive) value
  tauHat2 <- max(t.res/length(ccList), TAU2_MIN_CONST)

  
  if(!(all(peptides[,"nrNeiProt"] == 1) &&
      all(peptides[,"pepObs"] == mean(peptides[,"pepObs"])))) {
    ## mu is approximated with the slope of di vs Ui and betaHat
    if (verbose) message("    * Compute mu")
    muHat <- betaMuHat/sqrt(betaHat2)
  }
  
  if (verbose) cat("  => run finished: alpha =", alphaHat,"beta =",
                   sqrt(betaHat2), "mu = ", muHat, "and tau =",
                   sqrt(tauHat2), "\n")
  
  ## return parameter estimates
  res <- c(alphaHat, sqrt(betaHat2), muHat, sqrt(tauHat2))
  names(res) <- c("alphaH", "betaH", "muH", "tauH")
  
  if (verbose)
    message("  --> LSE parameter estimates computed successfully")

  return(res)
}

getSampleCov <- function(cc, myMean, verbose=FALSE)
{
  ## Purpose: Compute the sample covariance of the peptide abundances
  ##          of all peptides in 'cc'
  ## ----------------------------------------------------------------------
  ## Arguments: cc:     Pre-processed connected component: the sample
  ##                    covariance matrix should be added to the list
  ##            myMean: Mean value of the peptide intensities
  ## ----------------------------------------------------------------------
  ## Author: Sarah Gerster, Date: 31 Aug 2011, 15:20

  if (verbose>1)
    message("  ** compute sample covariance matrix for the peptides in 'cc'")
  
  ## sample covariance matrix
  cc$covU <- (cc$pepQty - myMean[as.numeric(cc$peptides)]) %*%
    t(cc$pepQty - myMean[as.numeric(cc$peptides)])
  
  return(cc)
}

getBetaContribLSE <- function(cc)
{
  ## Purpose: compute the contribution of the connected component 'cc'
  ##          towards the estimate of beta
  ## ----------------------------------------------------------------------
  ## Arguments: cc:          pre-processed connected component: the sample
  ##                         covariance matrix should be added to the list
  ## ----------------------------------------------------------------------
  ## Author: Sarah Gerster, Date:  1 Sep 2011, 10:35

  resNum <- 0
  resDen <- 0

  if (length(cc$peptides)>1) {
    ## NOTE: only CCs with more than one peptide contribute!
    for (i in 1:(ncol(cc$covU)-1)) {
      for (j in (i+1):ncol(cc$covU)) {
        ## NOTE: Only summing over the upper triangular elements,
        ##       since 'covU' and 'distance' are symmetric
        xij <- cc$pepObs[i]*cc$pepObs[j]*cc$distance[i,j]
        resNum <- resNum + cc$covU[i,j]*xij
        resDen <- resDen + xij^2
      }
    }
  } 

  return(c(resNum,resDen))
}

getTauContribLSE <- function(cc, beta2)
{
  ## Purpose: compute the contribution of the connected component 'cc'
  ##          towards the estimate of tau.
  ## ----------------------------------------------------------------------
  ## Arguments: cc:          pre-processed connected component: the sample
  ##                         covariance matrix should be added to the list
  ##            beta2:       value of the parameter beta^2
  ## ----------------------------------------------------------------------
  ## Author: Sarah Gerster, Date:  1 Sep 2011, 11:33

  resZi <- 0
  for (i in 1:ncol(cc$covU)) {
    ## NOTE: all CCs contribute!
    if (ncol(cc$covU) > 1) {
      resZi <- resZi + cc$covU[i,i]-(cc$pepObs[i])^2*beta2*cc$distance[i,i]
    } else {
      resZi <- resZi + cc$covU[i,i]-(cc$pepObs[i])^2*beta2*cc$distance
    }
  }
  
  return(resZi/length(cc$peptides))
}


estimateMLEparam <- function(paramStart, ccList, verbose=FALSE)
{
  ## Purpose: Compute estimates of the model paramteres alpha, beta, mu 
  ##          and tau with the MLE, assuming that the peptide intensities
  ##          follow a multivariate normal distribution.
  ## ----------------------------------------------------------------------
  ## Arguments: paramStart:  starting values used for alphaHat, betaHat,
  ##                         muHat and tauHat
  ##            ccList:      list of pre-processed connected components
  ##            verbose:     if true, print information during function
  ##                         evaluation
  ## ----------------------------------------------------------------------
  ## Dependencies: function nlLtoMinimze()
  ## ----------------------------------------------------------------------
  ## Author: Sarah Gerster, Date: 14 Aug 2012, 15:00

  if (verbose)
    message("  ** compute parameter estimates with MLE approach")
  
  if (verbose) {
    cat("     - start parameters:\n")
    cat("       alphaHs =", paramStart[1], ", betaHs =", paramStart[2],
        ", muHs =", paramStart[3], "and tauHs =", paramStart[4], "\n")
  }
  
  res <- optim(paramStart, fn=nlLtoMinimize, method = "L-BFGS-B",
               lower = c(-Inf, 0, 0, 0),
               #lower = c(-Inf, -Inf, -Inf, 0),
               upper = c(Inf, Inf, Inf, Inf),
               control = list(factr=1e12, trace=1),
               ccList = ccList, verbose=verbose)
  out <- c(res$par[1], res$par[2], res$par[3], res$par[4], res$value)
  names(out) <- c("alphaH","betaH","muH","tauH", "fnval")
 
  if (verbose) {
    cat("  => run finished: alpha =", res$par[1], "beta =",
        res$par[2], "mu =", res$par[3], "and tau =", res$par[4], "\n")
  }

  if (verbose)
    message("  --> MLE parameter estimates computed successfully")
  
  return(out)
}

nlLtoMinimize <- function(paramAlphaBetaMuTau, ccList,
                          verbose=FALSE)
{
  ## Purpose: sum up the negative log-likelihood (nlL) of all connected
  ##          components
  ## ----------------------------------------------------------------------
  ## Arguments: paramAlphaBetaMuTau: Values used for alphaHat, betaHat,
  ##                                 muHat and tauHat
  ##            ccList:   list of pre-processed connected component
  ##            verbose:        if true, print information during function
  ##                            evaluation
  ## ----------------------------------------------------------------------
  ## Dependencies: function getCCnlL()
  ##               function getCovU()
  ## ----------------------------------------------------------------------
  ## Author: Sarah Gerster, Date: 14 Aug 2012, 15:00

  ## complete CC preprocessing with current parameter values
  ## --> covU
  ccList <- lapply(ccList, getCovU,
                         beta=paramAlphaBetaMuTau[2],
                         tau=paramAlphaBetaMuTau[4])
  
  ccVals <- unlist(lapply(ccList, getCCnlL, paramAlphaBetaMuTau))
  res <- sum(ccVals)

  if (verbose > 1) {
    cat("       Call to the function nlLtoMinimize() just finished:\n" )
    cat("       Negative log-likelihood =", res, "\n")
    cat("       Parameters: alphaH =", paramAlphaBetaMuTau[1],", betaH =",
        paramAlphaBetaMuTau[2],", muH =", paramAlphaBetaMuTau[3], "and tauH =",
        paramAlphaBetaMuTau[4], "\n\n")
  }

  return(res)
}


getCCnlL <- function(cc, paramAlphaBetaMuTau)
{
  ## Purpose: compute the contribution to the negative log-likelihood (nlL)
  ##          of the connected component 'cc'
  ## ----------------------------------------------------------------------
  ## Arguments: cc:                  pre-processed connected component 
  ##            paramAlphaBetaMuTau: values used for alphaHat, betaHat,
  ##                                 muHat and tauHat
  ## ----------------------------------------------------------------------
  ## Dependencies: library(mvtnorm)
  ##               global constant DENSITY_CONST
  ## ----------------------------------------------------------------------
  ## Author: Sarah Gerster, Date: 15 Aug 2012, 15:00
  
  ## compute the value of the probability density function
  if (length(cc$pepQty) > 1) {
      return(-log(max(dmvnorm(as.vector(cc$pepQty),
                              mean=paramAlphaBetaMuTau[1] + cc$pepObs*
                              paramAlphaBetaMuTau[2]*paramAlphaBetaMuTau[3]*
                              diag(cc$distance),
                              sigma=cc$covU),
                      DENSITY_CONST)))
  } else {
    return(-log(max(dnorm(cc$pepQty,
                          mean=paramAlphaBetaMuTau[1]+cc$pepObs*
                          paramAlphaBetaMuTau[2]*paramAlphaBetaMuTau[3]*
                          cc$distance,
                          sd=cc$covU),
                    DENSITY_CONST)))
  }  
}


getCovU <- function(cc, beta, tau)
{
  ## Purpose: compute the covariance matrix of the peptide intensities for
  ##          a connected component for given parameter values
  ## ----------------------------------------------------------------------
  ## Arguments: cc:          pre-processed connected component
  ##            beta:        parameter in modelChoice
  ##            tau:         parameter in modelChoice
  ## ----------------------------------------------------------------------
  ## Author: Sarah Gerster, Date: 2 Sep 2011, 15:00
  
  ## different handling depending if the cc holds 1 or more peptides
  if (length(cc$pepObs) > 1) {
    covU <- tcrossprod(cc$pepObs) * cc$distance * beta^2
    diag(covU) <- diag(covU) + tau^2
    cc$covU <- covU
    return(cc)
  } else {
    cc$covU <- (cc$pepObs*beta)^2*cc$distance + tau^2
    return(cc)
  }
}

sampleInitialVector <- function(alphaMin=0, alphaMax=4,
                                betaMin=0, betaMax=4,
                                muMin=0, muMax=4,
                                tauMin=0, tauMax=4)
{
  ## Purpose: generate a random initial vector for alpha, beta, mu and tau
  ## ----------------------------------------------------------------------
  ## Arguments: *Min, *Max: min and max values between which each starting
  ##                        parameter should be sampled
  ## ----------------------------------------------------------------------
  ## Author: Sarah Gerster, Date: 03 May 2012, 08:30

  return(c(runif(1,alphaMin, alphaMax), runif(1,betaMin, betaMax),
           runif(1,muMin, muMax), runif(1,tauMin, tauMax)))
}

replacePepNumbers <- function(edPepID, pepList){
  ## Purpose: return new id for peptide with old id 'edPepID'
  ## ----------------------------------------------------------------------
  ## Arguments: edPepID: old peptide ID
  ##            pepList: dataframe with peptides, must include columns
  ##                     named 'pepId' and 'newId'
  ## ----------------------------------------------------------------------
  ## Author: Sarah Gerster, Date: 16 Aug 2012, 16:43
  return(pepList[which(pepList[,"pepId"]==edPepID),"newId"])
}

replaceProtNumbers <- function(edProtID, protList){
  ## Purpose: return new id for protein with old id 'edProtID'
  ## ----------------------------------------------------------------------
  ## Arguments: edPepID: old protein ID
  ##            pepList: dataframe with proteins, must include columns
  ##                     named 'protId' and 'newId'
  ## ----------------------------------------------------------------------
  ## Author: Sarah Gerster, Date: 16 Aug 2012, 16:43
  return(protList[which(protList[,"protId"]==edProtID),"newId"])
}

plot.scampi <- function(x, ...) {
  ## histograms of measured peptide abundnace scores
  hist(x@peptides[,"pepQty"],
       freq=TRUE,
       xlab="measured peptide abundance",
       ylab="frequency",
       main="measured peptide abundances",
       ...)
}

plot.scampiVal <- function(x, ...) {
  plotRows <- length(names(x@parameters))
  if (all(paste(names(x@parameters), "Score", sep="") %in%
          colnames(x@peptides))) {
    plotCols <- 3
  } else {
    plotCols <- 1
  }
  op <- par(mfrow=c(plotRows, plotCols))
  
  for (methIter in names(x@parameters)){
    ## histograms of protein abundnace scores
    hist(x@proteins[,paste(methIter, "Score", sep="")],
         freq=TRUE,
         xlab="protein abundance score",
         ylab="frequency",
         main=paste("protein abundance SCAMPI(", methIter, ")",
           sep=""), ...)
    
    if (paste(methIter, "Score", sep="") %in%
        colnames(x@peptides)) {
      resid <- unlist(x@peptides[paste(methIter, "Resid", sep="")])
      fitted <- unlist(x@peptides[paste(methIter, "Score", sep="")])

      ## Tukey-Anscombe plot
      plot(fitted, resid, 
           main = paste("peptide reassessment - SCAMPI(",
             methIter, ")", sep=""),
           sub = "Tukey-Anscombe plot",
           xlab = "fitted",
           ylab = "residuals", ...)
      abline(h = 0, lwd = 2)
      
      ## normal Q-Q plot
      qqnorm(y = resid,
             main = paste("peptide reassessment - SCAMPI(",
               methIter, ")", sep=""),
             sub = "Normal Q-Q Plot",
             bty = "n",
             xlab = "theoretical quantiles",
             ylab = "residuals", ...)
      qqline(resid, lwd = 2)
    } 
  }
  par(op)
}

summary.scampi <- function(object, ...) {
  cat("\nObject of class 'scampi'.\n")
  cat("Data overview:\n",
      "==============\n")
  cat("Number of proteins       :", nrow(object@proteins), "\n")
  cat("Number of peptides       :", nrow(object@peptides), "\n")
  cat("Number of pep-prot edges :", nrow(object@edgespp), "\n")
}

summary.scampiVal <- function(object, ...) {
  cat("\nObject of class 'scampiVal', from Call: \n",
      deparse(object@call), "\n")
  cat("Input data:\n",
      "===========\n")
  cat("Number of proteins       :", nrow(object@proteins), "\n")
  cat("Number of peptides       :", nrow(object@peptides), "\n")
  cat("Number of pep-prot edges :", nrow(object@edgespp), "\n")
  cat("Percent. shared peptides :",
      sum(object@peptides[,"nrNeiProt"]>1)/nrow(object@peptides),
      "%\n")
  cat("Output data:\n",
      "===========\n")
  for (methIter in names(object@parameters)){
    cat("**", methIter, ":\n")
    cat("   protein abundance scores:\n")
    print(summary(object@proteins[,paste(methIter, "Score",
                                         sep="")]))
    if (paste(methIter, "Score", sep="") %in%
        colnames(object@peptides)) {
      cat("   peptide abundance scores:\n")
      print(summary(object@peptides[,paste(methIter, "Score",
                                           sep="")]))
      cat("   peptide residuals:\n")
      print(summary(object@peptides[,paste(methIter, "Resid",
                                           sep="")]))
    }
  }
}

show.scampi <- function(object) {
  cat("Object of class 'scampi'.\n")
  invisible(object)
}

show.scampiVal <- function(object) {
  cat("Object of class 'scampiVal', from Call: \n",
      deparse(object@call),"\n")
  invisible(object)
}


## ####################### ##
## SET METHODS FOR CLASSES ##
## ####################### ##
#setOldClass("scampiVal")
setMethod("show", "scampiVal", show.scampiVal)

setMethod("summary", "scampiVal", summary.scampiVal)

setMethod("plot", signature(x = "scampiVal", y = "missing"),
          function(x, y, ...) plot.scampiVal(x,...))

#setOldClass("scampi")
setMethod("show", "scampi", show.scampi)

setMethod("summary", "scampi", summary.scampi)

setMethod("plot", signature(x = "scampi", y = "missing"),
          function(x, y, ...) plot.scampi(x,...))
