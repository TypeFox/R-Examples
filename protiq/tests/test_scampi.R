library("protiq")
data(leptoSRM)

## reading data and checking input
dataTest <- new("scampi", peptides=leptoSRMpeptides,
                proteins=leptoSRMproteins, edgespp=leptoSRMedgespp)
if (nrow(dataTest@peptides) != 151 ||
    nrow(dataTest@proteins) != 39 ||
    nrow(dataTest@edgespp) != 151) {
  stop("Error by reading input data or building 'scampi' object.")
}

dataChecked <-
  checkInputData(new("scampi", peptides=leptoSRMpeptides,
                     proteins=leptoSRMproteins, edgespp=leptoSRMedgespp),
                     rescaling=FALSE)
if (nrow(dataChecked@peptides) != 151 ||
    nrow(dataChecked@proteins) != 39 ||
    nrow(dataChecked@edgespp) != 151) {
  stop("Error by input data checking.")
}

## building the graph structure, extracting information from graph
ppGraph <- protiq:::buildBipartiteGraph(nrow(dataChecked@peptides),
                                        nrow(dataChecked@proteins),
                                        dataChecked@edgespp)

protiq:::getMyProteins(1, dataChecked@edgespp)
protiq:::getDmat(c(1,2,3,7,8,24), dataChecked@edgespp)
protiq:::getMyNeighborhoodSize(8, group="peptides", dataChecked@edgespp)

## decompose graph in connected components and work with these
ppComponents <- RBGL:::connectedComp(ppGraph)
testCC <-
  protiq:::preprocessCC(ppComponents[[1]], dataChecked@peptides,
                        dataChecked@edgespp, nrow(dataChecked@peptides))

ppComponents <-
  lapply(ppComponents, protiq:::preprocessCC, dataChecked@peptides,
         dataChecked@edgespp, nrow(dataChecked@peptides))

protiq:::isInCC(3, 21, group="peptides", ppComponents)
protiq:::getMyCCNr(24, "peptides", ppComponents)

## preprocessing the input data
tmpPrepro <- preprocessInputData(dataChecked)
dataPrepro <- tmpPrepro[["dataPrepro"]]
ppGraph <- tmpPrepro[["ppGraph"]]
ppComponents <- tmpPrepro[["ccList"]]
rm("tmpPrepro")

## class for output
dataOutTest <-
  new("scampiVal", peptides=dataChecked@peptides,
      proteins=dataChecked@proteins,
      edgespp=dataChecked@edgespp, parameters=list(),
      ppGraph=ppGraph, ccList=ppComponents)

## parameter estimation
protiq:::getSampleCov(ppComponents[[4]],
                      myMean=rep(3,nrow(dataPrepro@peptides)))
lseparam <- protiq:::estimateLSEparam(dataPrepro@peptides, ppComponents)
getCovU(ppComponents[[4]], beta=1, tau=0.2)
protiq:::getCCnlL(getCovU(ppComponents[[4]], beta=1, tau=0.2),
                  c(0.5,1,2,0.2))
scampiParam <- estimateModelParameters(method="LSE", ppComponents,
                                       peptides=dataPrepro@peptides,
                                       numIter=3)
estimateModelParameters(method="LSE", ppComponents,
                        peptides=dataPrepro@peptides)

## quantify proteins and reassess peptides
ppComponents <- lapply(ppComponents, getCovU,
                       beta=scampiParam[["LSE"]]["betaH"],
                       tau=scampiParam[["LSE"]]["tauH"])
quantifyProtein(dataPrepro@proteins[1,c("protId","ccInd")], ppComponents,
                scampiParam[["LSE"]])
quantifyPeptide(dataPrepro@peptides[1,c("pepId","ccInd")], ppComponents,
                scampiParam[["LSE"]])
scampiOut <- quantifyProteins(dataPrepro, ppComponents, scampiParam,
                              quantifyPeptides=TRUE)

## test of the two main functions
showProc.time <- local({
    pct <- proc.time()
    function() { ## CPU elapsed __since last called__
        ot <- pct ; pct <<- proc.time()
        cat('Time elapsed: ', (pct - ot)[1:3],'\n')
    }
})

scampiTest <-
  runScampi(leptoSRMpeptides, leptoSRMproteins, leptoSRMedgespp,
            rescaling=FALSE, method="LSE", quantifyPeptides=TRUE,
            numIter=3, verbose=FALSE)
showProc.time()

scampiIterateTest <-
  iterateScampi(leptoSRMpeptides, leptoSRMproteins, leptoSRMedgespp,
                rescaling=FALSE, method="LSE", numIter=3, thresh=1.73,
                verbose=FALSE)
showProc.time()


## commented for time reasons (keep test brief):				       
#mleparam <- protiq:::estimateMLEparam(c(0.1,0.31,3,0.1), ppComponents)
#scampiParam <- estimateModelParameters(method="all", ppComponents,
#                                       peptides=dataPrepro@peptides,
#                                       numIter=3)
#estimateModelParameters(method="MLE", ppComponents,numIter=3)
#scampiTest <-
#  runScampi(leptoSRMpeptides, leptoSRMproteins, leptoSRMedgespp,
#            rescaling=FALSE, method="all", quantifyPeptides=TRUE,
#            numIter=3, verbose=FALSE)
#showProc.time()
#scampiIterateTest <-
#  iterateScampi(leptoSRMpeptides, leptoSRMproteins, leptoSRMedgespp,
#                rescaling=FALSE, method="MLE", numIter=3,
#                numMLEIter=3, thresh=1.73, verbose=FALSE)
#showProc.time()
