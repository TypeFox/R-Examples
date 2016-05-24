CVparallel <-
function(Dat, l1range, alpharange, K=5, nodes=2){
  # run CV in parallel!!
  # Choose regularisation parameters based on l1 and l2 based on CV error
  # Supplied vectors, l1range and l2range, should be quite short as this will be computationally expensive!
  #
  #
  # INPUT:
  #      - Dat: data used to fit MNS model
  #      - l1range: vector of potential l1 parameters
  #      - alpharange: vector of alpha parameters 
  #      - K: number of CV folds to perform
  #
  #
  
  n = nrow(Dat[[1]]) # number of rows
  p = ncol(Dat[[1]]) # number of ROIs
  N = length(Dat) # number of subjects
  CVseqs = floor(seq(1, n, length.out = K+1))
  CVseqs[1] = 0
  
  # build data list to give to parLapply:
  cvList = vector("list", K)
  for (i in 1:K){
    cvList[[i]] = vector("list", 2)
    testID = seq(CVseqs[i]+1, CVseqs[i+1])
    # first entry of list is the training data
    cvList[[i]][[1]] = lapply(Dat, FUN=function(subDat){
      subDat[-testID, ]
    })
    # second entry is the test data:
    cvList[[i]][[2]] = lapply(Dat, FUN=function(subDat){
      subDat[testID,]
    })
  }
  
  # no we iterate through each entry in cvList using parLapply:
  workers=makeCluster(2)
  registerDoParallel(workers,cores=nodes)
  
  #manually export things I will need
  clusterExport(cl=workers, varlist =list("l1range", "alpharange",  "MNS", "PrepareData", "PenalizedLMM_opt", "ginv", "BuildRandomDesign", "glmnet", "calculateNLL", "dmvnorm"),envir=environment())
  
  results = parLapply(cl = workers, X = cvList, fun = function(CVentry){
    trainDat = CVentry[[1]]
    testDat = CVentry[[2]]
    N = length(CVentry[[1]])
    
    CVerrorNode = matrix(0, nrow=length(alpharange), ncol=length(l1range))
    
    # iterate through reg parameter grid:
    for (l1 in 1:length(l1range)){
      for (alpha in 1:length(alpharange)){
        X = MNS(dat = trainDat, lambda_pop = l1range[l1]*alpharange[alpha], lambda_random = l1range[l1]*(1-alpharange[alpha]))
        
        # now we calculate CV error for each subject:
        for (sub in 1:N){
          Subresiduals = (testDat[[sub]] - t((X$PresPop + X$PresBLUP[,,1])%*%t(testDat[[sub]])))
          CVerrorNode[alpha, l1] = CVerrorNode[alpha, l1] + sum(Subresiduals**2)          
        }
      }
    }
    
    return(CVerrorNode)
  })
  
  stopCluster(workers)
  
  # now tidy up results and select reg parameters!
  #
  CVerror = t(Reduce('+', results))
  
  return(list(l1 = l1range[which(CVerror==min(CVerror), arr.ind = T)[,1][1]], alpha=alpharange[which(CVerror==min(CVerror), arr.ind = T)[,2][1]], CV=CVerror))
}
