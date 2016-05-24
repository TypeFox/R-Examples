CVsequential <-
function(Dat, l1range, alpharange, K=5, verbose=FALSE){
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
  
  SubjectCVs = array(0, c(K, N, p))#matrix(0, ncol=p, nrow=N)
  CVerror = matrix(0, nrow=length(l1range), ncol=length(alpharange))
  
  # have to go through each pair - computationally very expensive!!
  for (l1 in 1:length(l1range)){
    for (alpha in 1:length(alpharange)){
      if (verbose) cat("Running CV for parameters l1:", l1range[l1]*alpharange[alpha], "and l2:", l1range[l1]*(1-alpharange[alpha]), "\n")
      # for each subject we cross-validate errors:
      for (k in 1:K){
        testID = seq(CVseqs[k]+1, CVseqs[k+1])
        # define training and test datasets:
        trainDat = lapply(Dat, FUN=function(subDat){
          subDat[-testID, ]
        })
        testDat = lapply(Dat, FUN=function(subDat){
          subDat[testID, ]
        })
        
        X = MNS(dat = trainDat, lambda_pop = l1range[l1]*alpharange[alpha], lambda_random = l1range[l1]*(1-alpharange[alpha]))
        
        # now we calculate CV error for each subject:
        for (sub in 1:N){
          #SubjectCVs[sub,] = apply(testDat[[sub]] %*% (X$PresPop + X$PresBLUP[,,1]), 2, FUN=function(x){ sum(x**2)})
          Subresiduals = (testDat[[sub]] - t((X$PresPop + X$PresBLUP[,,1])%*%t(testDat[[sub]])))
          SubjectCVs[k,,][sub,] = apply(Subresiduals, 2, FUN=function(x){ sum(x**2)})
          
        }
        CVerror[l1, alpha] = sum(SubjectCVs)
        if (verbose) cat("Fold", k, "complete\n")
      }
    }
  }
  return(list(l1 = l1range[which(CVerror==min(CVerror), arr.ind = T)[,1][1]], alpha=alpharange[which(CVerror==min(CVerror), arr.ind = T)[,2][1]], CV=CVerror))
}
