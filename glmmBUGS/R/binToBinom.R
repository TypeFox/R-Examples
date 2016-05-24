binToBinom = function(obs, covariates) {
# make sure covariates is a matrix
if(is.vector(covariates))
  covariates = matrix(covariates, ncol=1)

# get rid of missings
notNA = !is.na(obs)
obs = obs[notNA]
covariates = covariates[notNA,]


  covString = apply(covariates, 1, toString)
  tableCov = c(table(covString))
  
  multCov = names(tableCov[tableCov>1])
  singleCov = names(tableCov[tableCov==1])
  
  singleCov = which(covString %in% singleCov)
  result = cbind(covariates[singleCov,], y=obs[singleCov], N=rep(1, length(singleCov)))



  if(length(multCov)) {
  whichMultCov = which(covString %in% multCov)
  covString = covString[whichMultCov]
  theorder = order(covString)
  covString = covString[theorder]
  whichMultCov = whichMultCov[theorder]
  obs = obs[whichMultCov]
  covariates = covariates[whichMultCov,]

  thediff = which(covString[-1] != covString[-length(covString)])
  thediff = c(thediff, length(covString))
  
  Nseq = (1:length(covString))[thediff]
  Nseq = diff(c(0, Nseq))
  Yseq = cumsum(obs)[thediff]
  Yseq = diff(c(0, Yseq))
  
  result = rbind(result, cbind(covariates[thediff,], y=Yseq, N=Nseq))
  

  }
  return(result)  
}