Estimation <- function(Data, method, iters = 10, breaking.method = NA, prior = 0, race = FALSE) {
  m <- ncol(Data)
  
  if(!is.na(breaking.method))
    Data.pairs <- Breaking(Data, method = breaking.method)
  
  if(suppressWarnings(!is.na(as.numeric(substring(method, 1, 1)))))
    K <- as.numeric(substring(method, 1, 1))
  
  if(method == "PL.GMM") Estimation.PL.GMM(Data.pairs, m, prior = prior)
  
  else if(method == "PL.MLE")   {
    if(breaking.method == "full") Estimation.PL.MLE(Data, iters)
    else Estimation.RUM.MLE(Data, iters, dist = "exp", race = race)
  }
  
  else if(method == "N.FV.GMM") Estimation.Normal.GMM(Data.pairs, m, iter = iters, prior = prior)
  else if(method == "N.DV.GMM") Estimation.Normal.GMM(Data.pairs, m, iter = iters, Var = TRUE, prior = prior)
  else if(method == "N.FV.MLE") Estimation.RUM.MLE(Data, iters, dist = "norm.fixedvariance", race = race)
  else if(method == "N.DV.MLE") Estimation.RUM.MLE(Data, iters, dist = "norm", race = race)
  
  else if(substring(method, 2) == "PL.MLE")   Estimation.RUM.MultiType.MLE(Data, K, iters, dist = "exp",  race = race)
  else if(substring(method, 2) == "N.FV.MLE") Estimation.RUM.MultiType.MLE(Data, K, iters, dist = "norm", race = race)
  
  else if(method == "Zemel") Estimation.Zemel.MLE(Data.pairs, m)
  else if(method == "Nonparametric") Estimation.RUM.Nonparametric(Data, m, iter = iters, race = race)
  
  else stop(paste("Method", method, "not recognized"))
}