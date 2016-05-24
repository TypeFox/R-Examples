plsr = function (Xtrain, Ytrain, Xtest, ncomp = NULL, type = c("simpls", "nipals"), 
                 unit.weights = TRUE, weight = FALSE, beta = 0.1) 
{
  if(type == "simpls"){
    if (unit.weights == TRUE) {
        return(unitsimpls(Xtrain, Ytrain, Xtest, ncomp, weight = weight, beta = beta))
    }
    else {
        return(simpls(Xtrain, Ytrain, Xtest, ncomp, weight = weight, beta = beta))
    }
  }
  else {
        return(nipals(Xtrain, Ytrain, Xtest, ncomp, weight = weight, beta = beta))
  }
}
