predict.elmNN <-
function(object,newdata=NULL,...) {
  if(is.null(newdata))
    predictions <- fitted(object)
  else {
    if(!is.null(object$formula)){
      ## model has been fitted using formula interface
      x <- model.matrix(object$formula, newdata)
    }
    else{
      x <- newdata
    }
  
    inpweight <- object$inpweight
    biashid <- object$biashid
    outweight <- object$outweight
    actfun <- object$actfun
    nhid <- object$nhid
  
    TV.P <- t(x)
  
    tmpHTest=inpweight %*% TV.P
  
    biasMatrixTE <- matrix(rep(biashid, ncol(TV.P)), nrow=nhid, ncol=ncol(TV.P), byrow = F)
    
    tmpHTest = tmpHTest + biasMatrixTE
  
    if(actfun == "sig") HTest = 1 / (1 + exp(-1*tmpHTest))
    else {
      if(actfun == "sin") HTest = sin(tmpHTest)
      else {
        if(actfun == "radbas") HTest = exp(-1*(tmpHTest^2))
        else {
          if(actfun == "hardlim") HTest = hardlim(tmpHTest)
          else {
            if(actfun == "hardlims") HTest = hardlims(tmpHTest)
            else {
              if(actfun == "satlins") HTest = satlins(tmpHTest)
              else {
                if(actfun == "tansig") HTest = 2/(1+exp(-2*tmpHTest))-1
                else {
                  if(actfun == "tribas") HTest = tribas(tmpHTest)
                  else {
                    if(actfun == "poslin") HTest = poslin(tmpHTest)
                    else {
                      if(actfun == "purelin") HTest = tmpHTest
                      else stop(paste("ERROR: ",actfun," is not a valid activation function.",sep=""))
                    }
                  }
                }
              }
            }
          }
        }
      }
    }
    TY = t(t(HTest) %*% outweight)
  
    predictions <- t(TY)
    }
  predictions
}
