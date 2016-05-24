elmtrain.default <-
function(x,y,nhid,actfun,...) {
  require(MASS)
  
  if(nhid < 1) stop("ERROR: number of hidden neurons must be >= 1")
  
  T <- t(y)
  P <- t(x)
  
  inpweight <- randomMatrix(nrow(P),nhid,-1,1)
  tempH <- inpweight %*% P
  biashid <- runif(nhid,min=-1,max=1)
  biasMatrix <- matrix(rep(biashid, ncol(P)), nrow=nhid, ncol=ncol(P), byrow = F) 
  
  tempH = tempH + biasMatrix
  
  if(actfun == "sig") H = 1 / (1 + exp(-1*tempH))
  else {
    if(actfun == "sin") H = sin(tempH)
    else {
      if(actfun == "radbas") H = exp(-1*(tempH^2))
      else {
        if(actfun == "hardlim") H = hardlim(tempH)
        else {
          if(actfun == "hardlims") H = hardlims(tempH)
          else {
            if(actfun == "satlins") H = satlins(tempH)
            else {
              if(actfun == "tansig") H = 2/(1+exp(-2*tempH))-1
              else {
                if(actfun == "tribas") H = tribas(tempH)
                else {
                  if(actfun == "poslin") H = poslin(tempH)
                  else {
                    if(actfun == "purelin") H = tempH
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
  
  outweight <- ginv(t(H), tol = sqrt(.Machine$double.eps)) %*% t(T)
  Y <- t(t(H) %*% outweight)
  model = list(inpweight=inpweight,biashid=biashid,outweight=outweight,actfun=actfun,nhid=nhid,predictions=t(Y))
  model$fitted.values <- t(Y)
  model$residuals <- y - model$fitted.values
  model$call <- match.call()
  class(model) <- "elmNN"
  model
}
