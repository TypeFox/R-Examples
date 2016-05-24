getBestModel <-
function(mogavs,nvar,method=c("AIC","BIC","knee","mse",NULL)){
  if(class(mogavs)!="mogavs") stop("Arg 1 is not of class 'mogavs'.")
  bestModel<-NULL
  if(missing(method) && missing(nvar)) stop("You have to supply nvar or method")
  if(missing(method)) method<-"mse"
  if(!is.null(method)){
    if(method=="AIC"){
      #return best model according to AIC
      z<-sapply(1:nrow(mogavs$nonDominatedSet),function(x) mogavs$n_obs*log(mogavs$MSE[x])+2*mogavs$numOfVariables[x])
      ind<-which(z==min(z))
    }
    else if(method=="BIC"){
      #return best method according to BIC
      z<-sapply(1:nrow(mogavs$nonDominatedSet),function(x) mogavs$n_obs*log(mogavs$MSE[x])+log(mogavs$n_obs)*mogavs$numOfVariables[x])
      ind<-which(z==min(z))
    }
    else if(method=="knee"){
      #return best method according to knee point
      x1<-min(mogavs$numOfVariables)
      ind1<-which(mogavs$numOfVariables==x1)
      y1<-mogavs$MSE[ind1]
      x2<-max(mogavs$numOfVariables)
      ind2<-which(mogavs$numOfVariables==x2)
      y2<-mogavs$MSE[ind2]
      
      dmax<-0
      for(i in 2:(nrow(mogavs$nonDominatedSet)-1)){
        x0<-mogavs$numOfVariables[i]
        y0<-mogavs$MSE[i]
        d<-abs((y2-y1)*x0-(x2-x1)*y0+x2*y1-y2*x1)/sqrt((y2-y1)^2+(x2-x1)^2)
        if(d>dmax){
          dmax<-d
          ind<-i
        }
      }
    }
    else if(method=="mse"){
      ind<-which(mogavs$numOfVariables==nvar)
    }
  }
  else {
    #return best model with nvar variables
  ind<-which(mogavs$numOfVariables==nvar)
  }
  bestModel<-mogavs$nonDominatedSet[ind,]
  return(bestModel)
}
