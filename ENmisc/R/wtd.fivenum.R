wtd.fivenum <-
function(x, weights=NULL, na.rm=TRUE)
{
    interpolatedindex<-function(myval,weights){
    indices<-1:length(weights)
    n<-sum(weights)
    weightsbelow<-rep(0,length(weights))
  for (i in 2:length(weights))
        weightsbelow[i] <- weightsbelow[i-1]+weights[i-1]
    weightsabove<-n-weightsbelow-weights
    lowcands<-weightsbelow<myval
    highcands<-weightsabove<n-myval
    (ifelse(any(lowcands),max(indices[lowcands]),1)+
     ifelse(any(highcands),min(indices[highcands]),length(x)))/2
    }
    if (is.null(weights)) weights<-rep(1,length(x))
    if (length(x)>1)
    equalweights<- all((weights[2:length(weights)]-
          weights[1:length(weights)-1])==0)
    else
    equalweights<-TRUE
    xna <- (is.na(x) | weights==0)
    if(na.rm) x <- x[!xna]
    else if(any(xna)) return(rep.int(NA,5))
    sortorder<-order(x)
    x <- x[sortorder]
    weights<-weights[sortorder]
    n <- sum(weights)
    if(n == 0) rep.int(NA,5)
    else {
    if (equalweights){
  d <- c(1, 0.5*floor(0.5*(n+3)), 0.5*(n+1),
       n+1-0.5*floor(0.5*(n+3)), n)
      }
  else {
    if(length(x)>1)
  d<-c(1,sapply(c(0.25*n,0.5*n,0.75*n),
       function(xxx)interpolatedindex(xxx,weights)),
       length(x))
  else
  d<-rep(1,5)
  }
    0.5*(x[floor(d)]+x[ceiling(d)])
  }
}

