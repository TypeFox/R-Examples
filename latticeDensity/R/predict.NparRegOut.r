predict.NparRegOut <-
function(object,new.pred = NULL,...){
#
Zk <- as.vector(object$NparRegDenom)
Pk <- as.vector(object$NparRegNum)
if(is.null(new.pred)){
  which.nodes <- as.vector(object$which.nodes)
  predictions <- (Zk/Pk)[which.nodes]
  return(predictions)
  }else{
    temp <- bbox(rbind(new.pred,object$nodes))
    bound.vect <- c(temp[1,1],temp[1,2],temp[2,1],temp[2,2])
    X <- as.ppp(new.pred,W = bound.vect)
    Y <- as.ppp(object$nodes,W=bound.vect)
    closest <- nncross(X,Y)$which
    Pk[Pk==0] <- 1
    Zk[Pk==0] <- mean(Zk)*length(Zk)
    out <- (Zk/Pk)[closest]
    return(out)
  }
  }