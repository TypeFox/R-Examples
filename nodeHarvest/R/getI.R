getI <-
function(Z,X,Y=NULL,mode="mean"){
  I <- matrix(0, nrow= nrow(X), ncol=length(Z))
  if(!is.null(Y)){
    if(mode=="outbag"){
      for (ll in 1:length(Z)){
        ind <- getsamples(Z[[ll]],X,levelvec=attr(Z,"levelvec"))
        
        if(length(ind)>1){
          I[ind,ll] <- (sum(Y[ind])-Y[ind])/(length(ind)-1)
        }else{
          I[ind,ll] <- mean(Y)
        }
        
        attr(Z[[ll]],"mean") <- mean(Y[ind])
      } 
    }else{
      for (ll in 1:length(Z)){
        ind <- getsamples(Z[[ll]],X,levelvec=attr(Z,"levelvec"))
        if(length(ind) >1) I[ind,ll ] <- attr(Z[[ll]],"predict") <- mean(Y[ind]) else I[ind,ll ] <- mean(Y)
      }
    }
  }else{
    if(mode=="member"){
      for (ll in 1:length(Z)){
        ind <- getsamples(Z[[ll]],X,levelvec=attr(Z,"levelvec"))
        if(length(ind)>=1){
          I[ind,ll] <- 1
        }
      }
    }else{
      for (ll in 1:length(Z)){
        ind <- getsamples(Z[[ll]],X,levelvec=attr(Z,"levelvec"))
        if(length(ind)>=1){
          I[ind,ll] <- if(mode=="predict") attr(Z[[ll]],"predict") else attr(Z[[ll]],"mean")
        }
      }
    }
  }
  return(list(I=I,Z=Z))
}

