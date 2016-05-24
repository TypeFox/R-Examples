# bestfit defines the numer of different random 100% lenght dataset over which
# the models are fitted and tested before the final modelset is chosen 

bestfitmod <- function(indnr, paramnr, xv, yv, ch, zv, vv)
{ 
  if (indnr == 2)
  {
    nmodelterms = paramnr;  
    nterms = 17;
    # 100 % of the data is used, may be adjusted if the user would like to 
    # use different training and testing data
    frac = 1;
    
    SEtest <- matrix(100, nterms, choose(nterms, nmodelterms)) 
    # should be uncommented if the user intends to use differen training and testing data
    #SEtest_all <- matrix(100, nterms, choose(nterms, nmodelterms)) 
    
    # defining training and testing dataset
    # sample returns a row vector containing a random permutation of the integers form
    # one to length(x) inclusive
    ridx <- sample(length(xv))
    xv = xv[ridx]
    yv = yv[ridx]
    ch = ch[ridx]
    xtrain <- xv[1:ceiling(frac*length(xv))]
    ytrain <- yv[1:ceiling(frac*length(yv))]
    chtrain <- ch[1:ceiling(frac*length(ch))]
    # please uncomment if using separate training and testing set 
    #   xtest <- xv[ceiling(frac*length(xv)+1:length(xv))]
    #   ytest <- yv[ceiling(frac*length(yv)+1:length(yv))]
    #   chtest <- ch[ceiling(frac*length(ch)+1:length(ch))]
    
    # Fit and test every combination of 17 terms
    for (i in 1:nmodelterms)
    {
      M <- combs(1:nterms, i)
     
      for (j in 1:nrow(M)) 
      {
        tmp = polyfitreg(indnr, xtrain, ytrain, chtrain, M[j,])
        B <- tmp[[1]]
        err <- tmp[[2]]
        
        SEtest[i, j] <- err
      }  
    }
    return(SEtest) 
  }
  
  ############################# indnr == 2 ends, indnr == 3 begins #################################
  
  if (indnr == 3)
  {
    nmodelterms = paramnr;  
    nterms = 39;
    # 100 % of the data is used, may be adjusted if the user would like to 
    # use different training and testing data
    frac = 1;
    
    SEtest <- matrix(100, nterms, choose(nterms, nmodelterms)) 
    # should be uncommented if the user intends to use differen training and testing data
    #SEtest_all <- matrix(100, nterms, choose(nterms, nmodelterms)) 
    
    # defining training and testing dataset
    # sample returns a row vector containing a random permutation of the integers form
    # one to length(x) inclusive
    ridx <- sample(length(xv))
    xv = xv[ridx]
    yv = yv[ridx]
    zv = zv[ridx]
    ch = ch[ridx]
    xtrain <- xv[1:ceiling(frac*length(xv))]
    ytrain <- yv[1:ceiling(frac*length(yv))]
    ztrain <- zv[1:ceiling(frac*length(zv))]
    chtrain <- ch[1:ceiling(frac*length(ch))]
    # please uncomment if using separate training and testing set 
    #   xtest <- xv[ceiling(frac*length(xv)+1:length(xv))]
    #   ytest <- yv[ceiling(frac*length(yv)+1:length(yv))]
    #   ztest <- zv[ceiling(frac*length(zv)+1:length(zv))]
    #   chtest <- ch[ceiling(frac*length(ch)+1:length(ch))]
    
    # Fit and test every combination of 39 terms
    for (i in 1:nmodelterms)
    {
      M <- combs(1:nterms, i)
      
      for (j in 1:nrow(M)) 
      {
        tmp = polyfitreg(indnr, xtrain, ytrain, chtrain, M[j,], ztrain)
        B <- tmp[[1]]
        err <- tmp[[2]]
        
        SEtest[i, j] <- err
      }  
    }
    return(SEtest) 
  }
  
  ############################# indnr == 3 ends, indnr == 4 begins #################################
  
  if (indnr == 4)
  {
    nmodelterms = paramnr;  
    nterms = 97;
    # 100 % of the data is used, may be adjusted if the user would like to 
    # use different training and testing data
    frac = 1;
    
    SEtest <- matrix(100, nterms, choose(nterms, nmodelterms)) 
    # should be uncommented if the user intends to use differen training and testing data
    #SEtest_all <- matrix(100, nterms, choose(nterms, nmodelterms)) 
    
    # defining training and testing dataset
    # sample returns a row vector containing a random permutation of the integers form
    # one to length(x) inclusive
    ridx <- sample(length(xv))
    xv = xv[ridx]
    yv = yv[ridx]
    zv = zv[ridx]
    vv = vv[ridx]
    ch = ch[ridx]
    xtrain <- xv[1:ceiling(frac*length(xv))]
    ytrain <- yv[1:ceiling(frac*length(yv))]
    ztrain <- zv[1:ceiling(frac*length(zv))]
    vtrain <- vv[1:ceiling(frac*length(vv))]
    chtrain <- ch[1:ceiling(frac*length(ch))]
    # please uncomment if using separate training and testing set 
    #   xtest <- xv[ceiling(frac*length(xv)+1:length(xv))]
    #   ytest <- yv[ceiling(frac*length(yv)+1:length(yv))]
    #   ztest <- zv[ceiling(frac*length(zv)+1:length(zv))]
    #   chtest <- ch[ceiling(frac*length(ch)+1:length(ch))]
    
    # Fit and test every combination of 97 terms
    for (i in 1:nmodelterms)
    {
      M <- combs(1:nterms, i)
      
      for (j in 1:nrow(M)) 
      {
        tmp = polyfitreg(indnr, xtrain, ytrain, chtrain, M[j,], ztrain, vtrain)
        B <- tmp[[1]]
        err <- tmp[[2]]
        
        SEtest[i, j] <- err
      }  
    }
    return(SEtest) 
  }
}