#estimate bounds on processing
#RT should be a list of response times, with each entry corresponding to a single
#individual processing channel
#CR is the list of correct/incorrect logical indicators, with size matching RT

#if assume.ID=TRUE, then RT only needs to be a vector or the first entry of RT is used


estimate.bounds <- function (RT, CR = NULL, stopping.rule = c("OR","AND","STST"), assume.ID=FALSE,
                             numchannels=NULL, unified.space=FALSE) 
{
  rule <- match.arg(stopping.rule, c("OR","AND","STST"))
  
  nconds <- length(RT) #total number of input RTs
  n.channels <- numchannels
  #if numchannels for ID model not specified, what is number of channels in the UCIP model?
  if (is.null(numchannels) & !assume.ID) {
    n.channels <- nconds
  }
  
  #if n.channels <2, produce an error message that bounds cannot be 
  #estimated for less than 2 channels
  if (n.channels < 2){
    stop("Too few channels specified. Number of channels must be >=2 to estimate CDF bounds.")
  }
  #create correct/incorrect data if not provided
  if (is.null(CR) | length(CR) != nconds) {
    for (i in 1:nconds) {
      CR[[i]] <- rep(1, length(RT[[i]]))
    }
  }
  
  #find overall range of RTs
  times <- sort(unique(c(RT, recursive = TRUE)))
  nt <- length(times)
  
  #establish CDFs and survivor functions for each single target distribution
  G <- vector('list',nconds)
  Gmat <- matrix(rep(0,nconds*nt), ncol=nconds, nrow=nt)
  Smat <- matrix(rep(0,nconds*nt), ncol=nconds, nrow=nt)
  for (j in 1:nconds){
    RTx <- sort(RT[[j]], index.return = TRUE)
    RTj <- RTx$x
    CRj <- as.logical(CR[[j]])[RTx$ix]
    G[[j]] <- ecdf(RTj[CRj])
    Gmat[,j] <- G[[j]](times)
    Smat[,j] <- 1-Gmat[,j]
  }
  
  #given the stopping rule, estimate the CDF bounds
  if (rule == "OR"){  
    #two channel bounds
    if (n.channels == 2){
      if (assume.ID){
        if (unified.space){
          upper <- 2*Smat[,1]-1
          lower <- Smat[,1]
        }
        else {upper <- 2*Gmat[,1]
              lower <- Gmat[,1]}
      }
      #not ID, 2 channels
      else {
        if (unified.space){
          upper <- apply(Smat,1,sum)-1
          lower <- apply(Smat,1,min)
        }
        else{upper <- apply(Gmat, 1, sum)
             lower <- apply(Gmat, 1, max)} 
      }#end notID, 2 channels
    }#end (n.channels==2)
    #more than two channels
    if (n.channels > 2){
      if (assume.ID){
        #upper is 2Fn-1(t)-Fn-2(t)
        #lower is Fn-1(t)
        if (unified.space){
          upper <- 2*(Smat[,1]^(n.channels-1)) - Smat[,1]^(n.channels-2)
          lower <- Smat[,1]^(n.channels-1)
        }
        #note that min time Fn(t)=1-Gn(t) is 1-max time
        else {upper <- 2*(1-Gmat[,1]^(n.channels-1)) - (1-Gmat[,1]^(n.channels-2))
              lower <- 1-(Gmat[,1]^(n.channels-1)) }
      }  #end (assume.ID) 
      #not ID, >2 channels
      else {
        #find all n-1 distributions
        Fi.array <- matrix(rep(0,nt*nconds), nrow=nt, ncol=nconds)
        for (i.out in 1:nconds){
          Fi.array[,i.out] <- 1-apply(as.matrix(1-Gmat[,-i.out]), 1, prod)
        }
        i.out <- NULL
        
        num.cols <- (nconds*(nconds-1))/2
        Fij.array <- matrix(rep(0,nt*num.cols), nrow=nt, ncol=num.cols)
        Fij.array2 <- matrix(rep(0,nt*num.cols), nrow=nt, ncol=num.cols)
        column <- 1
        for (i.out in 1:(nconds-1)){
          for (j.out in (i.out+1):nconds){
            #CDF values for unified space
            Fij.array[,column] <- (1-apply(as.matrix(1-Gmat[,-i.out]), 1, prod)) + (1-apply(as.matrix(1-Gmat[,-j.out]), 1, prod)) - (1-apply(as.matrix(1-Gmat[,c(-i.out, -j.out)]), 1, prod))
            #survivor values for non-unified space
            Fij.array2[,column] <- apply(as.matrix(1-Gmat[,-i.out]), 1, prod) + apply(as.matrix(1-Gmat[,-j.out]), 1, prod) - apply(as.matrix(1-Gmat[,c(-i.out, -j.out)]), 1, prod)
            column <- column+1
          }
        if (unified.space){
          upper <- apply(Fij.array2, 1, max)
          lower <- apply(1-Fi.array, 1, min)
        }
        else{
          #min{i,j}[1-G^(i)_(n-1)(t) + 1-G^(j)_(n-1)(t) - 1-G^(i,j)_(n-2)(t)] for i!=j
          upper <- apply(Fij.array, 1, min)
          #max{i} F^(i)_(n-1)(t) = 1-G^(i)_(n-1)(t) -- max of all n-1 cdfs
          lower <- apply(Fi.array, 1, max)
          }
        }#end else [ie. not unified.space]
      }#end else [ie. not assume.ID]
     }#end (n.channels >2)
      if (unified.space){
        ucip <- apply(Smat, 1, prod)
        upper <- log(upper)/log(ucip)
        lower <- log(lower)/log(ucip)
      }
    }#end (rule=="OR")
    
    else if (rule == "AND"){
      if (n.channels == 2) {
        if (assume.ID){
          upper <- Gmat[,1]
          lower <- 2*Gmat[,1] - 1
        }
        else{upper <- apply(Gmat, 1, min)
             lower <- apply(Gmat, 1, sum) - 1
        }
      } #end (n.channels==2)
      if (n.channels > 2){
        if (assume.ID){
          #upper is G_(n-1)(t)
          #lower is 2G_(n-1)(t)-G_(n-2)(t)
          upper <- Gmat[,1]^(n.channels-1)
          lower <- (2*(Gmat[,1]^(n.channels-1))) - Gmat[,1]^(n.channels-2)
        }
        else {
          #min{i} G^(i)_(n-1)(t)  -- min of all n-1 cdfs
          Gi.array <- matrix(rep(0,nt*nconds), nrow=nt, ncol=nconds)
          for (i.out in 1:nconds){
            Gi.array[,i.out] <- apply(as.matrix(Gmat[,-i.out]), 1, prod)
          }
          upper <- apply(Gi.array, 1, min)
          i.out <- NULL
          
          #max{i,j}[G^(i)_(n-1)(t) + G^(j)_(n-1)(t) - G^(i,j)_(n-2)(t)] for i!=j
          num.cols <- (nconds*(nconds-1))/2
          Gij.array <- matrix(rep(0,nt*num.cols), nrow=nt, ncol=num.cols)
          column <- 1
          for (i.out in 1:(nconds-1)){
            for (j.out in (i.out+1):nconds){
              Gij.array[,column] <- apply(as.matrix(Gmat[,-i.out]), 1, prod) + apply(as.matrix(Gmat[,-j.out]), 1, prod) - apply(as.matrix(Gmat[,c(-i.out, -j.out)]), 1, prod)
              column <- column+1
            }
          }
          lower <- apply(Gij.array, 1, max)
        }
      } #end (n.channels>2)
      if (unified.space){
        ucip <- apply(Gmat, 1, prod)
        upper <- log(ucip)/log(upper)
        lower <- log(ucip)/log(lower)
      }
    }
    else if (rule == "STST"){
      if (assume.ID){
        upper <- n.channels * Gmat[,1]
        lower <- Gmat[,1]^(n.channels)
      }
      else {upper <- apply(Gmat, 1, sum)
            lower <- apply(Gmat, 1, prod)     
      }
      
      if (unified.space){
        ucip <- Gmat[,1]
        upper <- log(ucip)/log(upper)
        lower <- log(ucip)/log(lower)
      }
    }
    #create an interpolation function for each bound
    upper <- approxfun(times, upper)
    lower <- approxfun(times, lower)
    #return the bounds as a list
    return(list(Upper.Bound=upper, Lower.Bound=lower))
  }
