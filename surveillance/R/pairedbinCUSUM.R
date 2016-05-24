######################################################################
# Compute ARL for paired binary CUSUM charts as introducted in Steiner,
# Cook and Farefwell, 1999, Monitoring paired binary surgical outcomes,
# Stats in Med, 18, 69-86.
#
# This code is an R implementation of Matlab code provided by
# Stefan H. Steiner, University of Waterloo, Canada.
#
# Params:
#  p       - vector giving the probability of the four different possibilities
#            c((death=0,near-miss=0),(death=1,near-miss=0),
#            (death=0,near-miss=1),(death=1,near-miss=1))
#  w1, w2  - w1 and w2 are the sample weights vectors for the two CUSUMs.
#            (see (2)). We have w1 is equal to deaths
#            (according to paper it being 2 would be more realistic)
#  h1, h2  - decision barriers for the individual cusums (see (3))
#  h11,h22 - joint decision barriers (see (3))
#  sparse - use Matrix package
######################################################################

pairedbinCUSUM.runlength <- function(p,w1,w2,h1,h2,h11,h22, sparse=FALSE) {
  #Size of the sparse matrix -- assumption h1>h11 and h2>h22
  mw <- h1*h22+(h2-h22)*h11;

  cat("g =",mw+3,"\n")
    
  #build transition matrix; look at current state as an ordered pair (x1,x2)
  #the size of the matrix is determined by h1, h2, and h11 and h22
  #Look at all 3 possible absorbing conditions
  transm <- matrix(0, mw+3, mw+3)

  #the last row/column is the absorbing state, I_{3\times 3} block
  #Is this ever used??
  transm[mw+1,mw+1] <- 1
  transm[mw+2,mw+2] <- 1
  transm[mw+3,mw+3] <- 1

  #go over each row and fill in the transition probabilities
  for (i in 1:mw) {
#        cat(i," out of ", mw,"\n")
        
    #find the corresponding state
    if (i>h1*h22) {
      temp <- floor((i-h1*h22-1)/h11)
      x1 <- i-h1*h22-1-temp*h11
      x2 <- temp+h22
    } else {
      x2 <- floor((i-1)/h1);
      x1 <- i-x2*h1-1;
    }
    
    #go over the four different weight combinations
    for (j in 1:2) {
      for (k in 1:2) {
        x1n <- x1+w1[j+2*(k-1)]   #death chart
        x2n <- x2+w2[k]  #look at all possible combinations of weights
        #we cant go below zero
        if (x1n<0) { x1n <- 0 }
        if (x2n<0) { x2n <- 0 }
        newcond=0;
   	#try to figure out what condition index the new CUSUM values correspond to
        if (x1n>=h1) {
          newcond <- mw+1;     #absorbing state on x1
        } else {
          if (x2n>=h2) {
            newcond <- mw+2  #absorbing state on x2
          } else {
            if ((x1n>=h11)&(x2n>=h22)) {    
              #only register this if other two conditions are not satisfied
              newcond <- mw+3
            }
          }
        }
        
        if (newcond==0) {   #transition is not to an absorbing state
	  #translate legal ordered pair to state number
          if (x2n<h22) {
            newcond <- x2n*h1+x1n+1
          } else {
            newcond <- h1*h22+(x2n-h22)*h11+x1n+1
          }
        }
        
        #Update transition matrix value
        transm[i,newcond] <- transm[i,newcond]+p[(k-1)*2+j];
      }
    }
  }

  #remove rows and columns corresponding to absorbing boundary
  r <- transm[1:mw,1:mw]

  #find the average run length according to (8) in paper
  id <- diag(rep.int(1,mw))

  #Use sparse matrix computations (or just ordinary matrix
  #computations) for the averafe run-length
  mom <- if (!sparse) { # non-sparse computing
    .rowSums(solve(id-r), mw, mw)
  } else { #sparse-computing
    rowSums(solve(Matrix(id-r)))  # we import the required Matrix methods
  }
  arl <- mom[1]
  
  #Dummy return, as calculation of p-vector crashes on Mac OS X
  return(arl)
}

######################################################################
# Paired binary CUSUM method as described by Steiner et al (1999).
# This is the workhorse and is wrapped into a nice function
# being able to handle input as an S4 sts object in pairedbinCUSUM.
#
# Parameters:
#   x      - data
#   theta0 - in-control parameters of the paired-binary model
#   theta1 - out-of-control parameters of the paired-binary model
#   h1     - primary limit for 1st series
#   h2     - primary limit for 2nd series
#   h11    - secondary limit for 1st series
#   h22    - secondary limit for 2nd series
######################################################################

pairedbinCUSUM.LLRcompute <- function(x,theta0, theta1, h1,h2,h11,h22) {
  #Initialize variables
  t <- 0
  alarm <- c(FALSE,FALSE)
  stopped <- FALSE
  S <- matrix(0,nrow(x)+1,2)

  #Run the paired binomial CUSUM
  while (!stopped) {
    #Increase time
    t <- t+1

    #Compute log likelihood ratios
    lr1 <- (theta1[1] - theta0[1])*x[t,1] + log(1+exp(theta0[1])) - log(1+exp(theta1[1]))
    lr2 <- (theta1[2]-  theta0[2])*x[t,2] + log(1+exp(theta0[2]+theta0[3]*x[t,1])) - log(1+exp(theta1[2]+theta1[3]*x[t,1]))
    #Run CUSUMs
    S[t+1,1] <- max(0, S[t,1] + lr1)
    S[t+1,2] <- max(0, S[t,2] + lr2)

    #Do we have an alarm? Could be caused by primary limits or by
    #secondary limits
    alarm <- c(S[t+1,1] > h1, S[t+1,2] > h2)
    if ((S[t+1,1] > h11) & (S[t+1,2] > h22)) { alarm <- c(TRUE,TRUE) }
#    alarm <- (S[t+1,1] > h1) | (S[t+1,2] > h2) |
#             ((S[t+1,1] > h11) & (S[t+1,2] > h22))
    #If one or both of the CUSUMs produced an alarm then stop
    if ((sum(alarm)>0) | (t==nrow(x))) { stopped <- TRUE}
  }
  return(list(N=t,val=S[-1,],alarm=alarm))
}

######################################################################
# STS wrapper for the Paired binary CUSUM method. This follows in
# style the categoricalCUSUM method.
######################################################################

pairedbinCUSUM <- function(stsObj, control = list(range=NULL,theta0,theta1,h1,h2,h11,h22)) {
   # Set the default values if not yet set
  if(is.null(control[["range",exact=TRUE]])) { 
    control$range <- 1:nrow(observed(stsObj))
  }
  if(is.null(control[["theta0",exact=TRUE]])) { 
    stop("Error: No specification of in-control parameters theta0!")
  }
  if(is.null(control[["theta1",exact=TRUE]])) { 
    stop("Error: No specification of out-of-control parameters theta1!")
  }
  if(is.null(control[["h1",exact=TRUE]])) { 
    stop("Error: No specification of primary threshold h1 for first series.")
  }
  if(is.null(control[["h2",exact=TRUE]])) { 
    stop("Error: No specification of primary threshold h2 for 2nd series.")
  }
  if(is.null(control[["h11",exact=TRUE]])) { 
    stop("Error: No specification of secondary limit h11 for 1st series.")
  }
  if(is.null(control[["h22",exact=TRUE]])) { 
    stop("Error: No specification of secondary limit h11 for 2nd series.")
  }

  #Extract the important parts from the arguments
  range <- control$range
  y <- stsObj@observed[range,,drop=FALSE]
  theta0 <- control[["theta0",exact=TRUE]]
  theta1 <- control[["theta1",exact=TRUE]]
  h1 <- control[["h1",exact=TRUE]]
  h2 <- control[["h2",exact=TRUE]]
  h11 <- control[["h11",exact=TRUE]]
  h22 <- control[["h22",exact=TRUE]]
  
  #Semantic checks.
  if (ncol(y) != 2) {
    stop("Error: The number of columns in the sts object needs to be two.")
  }

  #Reserve space for the results. Contrary to the categorical CUSUM
  #method, each ROW represents a series.
  alarm <- matrix(data = 0, nrow = length(range), ncol = ncol(y))
  upperbound <- matrix(data = 0, nrow = length(range), ncol = ncol(y))
  
  #Setup counters for the progress
  doneidx <- 0
  N <- 1
  noofalarms <- 0
  noOfTimePoints <- length(range)

  #######################################################
  #Loop as long as we are not through the entire sequence
  #######################################################
  while (doneidx < noOfTimePoints) {
     #Run paired binary CUSUM until the next alarm
    res <- pairedbinCUSUM.LLRcompute(x=y, theta0=theta0, theta1=theta1, h1=h1, h2=h2, h11=h11, h22=h22)
  
    #In case an alarm found log this and reset the chart at res$N+1
    if (res$N < nrow(y)) {
      #Put appropriate value in upperbound
      upperbound[1:res$N + doneidx,]  <- res$val[1:res$N,]
      alarm[res$N + doneidx,] <- res$alarm
    
     #Chop & get ready for next round
      y <- y[-(1:res$N),,drop=FALSE]
#      theta0 <- pi0[,-(1:res$N),drop=FALSE]
#      theta1 <- pi1[,-(1:res$N),drop=FALSE]
#      n <- n[-(1:res$N)]
      
      #Add to the number of alarms
      noofalarms <- noofalarms + 1
    }
    doneidx <- doneidx + res$N
  }

  #Add upperbound-statistic of last segment, where no alarm is reached
  upperbound[(doneidx-res$N+1):nrow(upperbound),]  <- res$val
  
  # Add name and data name to control object
  control$name <- "pairedbinCUSUM"
  control$data <- NULL #not supported anymore

  #New direct calculations on the sts object
  stsObj@observed <- stsObj@observed[control$range,,drop=FALSE]
  stsObj@state <- stsObj@state[control$range,,drop=FALSE]
  stsObj@populationFrac <- stsObj@populationFrac[control$range,,drop=FALSE]
  stsObj@alarm <- alarm
  stsObj@upperbound <- upperbound

  #Fix the corresponding start entry
  start <- stsObj@start
  new.sampleNo <- start[2] + min(control$range) - 1
  start.year <- start[1] + (new.sampleNo - 1) %/% stsObj@freq 
  start.sampleNo <- (new.sampleNo - 1) %% stsObj@freq + 1
  stsObj@start <- c(start.year,start.sampleNo)

  #Done
  return(stsObj)
}

