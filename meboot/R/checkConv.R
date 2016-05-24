
checkConv <- function(y, bigx, trueb = 1, n999 = 999, nover = 5, 
  seed1 = 294, key = 0, trace = FALSE) 
{
  p <- NCOL(bigx)
  np <- p + 1

  bigx <- as.matrix(bigx) # 'p' columns

  if (trace)
  {
    cat("Number of coefficients including the intercept.\n") 
    print(np)
  }

  # basic regression 

  reg1 <- lm(y ~ bigx)

  if (trace)
  {
    print(reg1)
  }

  if (all(trueb == 0))
    trueb <- round(coef(reg1), 2)

  if (trace)
  {
    cat("True regression coefficients for simulation.\n")
    print(trueb)
  }

  if (length(trueb) != np) 
    stop("Inappropriate number of assumed true coefficient values.")

  if (nover == 0) 
    nover <- round(0.2*length(y), 0)

  bigt <- length(y) # sample size in data  
  ane <- rep(1, bigt) # column of ones for intercept
  bigX <- cbind(ane, bigx) # Note upper case X
  nstart <- bigt - nover + 1

  if (nstart < 0)
    stop("Error in input 'nover' to checkConv.")

  if (trace)
  {
    cat("Range of simulations over following sample sizes.\n")
    print(c(nstart, bigt))
    cat("Following coefficients will receive convergence evaluations.\n")
  }

  if (key == 0) 
    key <- seq(1, np)
  nk <- length(key)

  if (trace)
  {
    cat("key:\n")
    print(key)
  }

  if(nk > 5)
    warning("key regressors > 5 for detailed analysis.")

  # =matrix(rep(NA, nk*n999), nrow=n999, ncol=nk)  
  XnmX <- array(NA, dim = c(n999, nover, nk))  # place to store output
  #needs n999 rows for the CC package 
  #above is 3 dimensional matrix array

  sd1 <- sd(resid(reg1)) # to be used  
  #in the definition of sd of normal errors  
  #in the context of simulated errors
  set.seed(seed1)  
  epsn <- rnorm(bigt,sd=sd1)  
  ytrue  <-  bigX %*% trueb + epsn
  reg2 <- lm(ytrue ~ bigx) # OLS  
  olsb <- coef(reg2)

  if (trace)
  {
    print(olsb) # basic regression 
  }  

  if (trace)
  {
    cat("Compare OLS ytrue~bigx with assumed true coefficients.\n")
    print(cbind(olsb, trueb))
  }

  #prepare to simulate for np coefficients  
  #dat=array(NA, dim=c(np,n999,nover) )
  #set.seed(seed1)  #now lower case meboot data create
  meby <- meboot(x=ytrue, reps=n999)$ensem  

  mebigx <- array(NA, dim=c(bigt,n999,p) )
  for (k in seq(1, p))
  {
    mebigx[,,k] <- meboot(x=bigx[,k], reps=n999)$ensem  
  } # end creation of meboot resamples for all regressors

  allcoef <- matrix(rep(NA, np*n999), nrow=n999, ncol=np)
  nout <- 0 #initial n for the output Xn-X

  for (n in seq(nstart, bigt))
  { #begin n loop
    nout <- nout + 1 #start at 1 till nover  
    reg3 <- lm(ytrue[1:n] ~ bigx[1:n,]) #OLS  
    olsb <- coef(reg3)

    if (trace)
      print(round(olsb, 2))

    Meby <- meby[1:n,]#choose n rows starting at 2  
    #Meby is n by 999,  Mebigx is n by 999 by p
    Mebigx <- mebigx[1:n,,] #choose n rows of regressors
    #note first character is upper case M here

    Mebigx <- as.array(Mebigx)
    dim(Mebigx) <- c(n, n999, p)

    for (i in seq(1, n999))
    { # begin i loop to use meboot resamples  
      # output rows are resamples, hence notation i is used  
      allcoef[i,] <- coef(lm(Meby[,i] ~ Mebigx[,i,]))  
      #store the coefficients  all variables with upper case M in Meb
    }
 
    for (j in key)
    { #loop j for slope  coefficients
      #print(c("convergence for coefficient No.",j),quote=FALSE)  
      XnmX[,nout,j] <- allcoef[,j] - olsb[j]
    } # end j loop key coefficients selected for convergence study 
  } # end n loop

  #list(allcoef=allcoef, np=np, olsb=olsb, Xn=Xn, X=X)
  return(XnmX)  #returns (Xn - X) of ConvergenceConcepts pkg
}

