###########################################################################
####                         cclcda2.default                           ####
####  ===============================================================  ####
####  - estimation of one Latent-Class-Model for all classes           ####
####    with class conditional mixture weights                         ####
####  - determination of model selection criteria                      ####
###########################################################################

cclcda2.default <- function(       x,                   # dataframe containing the observations of the manifest variables
                                   grouping=NULL,       # vector containing the class membership
                                   prior=NULL,          # vector containing the prior probabilites for each class
                                   probs.start=NULL,    # list of matrices containing the starting values for the EM Algorithm
                                   wmk.start=NULL,      # matrix containing the starting values for the EM Algorithm
                                   nrep=1,              # number of repitions for the EM Algorithm
                                   m=3,                 # number of subclasses per class (integer or vector)
                                   maxiter = 1000,      # maximum number of iterations of the EM Algorithm
                                   tol=1e-10,           # limit of difference of likelihoods for stopping the EM Algorithm
                                   subset=1:nrow(x),    # a subset for the training data
                                   na.rm = FALSE,       # remove NA values
                                   ...
                                  )
{
  x <- x[subset,]
  if (na.rm==TRUE) {x <- x[apply(x[1:10,],1,function(z) !any(is.na(z))),]
                    grouping <- grouping[apply(x[1:10,],1,function(z) !any(is.na(z))),]
                   }
  if (any(is.na(grouping))) { warning("At least one of the grouping elements is NA! Observation(s) will be omitted!")
                              x <- x[which(is.na(grouping)),]
                              grouping <- grouping[which(is.na(grouping)),] 
                            }

  # if `grouping' vector not given, first data column is taken. 
  if (is.null(grouping)) {
    grouping <- data[,1]
    data <- data[,-1]
  }

  varnames <- colnames(x)
  k <- max(grouping, na.rm=TRUE) # number of groups 
  d <- ncol(x)                   # number of variables 
  n <- length(grouping)          # number of observations

  # number of categories per variable
  r <- as.numeric(matrix(apply(x,2,max, na.rm=TRUE)))

  npar <- m*(sum(r)-d+1)-1
  if (min(npar, n-1) > prod(r)) warning("LCM is not local identifiable! Please respecify model!")

  nk <- as.numeric(table(grouping))
  if (is.null(prior)) prior <- nk/n


  # all possible outcomes for each variable
  comb <- matrix(0, ncol=2, nrow=sum(r))
  comb[,1] <- rep(1:d,r)
  temp <- NULL
  for (i in 1:d) {temp <- c(temp, 1:r[i])}
  comb[,2] <- temp

  # number of times, that each of these combinations occurs
  ncomb <- list()
  for (i in 1:d)
  {
  ncomb[[i]] <- list()
  for (j in 1:r[i])
  {
   ncomb[[i]][[j]] <- which(x[,i]==j)
  }
  }

  # log-Likelihood
  llik <- function(x)
  {
   log(prior%*%apply(t(t(lca.wmk)*apply(theta^x.bin(x), 2, prod)),1,sum))
  }

  llik.f <- numeric()

  # function to get a binary x
  x.bin <- function(x=x)
  {
  res <- numeric(sum(r))
  for (i in 1:d)
  {
  for (j in 1:r[i])
  {
  if((is.na(x[i])==FALSE & x[i]==j))
  res[(sum(r[1:i]))-(r[i]-j)] <- 1
  }
  }
  return(res)
  }

  # function for computation of tau
  tau <- function(z)
  {
  temp <- t(t(lca.wmk)*apply(theta^x.bin(z), 2, prod))
  return(temp/apply(temp,1,sum))
  }

  lca.theta.list <- list()
  lca.wmk.list <- list()
  maxllik.list <- numeric()

  for (nr in 1:nrep)  # nrep repitions of the EM-Algorithm
  {

  if (is.null(probs.start)) 
     {
      probs <- list()
      for (j in 1:d) 
      { 
      probs[[j]] <- matrix(runif(m*r[j]),nrow=m,ncol=r[j])
      probs[[j]] <- probs[[j]]/rowSums(probs[[j]]) 
      }
     }else(probs <- probs.start)

  if (is.null(wmk.start)) 
     {
      wmk  <- matrix(runif(m*k),nrow=k,ncol=m)
      wmk <- wmk/rowSums(wmk)
     }else(wmk <- wmk.start)

  # convert lca.theta to a matrix
  lca.theta <- probs
  lca.wmk <- wmk
  counter <- 0

repeat{

  counter <- counter+1

  theta <- unlist(lca.theta, use.names=FALSE)
  theta <- matrix(theta, ncol=m, byrow=TRUE)

  tau.est <- apply(x,1,tau)  #estimation of tau

  wmk.est <- matrix(0,nrow=k, ncol=m)
  for (i in 1:k)
  {
  wmk.est[i,] <- 1/nk[i] * apply(tau.est[seq(i,m*k,k),grouping==i],1,sum)
  }

  # estimation of theta
  theta.est <- list()
  for (i in 1:d)
  {
  theta.est[[i]] <- matrix(0,nrow=m, ncol=r[i])
  }

  for (j in 1:d)
  {
  for (l in 1:r[j])
  {
   theta.est[[j]][,l] <- apply(matrix(apply(tau.est[,ncomb[[j]][[l]]],1,sum),nrow=k),2,sum)/apply(matrix(apply(tau.est,1,sum),nrow=k),2,sum)
  }
  }

  lca.wmk <- wmk.est
  lca.theta <- theta.est

  # calculation of the log-likelihood
  llik.f[counter] <- sum(apply(x,1,llik))

  if (counter > 1)  { if(((llik.f[counter]-llik.f[counter-1])<tol) | counter > maxiter) break }
}

  maxllik <- llik.f[counter]




  lca.theta.list[[nr]] <- lca.theta
  lca.wmk.list[[nr]] <- lca.wmk
  maxllik.list[nr] <- maxllik
  best <- which.max(maxllik.list)  # best model

  # output like poLCA
  if(nrep>1) cat("Model",nr ,": llik =" , maxllik, " ... best llik = ", maxllik.list[best], "\n", sep="") 
  probs.start <- NULL
  wmk.start <- NULL
}



  lca.theta <- lca.theta.list[[best]]
  lca.wmk <- lca.wmk.list[[best]]
  maxllik <- maxllik.list[best]

  # computation of AIC and BIC
  aic <- -2*maxllik + 2*npar
  bic <- -2*maxllik + npar*log(n)

  lca.wkm <- t(t(lca.wmk * prior) / apply(lca.wmk * prior, 2, sum))
  lca.w <- prior%*%lca.wmk

  # sort x for computation of Gsq and Chisq
  xsorted <- x[do.call(order,data.frame(x)),]

  datamat <- xsorted[1,]
  freq <- numeric()
  freq[1] <- 1
  curpos <- 1

  for (i in 2:n) {
        if (sum(xsorted[i,] == xsorted[i-1,], na.rm=TRUE)==d) {
            freq[curpos] <- freq[curpos]+1
        } else {
            datamat <- rbind(datamat,xsorted[i,])
            freq <- c(freq,1)
            curpos <- curpos+1
        }
    }

 
  # Gsq und Chisq
  # funtion to calculate the probability of each outcome unter the estimated model
  expn <- function(x)
  {
  return(n*lca.w%*%apply(theta^x.bin(x), 2, prod))
  }

  fit <- apply(datamat, 1, expn)
  lca.Chisq <- sum((freq-fit)^2/(fit)) + (n-sum(fit))
  lca.Gsq <- 2*sum(freq*log(freq/fit))

  # entropy
  entro <- lca.w %*% apply(lca.wkm, 2, entropy)

  # computation of the gini-coefficient
  gini <- lca.w %*% (1-apply(lca.wkm^2, 2, sum))

  #chi-square-test
  chi <- chisq.test(n*(lca.wmk * prior))
  chi.stat <- chi$statistic
  chi.p <- chi$p.value

  # output of result
  result <- list(call=match.call(),    # function call
                 lca.theta=lca.theta,  # estimates of class conditional response probabilities
                 lca.w=lca.w,          # estimates of mixing proportions
                 lca.wmk=lca.wmk,      # estimates of the class conditional mixing proportions
                 prior=prior,          # class prior probabilites
                 m=m,                  # number of subclasses per class
                 r=r,                  # number of possible outcomes per variable
                 k=k,                  # number of classes
                 d=d,                  # number of variables
                 aic=aic,              # AIC
                 bic=bic,              # BIC
                 Gsq=lca.Gsq,          # Gsq
                 Chisq=lca.Chisq,      # Chisq
                 entropy=entro,        # Entropy
                 gini=gini,            # Gini
                 chi.stat=chi.stat,    # Chi^2-statistic
                 chi.p=chi.p           # Chi^2 p-value
)
  class(result) <- "cclcda2"
  return(result)

}


## print.cclcda2: readable output

print.cclcda2 <- function(x, ...)
{
    cl <- paste("class", 1:x$k, sep=" ")
    cat("Call:\n")
    print(x$call, ...)
    cat("\nNumber of classes: ", x$k, "\n")
    cat("\nPrior probabilites:\n")
    for (i in 1:x$k)
    cat(paste(cl[i], round(x$prior[i],2), sep=": "), "\n")
    cat("\nNumber of latent classes: ", x$m, "\n")
    cat("\nAIC: ", round(x$aic,2), "\n")
    cat("\nBIC: ", round(x$bic,2), "\n")
    cat("\nGsq: ", round(x$Gsq,2), "\n")
    cat("\nChisq: ", round(x$Chisq,2), "\n")
    cat("\nEntropy: ", round(x$entropy,2), "\n")
    cat("\nGini: ", round(x$gini,2), "\n")
    cat("\nChisq-statistic: ", round(x$chi.stat,2), "\n")
    cat("\nChisq p-value: ", round(x$chi.p,4), "\n")
    invisible(x)
} 















