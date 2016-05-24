###########################################################################
####                         cclcda.default                            ####
####  ===============================================================  ####
####  - estimation of one Latent-Class-Model for all classes           ####
####  - estimation of probabilities for class membership               ####
####  - determination of model selection criteria                      ####
###########################################################################



cclcda.default <- function(
                                   x,                   # dataframe containing the observations of the manifest variables
                                   grouping=NULL,       # vector containing the class membership
                                   prior=NULL,          # vector containing the prior probabilites for each class
                                   probs.start=NULL,    # list of matrices containing the starting values for the EM Algorithm
                                   nrep=1,              # number of repitions for the EM Algorithm
                                   m=3,                 # number of subclasses per class (integer or vector)
                                   maxiter = 1000,      # maximum number of iterations of the EM Algorithm
			 	   tol = 1e-10,         # tolerance for aborting the EM Algorithm
                                   subset=1:nrow(x),    # a subset for the training data
                                   na.rm = FALSE,       # remove NA values
                                   ...
                                  )
{

  data <- x[subset,]  # define the data
  rm(x)

  # save variable names
  varnames <- colnames(data)

  # if `grouping' vector not given, first data column is taken. 
  if (is.null(grouping)) {
    grouping <- data[,1]
    data <- data[,-1]
  }

  k <- max(grouping, na.rm=TRUE) # number of groups 
  d <- ncol(data)                # number of variables 
  n <- length(grouping)          # number of observations

  # number of categories per variable
  r <- as.numeric(apply(data, 2, function(y) length(table(y))))  

  # if one class does not contain all possible observations stop procedure
  r.temp <- list()
  for (i in 1:k)
  {
  r.temp[[i]] <- as.numeric(apply(data[grouping==i,], 2, function(y) length(table(y))))
  }
  if(any(sapply(r.temp, function(z) all(z==r))==FALSE)) 
  stop("There is at least one class that does not contain all possible observations!", call.=FALSE)

  # in the commom components case m must be an integer
  if(length(m)!=1)
  stop("'m' must be an integer!")

  # in the commom components case m must be > 1
  if(m==1)
  stop("'m' must be > 1")


  # frequencies as prior
  if (is.null(prior)) { 
    prior <- tabulate(grouping)
    prior <- prior / sum(prior) 
  } 

  # formula for poLCA
  vars <- paste(varnames, sep="", collapse=",")
  f <- as.formula(paste("cbind(",vars,")~1",sep=""))

  # result of the LCA
  lca.res <- poLCA(f, data=data, nclass=m, maxiter = maxiter, na.rm = na.rm, probs.start = probs.start, nrep = nrep, verbose=FALSE)
  # estimates of subclass weights
  lca.w <- lca.res$P
  # estimates of parameters
  lca.theta <- lca.res$probs
  # values of the aic-criterion
  lca.aic <- lca.res$aic
  # values of the bic-criterion
  lca.bic <- lca.res$bic
  # Likelihood ratio/deviance statistic
  lca.Gsq <- lca.res$Gsq
  # Pearson Chi-square goodness of fit statistic for fitted vs. observed multiway tables
  lca.Chisq <- lca.res$Chisq

  lca.wmk <- matrix(0, ncol=m, nrow=k)
  for (i in 1:k)
  {
   lca.wmk[i,] <- (1/sum(grouping==i))*apply(lca.res$posterior[grouping==i,], 2, sum)
  }


  lca.wkm <- t(t(lca.wmk * prior) / apply(lca.wmk * prior, 2, sum))

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
                 lca.wmk=lca.wmk,       # estimates of the class conditional mixing proportions
                 prior=prior,          # class prior probabilites
                 m=m,                  # number of subclasses per class
                 r=r,                  # number of possible outcomes per variable
                 k=k,                  # number of classes
                 d=d,                  # number of variables
                 aic=lca.aic,          # AIC
                 bic=lca.bic,          # BIC
                 Gsq=lca.Gsq,          # Gsq
                 Chisq=lca.Chisq,      # Chisq
                 entropy=entro,        # Entropy
                 gini=gini,            # Gini
                 chi.stat=chi.stat,     # Chi^2-statistic
                 chi.p=chi.p           # Chi^2 p-value
)
  class(result) <- "cclcda"
  return(result)
}


## print.lcda: readable output

print.cclcda <- function(x, ...)
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
