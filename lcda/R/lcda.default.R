###########################################################################
####                         lcda.default                              ####
####  ===============================================================  ####
####  - estimation of one Latent-Class-Model per class                 ####
####  - determination of model selection criteria                      ####
###########################################################################



lcda.default <- function(
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

  varnames <- colnames(data)  # save variable names

  # if `grouping' vector not given, first data column is taken. 
  if (is.null(grouping)) {
    grouping <- data[,1]
    data <- data[,-1]
  }

  k <- max(grouping, na.rm=TRUE) # number of groups 
  d <- ncol(data)                # number of variables 
  n <- length(grouping)          # number of observations

  # if m is given as integer make a vector of m.
  if(length(m)==1){m <- rep(m, k)}
  # if m i not of length 1 or m stop procedure.
  if(length(m)!=1 & length(m)!=k) stop("'m' must be of length 1 or 'k'", call.=FALSE)
  # minimum for m is 2.
  if(any(m<2)) na.rm <- TRUE

  # number of categories per variable
  r <- as.numeric(apply(data, 2, function(y) length(table(y))))

  # if one class does not contain all possible observations stop procedure.
  r.temp <- list()
  for (i in 1:k)
  {
  r.temp[[i]] <- as.numeric(apply(data[grouping==i,], 2, 
                            function(y) length(table(y))))
  }
  if(any(sapply(r.temp, function(z) all(z==r))==FALSE)) 
  stop("There is at least one class that does not contain all 
       possible observations!"
       , call.=FALSE)

  # frequencies as prior
  if (is.null(prior)) { 
    prior <- tabulate(grouping)
    prior <- prior / sum(prior) 
  } 

  # formula for poLCA
  vars <- paste(varnames, sep="", collapse=",")
  f <- as.formula(paste("cbind(",vars,")~1",sep=""))

  # data and m per class as a list
  datam <- list()
  for (i in 1:k)
  {
  datam[[i]] <- list()
  datam[[i]][[1]] <- data[grouping==i,]
  datam[[i]][[2]] <- m[i]
  }

  # result of the LCA
  lca.res <- lapply(datam, function(z) 
                    poLCA(f, data = z[[1]], nclass=z[[2]], 
                          maxiter = maxiter, na.rm = na.rm, 
                          probs.start = probs.start, nrep = nrep, 
                          verbose = FALSE, tol=tol))


  # estimates of subclass weights
  lca.w <- lapply(lca.res, function(z) z$P)
  # estimates of parameters
  lca.theta <- lapply(lca.res, function(z) z$probs)
  # values of the aic-criterion
  lca.aic <- sapply(lca.res, function(z) z$aic)
  # values of the bic-criterion
  lca.bic <- sapply(lca.res, function(z) z$bic)
  # Likelihood ratio/deviance statistic
  lca.Gsq <- sapply(lca.res, function(z) z$Gsq)
  # Pearson Chi-square goodness of fit statistic for fitted vs. observed multiway tables
  lca.Chisq <- sapply(lca.res, function(z) z$Chisq)

  # output of resultdata(glass)
  result <- list(call=match.call(),    # function call
                 lca.theta=lca.theta,  # estimates of class conditional response probabilities
                 lca.w=lca.w,          # estimates of mixing proportions
                 prior=prior,          # class prior probabilites
                 m=m,                  # number of subclasses per class
                 r=r,                  # number of possible outcomes per variable
                 k=k,                  # number of classes
                 d=d,                  # number of variables
                 aic=lca.aic,          # AIC
                 bic=lca.bic,          # BIC
                 Gsq=lca.Gsq,          # Gsq
                 Chisq=lca.Chisq)      # Chisq
  class(result) <- "lcda"
  return(result)
}




## print.lcda: readable output


print.lcda <- function(x, ...)
{
    cl <- paste("class", 1:x$k, sep=" ")
    cat("Call:\n")
    print(x$call, ...)
    cat("\nNumber of classes: ", x$k, "\n")
    cat("\nPrior probabilites:\n")
    for (i in 1:x$k)
    cat(paste(cl[i], round(x$prior[i],2), sep=": "), "\n")
    cat("\nNumber of subclasses:\n")
    for (i in 1:x$k)
    cat(paste(cl[i], x$m[i], sep=": "), "\n")
    cat("\nAIC:\n")
    for (i in 1:x$k)
    cat(paste(cl[i], round(x$aic[i], 2), sep=": "), "\n")
    cat("\nBIC:\n")
    for (i in 1:x$k)
    cat(paste(cl[i], round(x$bic[i], 2), sep=": "), "\n")
    cat("\nGsq:\n")
    for (i in 1:x$k)
    cat(paste(cl[i], round(x$Gsq[i], 2), sep=": "), "\n")
    cat("\nChisq:\n")
    for (i in 1:x$k)
    cat(paste(cl[i], round(x$Chisq[i], 2), sep=": "), "\n")
    invisible(x)
}

