dataXY <- function(formula, data,
                         family = "gaussian",
                         trialscol = NULL,
                         trans.power = NULL,
                         trans.shift = 0,
                         distord,
                         CorModels
                         )
{
  ##start function here, pass distord to function
    ## data.orig <- data ## not used again so commented out
  ##make sure the order of the data matches order in distance matrix
  data <- data[distord,]

  # get a list of response and covariate names
  mod.names <- as.character(attr(terms(formula, data = data),"variables"))

  tt <- try(model.matrix(formula,data))
  if(class(tt)=="try-error") stop("The specified formula contains invalid variable names\n",tt)

  if(!all(mod.names[-1] %in% colnames(data))) stop("The specified formula causes problems in the dataXY function, try to use a simple formula, avoid indicator functions if you can.")

  ## Check random and fixed effects are different
  if(any(mod.names[-1] %in% CorModels)) stop("Random effects and fixed effects overlap")

  ## get the number of names ( + 1, as the first is always "list")
  nc.tmp <- length(mod.names)
  # name of the response variable
  response.col <- mod.names[2]

  # total number of observations
  n.all <- length(data[,1])
  # create a vector of all TRUE values
  ind.all <- rep(TRUE, times = n.all)
  ind.allcov <- ind.all
  # if there are any covariates ...
  if(nc.tmp > 2) {
    # create a FALSE for a record with missing values of the covariates
    for(i in 3:nc.tmp) ind.allcov <- ind.allcov & !is.na(data[,mod.names[i]])
  }
  ## Check the random effects also
  REind <- which(names(data) %in% CorModels)
  if(length(REind)) {
      for(ii in REind) ind.allcov <- ind.allcov & !is.na(data[,ii])
  }

  # sample size without missing covariates
  n.allcov <- sum(ind.allcov)
  # remove records that had any missing values for any of the covariates
  data1 <- data[ind.allcov,]
  #standardize numeric covariates and turn character fields to factors
  covnames <- NULL
  if(nc.tmp > 2) covnames <- mod.names[3:nc.tmp]
#  StdXDataFrame <- NULL
  if(!is.null(covnames)) {
    for(i in 1:length(covnames)) {
#      if(is.numeric(data1[,covnames[i]])) {
#        xmean <- mean(data1[,covnames[i]])
#        xsdev <- sqrt(var(data1[,covnames[i]]))
#        data1[,covnames[i]] <- (data1[,covnames[i]] - xmean)/xsdev
#        StdXDataFrame <- rbind(StdXDataFrame,
#                               data.frame(VariableName = covnames[i], Mean = xmean, StdDev = xsdev))
#      }
      if(is.character(data1[,covnames[i]])) {
        data1[,covnames[i]] <- as.factor(data1[,covnames[i]])
      }
    }
  }
  # replace response column with all 1's for design matrix of all records
  z1 <- data1[,response.col]
  data1[,response.col] <- rep(1, times = n.allcov)
  mf <- model.frame(formula, data = data1)
  mt <- attr(mf, "terms")
  X1 <- model.matrix(mt, mf, contrasts)
  # get names for all terms and interactions, including those set to zero
  terms.list <- attr(mt,"term.labels")
  if(attr(mt,"intercept") == 1) effnames <- "(Intercept)"
  else effnames <- "NULL"
  if(!is.null(covnames)) {
    for (i in 1:length(terms.list)) {
      form1 <- formula(paste("~ ", terms.list[i], " - 1"))
      Xt <- model.matrix(form1, data = data1)
      effnames <- c(effnames, colnames(Xt))
    }
  }
  setzero <-  !(effnames %in% colnames(X1))
  # check for estimability
  X1re <- rref(X1)
  # create X1 with estimable functions
  cutX1toX2 <- apply(abs(X1re),2,sum) == 1
  if(any(cutX1toX2 == FALSE)) X1 <- X1[,cutX1toX2]
  setNA <- !(effnames %in% colnames(X1)) & !setzero

  data1[,response.col] <- z1

  #indicator vector of all records without missing covariates or response
  ind.allxy <- ind.allcov & !is.na(data[,response.col])
  #sample size of all records without missing covariates or response
  n.allxy <- sum(ind.allxy)
  # indicator vector of data1 without missing covariates
  ind.ysubx <- !is.na(data1[response.col])
  # data set of all records without missing covariates or response
  data2 <- data[ind.allxy,]
  # create a working column of the response variable
  data2[,"work.col"] <- data2[,response.col]
  # transformations on the working column
  if(!is.null(trans.power)) {
    if(family != "gaussian")
      return(print("Power transformation can only be used with gaussian family"))
    if(trans.power < 0) return(print("Power transformation must be > 0"))
    if(trans.power > .5) return(print("Power transformation must be < 0.5"))
    if(trans.power == 0) {
	  if(any((data2[,response.col] + trans.shift) <= 0))
		  return(print("Data must be > 0 to use log transformation"))
      data2[,"work.col"] <- log(data2[,response.col] + trans.shift)
    } else
    data2[,"work.col"] <- (data[,response.col] + trans.shift)^trans.power
  }
  X1 <- as.matrix(X1)
  X2 <- X1[ind.ysubx,]
  X2 <- as.matrix(X2)
  # check for estimability
  X2re <- rref(X2)
  # create X1 with estimable functions
  cutX1toX2 <- apply(abs(X2re),2,sum) == 1
  if(any(cutX1toX2 == FALSE)) X2 <- X2[,cutX1toX2]
  setNA2 <- !(effnames %in% colnames(X2))
  # data set of observed data only
  if(is.factor(data2[,"work.col"]))
		data2[,"work.col"] <- as.numeric(as.character(data2[,"work.col"]))
  if(is.character(data2[,"work.col"]))
		data2[,"work.col"] <- as.numeric(data2[,"work.col"])
  # vector of observed response variable
  z <- as.matrix(data2[,"work.col"], ncol = 1)
	attr(z,"pid") <- data2[,"pid"]
  trialsvec <- NULL
  # if response variable is binomial with n trials change z to proportion
  if(!is.null(trialscol) & family == "binomial"){
    if(is.factor(data2[,"trialscol"]))
	  data2[,"trialscol"] <- as.numeric(as.character(data2[,"trialscol"]))
	if(is.character(data2[,"trialscol"]))
	  data2[,"trialscol"] <- as.numeric(data2[,"trialscol"])
	trialsvec <- data2[,"trialscol"]
	z <- as.matrix(z/trialsvec, ncol = 1)
  }
  # else if response is Bernoulli, set trialsvec to all ones
  if(is.null(trialscol) & family == "binomial")
    trialsvec <- rep(1, times = n.allxy)

  # if any missing response values, design matrix of missing records
  data.na <- NULL
  ind.RespNA <- NULL
  if((n.allcov - n.allxy) > 0) {
	data.na <- data1[!ind.ysubx,]
	ind.RespNA <- data[,"pid"] %in% data.na[,"pid"]
  }
  # get the rank of the design matrix
  p <- sum(svd(X1)$d>1e-10)

  ## GET REs
  REs <- NULL
  REmodelmatrices <- NULL
  REind <- which(names(data) %in% CorModels)
  if(length(REind)) {
      REs <- list()
      REnames <- sort(names(data)[REind])
      ## model matrix for a RE factor
      for(ii in 1:length(REind)) REmodelmatrices[[REnames[ii]]] <- model.matrix(~data[ind.allxy,REnames[ii]] - 1)
      ## corresponding block matrix
      for(ii in 1:length(REind)) REs[[ii]] <- tcrossprod(REmodelmatrices[[ii]])
      names(REs) <- REnames
      names(REmodelmatrices) <- REnames
  }

  outpt <- list(
                datasets = list(
                data0 = data,
                data1 = data1,
                data2 = data2,
                data.na = data.na
                ),
                indvecs = list(
                ind.allcov = ind.allcov,
                ind.allxy = ind.allxy,
                ind.ysubx = ind.ysubx,
                ind.RespNA = ind.RespNA,
                setzero = setzero,
                setNA = setNA,
                setNA2 = setNA2,
                cutX1toX2 = cutX1toX2
                ),
                sampsizes = list(
                n.all = n.all,
                n.allcov = n.allcov,
                n.allxy = n.allxy
                ),
                Xmats = list(
                X1 = X1,
                X2 = X2,
                p = p,
                ##	  StdXDataFrame = StdXDataFrame,
                effnames = effnames
                ),
                respvecs = list(
                response.col = response.col,
                z = as.matrix(z, ncol = 1),
                trialsvec = trialsvec
                ),
                REs = REs,
                REmodelmatrices = REmodelmatrices
                )
  outpt
}

