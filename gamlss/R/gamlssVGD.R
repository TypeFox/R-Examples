# ------------------------------------------------------------------------------
# this work is a rethinking on Training, Validation and Test data sets analysis
# gamlssVGD() is a combination of the original TGD()  TGD1() TGD2() functions
# it fits a gamlss model at the traning  data (rand==1) and calculates the global
# deviance for the validation set (rand==2)
# TO DO: a) It should also save the residuals for the validation test 
#         Now done the validated residuals are saved as "residVal"
#        b) gamlssCV() needs parallelization
#        c) need a function for comaring CV models    
#------------------------------------------------------------------------------
# functions
#-------------------------------------------------------------------------------
#   i) gamlssVGD() for fiting a model and then calculate the deviance for the extra data
#  ii) VGD() for comparing fitted gamlssVHD models
# iii) getTGD()  after fitting to training data  get the global deviance for the test data
#  iv) TDG() comparing fitted TGD objects
#   v) gamlssCV() for fitting k-folds cross validation
#  vi) CV() comparing fitted CV objects
#-------------------------------------------------------------------------------
#rand <- sample(2, 610, replace=T, prob=c(0.6,0.4)
#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------
# VALIDATION: this is fitting a gamlss model on a sample from the original data 
#  and calculated the Validated Global Deviance from the new Validated data
gamlssVGD <-function(formula = NULL, 
        sigma.formula =~1, 
           nu.formula =~1, 
          tau.formula =~1, 
                 data = NULL, # original data 
               family = NO,  
              control = gamlss.control(trace=FALSE),
                 rand = NULL, # 1 for training 2 for validation
              newdata = NULL, 
                 ...)
 {
#-------------------------------------------------------------------------------
# main functin starts here
#-------------------------------------------------------------------------------
if (is.null(data))   stop("data should be set here")
if (is.null(rand)&&is.null(newdata)) stop("rand or newdata should be set")
if (!is.null(rand))
  {
  if ( any(!rand%in%c(1,2))) stop("rand values should be 1 or 2")
  dataTraining <- subset(data, rand==1)
  dataValid <- subset(data, rand==2)   
  }  
       fname <- as.gamlss.family(family)
        dfun <- paste("d", fname$family[[1]],sep="")
        pfun <- paste("p", fname$family[[1]],sep="")
        lpar <- length(fname$parameters)
        if (is.null(formula)) stop("no formula is set in gamlssVGD")        
        if (is.null(data)) stop("the data argument is needed in gamlssVGD")
 if (!is.null(rand))#  if rand is set to this ----------------------------------
 {
   m1 <- gamlss(formula=formula, sigma.formula=sigma.formula, nu.formula=nu.formula, 
                tau.formula=tau.formula, data=dataTraining, family=family, control=control, ...) 
   nfitted <- predictAll(m1, newdata=dataValid, data=dataTraining)
   if (fname$family[1] %in% .gamlss.bi.list)# if binomial
   {
     if (NCOL(nfitted$y) == 1) 
     {
       y1 <- nfitted$y 
     }
     else                 
     {
       bd <- nfitted$y[,1] + nfitted$y[,2]
       y1 <- nfitted$y[,1]
     }
   } else
   {
     y1 <- nfitted$y 
   }
   if(lpar==1) 
   {
     if (fname$family[[1]] %in% .gamlss.bi.list)
     {
       devi <-  call(dfun, x=y1, mu =  nfitted$mu, bd=bd, log=TRUE) 
       ures <-  call(pfun, q=y1, mu =  nfitted$mu, bd=bd)  
     } else
     {
       devi <-  call(dfun, x=y1, mu =  nfitted$mu,  log=TRUE) 
       ures <-  call(pfun, q=y1, mu =  nfitted$mu)   
     } 
   }
   else if(lpar==2)
   {
     if (fname$family[[1]] %in% .gamlss.bi.list)
     {
       devi <-  call(dfun, x=y1, mu =  nfitted$mu, sigma= nfitted$sigma,  bd=bd, log=TRUE) 
       ures <-  call(pfun, q=y1, mu =  nfitted$nmu, sigma= nfitted$sigma,, bd=bd)  
     } else
     {
       devi <-  call(dfun, x=y1, mu =  nfitted$mu, sigma =  nfitted$sigma, log=TRUE) 
       ures <-  call(pfun, q=y1, mu =  nfitted$mu, sigma =  nfitted$sigma) 
     } 
   }
   else if(lpar==3)
   {
     if (fname$family[[1]] %in% .gamlss.bi.list)
     {
       devi <-  call(dfun, x=y1, mu =  nfitted$mu, sigma =  nfitted$sigma, nu =  nfitted$nu, bd=bd, log=TRUE)
       ures <-  call(pfun, q=y1, mu =  nfitted$mu, sigma =  nfitted$sigma, nu =  nfitted$nu, bd=bd) 
     } else
     {
       devi <-  call(dfun, x=y1, mu =  nfitted$mu, sigma =  nfitted$sigma, nu =  nfitted$nu, log=TRUE)
       ures <-  call(pfun, q=y1, mu =  nfitted$mu, sigma =  nfitted$sigma, nu =  nfitted$nu)
     } 
   }
   else 
   {
     if (fname$family[[1]] %in% .gamlss.bi.list)
     {
       devi <-  call(dfun, x=y1, mu =  nfitted$mu, sigma =  nfitted$sigma, nu =  nfitted$nu, tau= nfitted$tau, bd=bd, log=TRUE)
       ures <-  call(pfun, q=y1, mu =  nfitted$mu, sigma =  nfitted$sigma, nu =  nfitted$nu, tau= nfitted$tau, bd=bd)
     } else
     {
       devi <-  call(dfun, x=y1, mu =  nfitted$mu, sigma =  nfitted$sigma,nu =  nfitted$nu,tau =  nfitted$tau, log=TRUE)
       ures <-  call(pfun, q=y1, mu =  nfitted$mu, sigma =  nfitted$sigma,nu =  nfitted$nu,tau =  nfitted$tau)
     } 
   }
   Vresid <- qNO(eval(ures))
   dev <- -2*sum(eval(devi))
   m1$VGD <- dev
   m1$predictError  <- dev/dim(dataValid)[1]#  end of if rand-------------------
   m1$residVal <-  Vresid 
 } else # if the newdata is set do this-----------------------------------------
 { # 
   m1 <- gamlss(formula = formula, sigma.formula = sigma.formula, nu.formula = nu.formula, 
                tau.formula = tau.formula, data = data, family = family, control = control, 
                ...)
   nfitted <- predictAll(m1, newdata=newdata, data=data)
   if (fname$family[1] %in% .gamlss.bi.list)# if binomial
   {
     if (NCOL(nfitted$y) == 1) 
     {
       y1 <- nfitted$y 
     }
     else                 
     {
       bd <- nfitted$y[,1] + nfitted$y[,2]
       y1 <- nfitted$y[,1]
     }
   } else
   {
     y1 <- nfitted$y 
   }
   if(lpar==1) 
   {
     if (fname$family[[1]] %in% .gamlss.bi.list)
     {
       devi <-  call(dfun, x=y1, mu =  nfitted$mu, bd=bd, log=TRUE) 
       ures <-  call(pfun, q=y1, mu =  nfitted$mu, bd=bd)  
     } else
     {
       devi <-  call(dfun, x=y1, mu =  nfitted$mu,  log=TRUE) 
       ures <-  call(pfun, q=y1, mu =  nfitted$mu)   
     } 
   }
   else if(lpar==2)
   {
     if (fname$family[[1]] %in% .gamlss.bi.list)
     {
       devi <-  call(dfun, x=y1, mu =  nfitted$mu, sigma= nfitted$sigma,  bd=bd, log=TRUE) 
       ures <-  call(pfun, q=y1, mu =  nfitted$nmu, sigma= nfitted$sigma,, bd=bd)  
     } else
     {
       devi <-  call(dfun, x=y1, mu =  nfitted$mu, sigma =  nfitted$sigma, log=TRUE) 
       ures <-  call(pfun, q=y1, mu =  nfitted$mu, sigma =  nfitted$sigma) 
     } 
   }
   else if(lpar==3)
   {
     if (fname$family[[1]] %in% .gamlss.bi.list)
     {
       devi <-  call(dfun, x=y1, mu =  nfitted$mu, sigma =  nfitted$sigma, nu =  nfitted$nu, bd=bd, log=TRUE)
       ures <-  call(pfun, q=y1, mu =  nfitted$mu, sigma =  nfitted$sigma, nu =  nfitted$nu, bd=bd) 
     } else
     {
       devi <-  call(dfun, x=y1, mu =  nfitted$mu, sigma =  nfitted$sigma, nu =  nfitted$nu, log=TRUE)
       ures <-  call(pfun, q=y1, mu =  nfitted$mu, sigma =  nfitted$sigma, nu =  nfitted$nu)
     } 
   }
   else 
   {
     if (fname$family[[1]] %in% .gamlss.bi.list)
     {
       devi <-  call(dfun, x=y1, mu =  nfitted$mu, sigma =  nfitted$sigma, nu =  nfitted$nu, tau= nfitted$tau, bd=bd, log=TRUE)
       ures <-  call(pfun, q=y1, mu =  nfitted$mu, sigma =  nfitted$sigma, nu =  nfitted$nu, tau= nfitted$tau, bd=bd)
     } else
     {
       devi <-  call(dfun, x=y1, mu =  nfitted$mu, sigma =  nfitted$sigma,nu =  nfitted$nu,tau =  nfitted$tau, log=TRUE)
       ures <-  call(pfun, q=y1, mu =  nfitted$mu, sigma =  nfitted$sigma,nu =  nfitted$nu,tau =  nfitted$tau)
     } 
   }
   Vresid <- qNO(eval(ures))
   dev <- -2 * sum(eval(devi))
   m1$VGD <- dev
   m1$predictError  <- dev/dim(newdata)[1]
   m1$residVal <-  Vresid 
   }#----------END of newdata -------------------------------------------------
 class(m1) <- c("gamlssVGD", "gamlss", "gam",    "glm",    "lm"  )
     m1
 } # end of function
#------------------------------------------------------------------------------
#------------------------------------------------------------------------------
# ------------------------------------------------------------------------------
#-------------------------------------------------------------------------------
VGD <- function(object,...) #UseMethod("AIC")
{
  # local function
  is.gamlssVGD <-  function (x)   inherits(x, "gamlssVGD")
  if (length(list(...))) 
  {
      object <- list(object, ...)
    isgamlss <- unlist(lapply(object, is.gamlssVGD))
    if (!any(isgamlss)) stop("some of the objects are not gamlssVGD")
         VGD <- as.numeric(lapply(object, function(x) x$VGD)) 
         val <- cbind(VGD)
        Call <- match.call()
   row.names <- as.character(Call[-1])
       o.val <- order(val)
         val <-  as.data.frame(val[o.val])
       o.r.n <-row.names[o.val]
rownames(val) <- o.r.n 
  val
  }
  else 
  { val <- if (is.gamlssVGD(object)) object$VGD 
    else stop(paste("this is not a gamlssVGD object"))
    val 
  }
}
#-------------------------------------------------------------------------------    
#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------
# this requires a fitted model and the new data only
getTGD<- function (object,   newdata = NULL, ...) 
{
  if (!is.gamlss(object)) stop("not a gamlss object")
  fname <- as.gamlss.family(object$family[[1]])
  dfun <- paste("d", fname$family[[1]], sep = "")
  pfun <- paste("p", fname$family[[1]],sep="")
  lpar <- length(fname$parameters)
  if (is.null(newdata)) 
    stop("no newdata is set in VGD")
  nfitted <- predictAll(object, newdata=newdata, ...)
  if (is.null(nfitted$y)) stop("the response variables is missing in the newdata")
  
  if (fname$family[1] %in% .gamlss.bi.list)# if binomial
  {
    if (NCOL(nfitted$y) == 1) 
    {
      y1 <- nfitted$y 
    }
    else                 
    {
      bd <- nfitted$y[,1] + nfitted$y[,2]
      y1 <- nfitted$y[,1]
    }
  } else
  {
    y1 <- nfitted$y 
  }
  
  if(lpar==1) 
  {
    if (fname$family[[1]] %in% .gamlss.bi.list)
    {
      devi <-  call(dfun, x=y1, mu =  nfitted$mu, bd=bd, log=TRUE) 
      ures <-  call(pfun, q=y1, mu =  nfitted$mu, bd=bd)  
    } else
    {
      devi <-  call(dfun, x=y1, mu =  nfitted$mu,  log=TRUE) 
      ures <-  call(pfun, q=y1, mu =  nfitted$mu)   
    } 
  }
  else if(lpar==2)
  {
    if (fname$family[[1]] %in% .gamlss.bi.list)
    {
      devi <-  call(dfun, x=y1, mu =  nfitted$mu, sigma= nfitted$sigma,  bd=bd, log=TRUE) 
      ures <-  call(pfun, q=y1, mu =  nfitted$nmu, sigma= nfitted$sigma,, bd=bd)  
    } else
    {
      devi <-  call(dfun, x=y1, mu =  nfitted$mu, sigma =  nfitted$sigma, log=TRUE) 
      ures <-  call(pfun, q=y1, mu =  nfitted$mu, sigma =  nfitted$sigma) 
    } 
  }
  else if(lpar==3)
  {
    if (fname$family[[1]] %in% .gamlss.bi.list)
    {
      devi <-  call(dfun, x=y1, mu =  nfitted$mu, sigma =  nfitted$sigma, nu =  nfitted$nu, bd=bd, log=TRUE)
      ures <-  call(pfun, q=y1, mu =  nfitted$mu, sigma =  nfitted$sigma, nu =  nfitted$nu, bd=bd) 
    } else
    {
      devi <-  call(dfun, x=y1, mu =  nfitted$mu, sigma =  nfitted$sigma, nu =  nfitted$nu, log=TRUE)
      ures <-  call(pfun, q=y1, mu =  nfitted$mu, sigma =  nfitted$sigma, nu =  nfitted$nu)
    } 
  }
  else 
  {
    if (fname$family[[1]] %in% .gamlss.bi.list)
    {
      devi <-  call(dfun, x=y1, mu =  nfitted$mu, sigma =  nfitted$sigma, nu =  nfitted$nu, tau= nfitted$tau, bd=bd, log=TRUE)
      ures <-  call(pfun, q=y1, mu =  nfitted$mu, sigma =  nfitted$sigma, nu =  nfitted$nu, tau= nfitted$tau, bd=bd)
    } else
    {
      devi <-  call(dfun, x=y1, mu =  nfitted$mu, sigma =  nfitted$sigma,nu =  nfitted$nu,tau =  nfitted$tau, log=TRUE)
      ures <-  call(pfun, q=y1, mu =  nfitted$mu, sigma =  nfitted$sigma,nu =  nfitted$nu,tau =  nfitted$tau)
    } 
  }
  Vresid <- qNO(eval(ures))
  dev <- -2*sum(eval(devi))
  out <-list()
  out <- list(TGD=dev,  predictError=dev/dim(newdata)[1], resid=Vresid)
  class(out) <- "gamlssTGD"
  out
}
#----------------------------------------------------------------------------------------    
TGD <- function(object,...) #UseMethod("AIC")
{
  # local function
  is.TGD <-  function (x)   inherits(x, "gamlssTGD")
  if (length(list(...))) 
  {
    object <- list(object, ...)
    isTGD <- unlist(lapply(object, is.TGD))
    if (!any(isTGD)) stop("some of the objects are not TGD")
    TGD <- as.numeric(lapply(object, function(x) x$TGD)) 
    val <- cbind(TGD)
    Call <- match.call()
    row.names <- as.character(Call[-1])
    o.val <- order(val)
    val  <-  as.data.frame(val[o.val])
    o.r.n <-row.names[o.val]
    rownames(val) <- o.r.n 
    val
  }
  else 
  { val <- if (is.TGD(object)) object$TGD 
    else stop(paste("this is not a TGD object"))
    val 
  }
}

#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------
# gamlssCV0 <- function(formula = NULL, 
#                      sigma.formula =~1, 
#                         nu.formula =~1, 
#                        tau.formula =~1, 
#                               data = NULL, # original data 
#                             family = NO,  
#                            control = gamlss.control(trace=FALSE),
#                             K.fold = 10,
#                           set.seed = 123,
#                               rand = NULL) 
# {
#   if (is.null(data))   stop("data should be set here")
#         N <- dim(data)[1]
#   #RNAMES <- rownames(data)
#   residCV <- rep(0, N)
#  set.seed <- set.seed
#      rand <- if (is.null(rand)) sample (K.fold , N, replace=TRUE)
#            else rand
#   if (length(rand)!=N) stop("the length of the rand should be equal to data")
#        CV <- rep(0, K.fold)
#   for (i in sort(unique(rand)))
#   {
#     cat("fold ", i, "\n",sep="")
#     learn <- gamlssVGD(formula=formula, sigma.formula=sigma.formula, nu.formula=nu.formula, tau.formula=tau.formula,  family=family, control=control, data=data[rand!=i,], newdata=data[rand==i,])
#     residCV[rand==i] <-  learn$residVal
#                CV[i] <- learn$VGD
#   }
# out <- list(VGD=sum(CV), CVGD=CV, residCV=residCV) 
# class(out) <- "gamlssCV"
# out
# }
#------------------------------------------------------------------------------
gamlssCV <- function(formula = NULL, 
                     sigma.formula =~1, 
                     nu.formula =~1, 
                     tau.formula =~1, 
                     data = NULL, # original data 
                     family = NO,  
                     control = gamlss.control(trace=FALSE),
                     K.fold = 10,
                     set.seed = 123,
                     rand = NULL,
                     parallel = c("no", "multicore", "snow"), 
                     ncpus = 1L, 
                     cl = NULL, ...) 
{
  #--------------- PARALLEL-------------------------------------------------------
  #----------------SET UP PART----------------------------------------------------
  if (missing(parallel)) 
    parallel <- "no"
  parallel <- match.arg(parallel)
  have_mc <- have_snow <- FALSE
  if (parallel != "no" && ncpus > 1L) 
  {
    if (parallel == "multicore") 
      have_mc <- .Platform$OS.type != "windows"
    else if (parallel == "snow") 
      have_snow <- TRUE
    if (!have_mc && !have_snow) 
      ncpus <- 1L
    loadNamespace("parallel")
  }
  # -------------- finish parallel------------------------------------------------
  #-------------------------------------------------------------------------------
  if (is.null(data))   stop("data should be set here")
  N <- dim(data)[1]
  #RNAMES <- rownames(data)
  set.seed <- set.seed
  rand <- if (is.null(rand)) sample(K.fold , N, replace=TRUE)
  else rand
  if (length(rand)!=N) stop("the length of the rand should be equal to data")
  CV <- rep(0, K.fold)
  residCV <- rep(0, N)
  i <- sort(unique(rand))
  #---------------------------------------
  fn <- function(i,...)
  {
    cat("fold ", i, "\n",sep="")
    learn <- gamlssVGD(formula=formula, sigma.formula=sigma.formula, nu.formula=nu.formula, tau.formula=tau.formula,  family=family, control=control, data=data[rand!=i,], newdata=data[rand==i,],...)
    # residCV[rand==i] <<-  learn$residVal
    list(CV=learn$VGD, resid=learn$residVal) 
  }
  #pp<-vapply(i, fn, list("gc", "res"))
  # ll<-lapply(i, fn)
  # ss<-sapply(i, fn)
  #========================================
  # --------  parallel -----------------------------------------------------------
  CV <- if (ncpus > 1L && (have_mc || have_snow)) 
  { 
    if (have_mc) 
    {# sapply(scope, fn)
      parallel::mclapply(i, fn, mc.cores = ncpus)
    }
    else if (have_snow) 
    {
      list(...)
      if (is.null(cl)) 
      {
        # make the cluster
        # cl <- parallel::makePSOCKcluster(rep("localhost", ncpus))
        cl <- parallel::makeForkCluster(ncpus)
        if (RNGkind()[1L] == "L'Ecuyer-CMRG") 
          parallel::clusterSetRNGStream(cl)
        res <- parallel::parLapply(cl, i, fn)
        parallel::stopCluster(cl)
        res
      }
      else t(parallel::parLapply(cl, i, fn))
    }
  } # end parallel ---------------------------------------------------------------
  else lapply(i, fn)
  cv <- rep(0, K.fold)
  for (i in 1:K.fold)
  {
    residCV[rand==i] <-  CV[[i]]$resid
    cv[i] <- CV[[i]]$CV
  }
  # CV <-  sapply(i, fn) 
  out <- list(CV=sum(cv), allCV=cv, residCV=residCV) 
  class(out) <- "gamlssCV"
  out
}
#-------------------------------------------------------------------------------
CV <- function(object,...) #UseMethod("AIC")
{
  # local function
  is.CV <-  function (x)   inherits(x, "gamlssCV")
  if (length(list(...))) 
  {
    object <- list(object, ...)
    isCV <- unlist(lapply(object, is.CV))
    if (!any(isCV)) stop("some of the objects are not gamlssCV")
    CV <- as.numeric(lapply(object, function(x) x$CV)) 
    val <- cbind(CV)
    Call <- match.call()
    row.names <- as.character(Call[-1])
    o.val <- order(val)
    val  <-  as.data.frame(val[o.val])
    o.r.n <- row.names[o.val]
    rownames(val) <- o.r.n 
    val
  }
  else 
  { val <- if (is.CV(object)) object$CV 
    else stop(paste("this is not a gamlssCV object"))
    val 
  }
}
