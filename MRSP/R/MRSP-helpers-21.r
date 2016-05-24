################################################################################
#####    Helper functions for MRSP.                                        #####
################################################################################
#####    Author: Wolfgang Pößnecker                                        #####
#####    Last modified: 21.10.2014, 17:40                                  #####
################################################################################


## "evalh(ere)": evaluate a function "fun" in the environment where evalh is used,
## with the value of the function arguments as found in the current environment
evalh <- function(fun, funargs = NULL){
 if(is.character(fun)) fun <- get(fun)
 if(!is.function(fun)) stop("fun is not a function or a valid function name in evalh")
 environment(fun) <- parent.frame()
 if(is.null(formals(fun))){
  fun()
 }else{
  if(is.null(funargs)) funargs <- mget(names(formals(fun)), envir = parent.frame())
  eval(as.call(c(fun, funargs)))
 }
}

## a function to extract a named element from each part of a list. name must be
## the quoted (!) name of the element to be extracted.
lextract <- function(object, name, simplify=FALSE, recursive=FALSE){
 out <- lapply(object, "[[", as.character(name))
 if(simplify) out <- unlist(out, recursive=recursive)
 out
}

## a vectorized version of all.equal
v.allequal <- function(x, y, tolerance = .Machine$double.eps ^ 0.5){
 if(length(x) != length(y) && !(length(x) == 1 | length(y) == 1))
  stop("length of x and y in v.allequal do not match")
 x > y - tolerance & x < y + tolerance
}


## a function for a more efficient computation of the column-wise mean in 
## treshsymx and similar functions. Background: this function is called many
## thousand times during a typical call to MRSP, but it has to deal both with
## vector arguments, which are internally represented as matrices with one
## column in MRSP (to allow one code to deal with both binomial and multinomial
## responses) and with "true matrix arguments". However, due to the thousands of
## calls to these functions, it is much quicker to use "sum/length" instead of 
## "mean" because this avoids overhead. for true matrix arguments though, there
## is no way to avoid the use of "colMeans". 
treshmean <- function(x){if(ncol(x) == 1) sum(x)/length(x) else colMeans(x)}

## Function for evaluating rowsums of logarithms in a numerically more stable way
## Input: matrix m of "log'd" values
log_row_sums <- function(m) {
  M <- maxRow(m)                 # maxRow is a function provided by package MRSP
  M + log(rowSums(exp(m - M)))
}

## a function for a numerically stable computation of log(1 + exp(x)), according
## to Maechler (2012)
log1pexp <- cmpfun(function(x){
 ind1 <- which(x <= 18)
 ind2 <- which(x > 18 && x <= 33.3)
 #ind3 <- which(x > 33.3)
 out <- x
 out[ind1] <- log1p(exp(x[ind1]))
 out[ind2] <- x[ind2] + exp(-x[ind2])
 #out[ind3] <- x[ind3]
 out
}) 
 

## a function for computing the number of knots used in the call to smooth.spline
## in 'plot'.
n.knots <- function(n) {
    ## Number of inner knots
    if(n < 50L) n
    else trunc({
        a1 <- log2( 50)
        a2 <- log2(100)
        a3 <- log2(140)
        a4 <- log2(200)
        if	(n < 200L) 2^(a1+(a2-a1)*(n-50)/150)
        else if (n < 800L) 2^(a2+(a3-a2)*(n-200)/600)
        else if (n < 3200L)2^(a3+(a4-a3)*(n-800)/2400)
        else  200 + (n-3200)^0.2
    })
}

## a function for computing the unpenalized ML estimate
## note that a small ridge penalty will automatically be added if the program 
## detects divergence of parameter estimates, so this "ML"-fit will always exist
getmlfit <- function(dat, coef.init = NULL, coef.stand.init = NULL, coef.pretres.init = NULL,
                          offset = rep(0, nrow(dat$y)), weights = rep(1, nrow(dat$y)), grpindex = NULL, penindex = NULL,
                          lambda = 123, lambdaR = lambda, lambdaF = lambda, gamma = 1, psi = 1,
                          indg = NULL, indcs = NULL, model = multinomlogit(),
                          constr = NULL, control = MRSP.control(), fista.control = NULL, Proximal.control = NULL,
                          Proximal.args = NULL, adaptive = FALSE, threshold = FALSE, refit = FALSE, nonneg = FALSE, ...)
{
 ## if any objects that are not formals of MRSP.fit are passed via ... , we need
 ## to explicitly assign the corresponding objects:
 dotlist <- list(...)
 for(i in seq_len(length(dotlist))) assign(names(dotlist)[i], dotlist[[i]])
 if(is.expression(model)) model <- eval(model)

 if(missing(constr) | is.null(constr))  constr <- model@constraint
 if(nonneg == T && constr == "symmetric"){
  constr <- 1
  warning()
  cat(" You cannot combine a nonnegativity constraint on the coefficients with a symmetric side constraint!", "\n", "The first response category is instead used as reference.","\n")
 }
 hasV <- !is.null(dat$V)
 P <- ncol(dat$x)
 if(hasV) L <- ncol(dat$V[[1]])
   
 if(missing(grpindex) | is.null(grpindex)){
  grpindex <- list()
  grpindex[[1]] <- seq_len(P)
  if(hasV){
   grpindex[[2]] <- seq(from = P+1, to = P+L)
  } 
 } 

 #lambda <- list(12345)
 lambda <- 12345
 #lambdaR <- list(23456)
 lambdaR <- 23456
 lambdaF <- 34567
 
 penweights <- list()
 penweights[[1]] <- list()
 penweights[[2]] <- list()
 
 penweights[[1]][[1]] <- rep(1, length(unique(grpindex[[1]])))
 if(hasV){
  penweights[[1]][[2]] <- rep(1, length(unique(grpindex[[2]])))
 }
 penweights[[2]][[1]] <- matrix(1, nrow = ncol(dat$y), ncol = length(unique(grpindex[[1]])))
 if(hasV){
  penweights[[2]][[2]] <- matrix(1, nrow = ncol(dat$y), ncol = length(unique(grpindex[[2]])))
 }
 
 if(any(penindex[[1]] %in% c(4, 40, 41))){
  sl.indg <- which(penindex[[1]] %in% c(4, 40, 41))
 }else{sl.indg <- NULL}
 #if(any(penindex[[1]] %in% c(15))){
 # which.Pspline <- which(penindex[[1]] == 15)
 # penindex <- list()
 # penindex[[1]] <- rep(10, P)
 # penindex[[1]][which.Pspline] <- 13
 # lambdaR <- control@lambdastabil
 #}else{
  penindex <- list()
  penindex[[1]] <- rep(10, P)
 #}
 penindex[[1]][sl.indg] <- 40
 if(hasV){
  penindex[[2]] <- numeric(L)
  penindex[[2]][indg] <- 20
  penindex[[2]][indcs] <- 30  
 } 

 control@doMLfit <- F

 if("scale.notpen.ML" %in% names(dotlist)){              ## these lines are necessary to outsource standardization of covariates from MRSP.fit
  if("scale.notpen.x" %in% names(dotlist)){              ## to MRSP while being able to backtransform coefs to the original scale.
   dotlist$scale.notpen.x <- dotlist$scale.notpen.ML
   dotlist$scale.notpen.ML <- NULL
  }
 }
 if("scale.notpen.VML" %in% names(dotlist)){
  dotlist$scale.notpen.V <- dotlist$scale.notpen.VML
  dotlist$scale.notpen.VML <- NULL
 }

 initialfit <- eval(as.call(c(as.symbol("MRSP.fit"), list(dat = dat, coef.init = coef.init, coef.stand.init = coef.stand.init,
                         coef.pretres.init = coef.pretres.init, offset = offset, weights = weights,
                         grpindex = grpindex, penindex = penindex, lambda = lambda, lambdaR = lambdaR,
                         lambdaF = lambdaF, gamma = gamma, psi = psi, indg = indg, indcs = indcs, model = model,
                         constr = constr, control = control, fista.control = fista.control,
                         Proximal.control = Proximal.control, Proximal.args = Proximal.args,
                         penweights = penweights, adaptive = FALSE, threshold = FALSE, refit = FALSE, fusion = FALSE,
                         nonneg = nonneg), dotlist)))
 
 out <- list(coef = initialfit@coef,
             coef.stand = initialfit@coef.stand,
             coef.pretres = initialfit@coef.pretres,
             penweights = initialfit@penweights)
 
 return(out)
}             



## a function for computing the weights if adaptive lasso type of penalties are
## used.
getpenweights <- function(dat, coef.init = NULL, coef.stand.init = NULL, coef.pretres.init = NULL,
                          offset = rep(0, nrow(dat$y)), weights = rep(1, nrow(dat$y)), grpindex = NULL, penindex = NULL,
                          lambda = list(123), lambdaR = lambda, lambdaF = lambda, gamma = 1, psi = 1, indg = NULL,
                          indcs = NULL, model = multinomlogit(),
                          constr = NULL, control = MRSP.control(), fista.control = NULL, Proximal.control = NULL,
                          Proximal.args = NULL, adaptive = FALSE, threshold = FALSE, refit = FALSE, mlfit, fusion = FALSE,
                          nonneg = nonneg, ...)
{
 if(is.expression(model)) model <- eval(model)
 if(missing(constr) | is.null(constr))  constr <- model@constraint
 if(nonneg == T && constr == "symmetric"){
  constr <- 1
  warning()
  cat(" You cannot combine a nonnegativity constraint on the coefficients with a symmetric side constraint!", "\n", "The first response category is instead used as reference.","\n")
 }
 hasV <- !is.null(dat$V)
 P <- ncol(dat$x)
 K <- ncol(dat$y)
 if(hasV) L <- ncol(dat$V[[1]])
 if(any(penindex[[1]] %in% c(4, 40, 41))){
  sl.indg <- which(penindex[[1]] %in% c(4, 40, 41))
 }else{sl.indg <- NULL}

 ## set up the difference matrix 'R' for GFL-fusion-type penalties
 if(fusion == F){
  R <- NULL
 }else{
  if(is.numeric(constr)) Kf <- K-1 else Kf <- K         ## Kf = relevant size of the problem for fusion penalties
  #if(Kf < 2) stop("a fusion penalty was requested for a too small parameter group")
  if(Kf <= 2) R <- matrix(c(-1,1), nrow=1)
  if(Kf > 2){
   if(fusion == "adja"){
    R <- c(-1, 1, rep(0, Kf-2))
    for(i in 1:(Kf-2)){
     R <- rbind(R, c(rep(0,i),-1,1,rep(0,Kf-2-i)))
    }
    R <- matrix(R, ncol=Kf)
   }else if(fusion == "all"){
    R <- cbind(rep(1,Kf-1), -diag(Kf-1))
    for(i in 1:(Kf-1)){
     add <- cbind(matrix(0,nrow=Kf-1-i,ncol=i),rep(1,Kf-1-i),-diag(Kf-1-i))
     R <- rbind(R,add)
    }
   }
  }
 }#else stop("option 'fusion' in MRSP.fit must be either FALSE or ``adja`` or ``all``")
   
 if(missing(grpindex) | is.null(grpindex)){
  grpindex <- list()
  grpindex[[1]] <- seq_len(P)
  if(hasV){
   grpindex[[2]] <- seq(from = P+1, to = P+L)
  } 
 } 
  
 grpindexunique <- lapply(grpindex, unique)

 if(missing(mlfit) | is.null(mlfit)){
  penweights <- list()
  penweights[[1]] <- list()
  penweights[[2]] <- list()
  penweights[[3]] <- list()
  penweights[[4]] <- list()
  penweights[[5]] <- list()
  penweights[[1]][[1]] <- rep(1, length(grpindexunique[[1]])) # for CATS-type penalties
  penweights[[2]][[1]] <- matrix(1, nrow = ncol(dat$y), ncol = length(grpindexunique[[1]])) # for ordinary lasso penalties
  if(fusion != F && ( (!is.numeric(constr) & K > 1) || (is.numeric(constr) & K > 2))){
   penweights[[3]][[1]] <- rep(1, length(grpindexunique[[1]])) # for grouped fusion penalties
   penweights[[4]][[1]] <- matrix(1, nrow = nrow(R), ncol = length(grpindexunique[[1]])) # for L1-type fusion penalties
  }
  if(any(penindex[[1]] == 18)){
   penweights[[5]][[1]] <- list()                                                        # for row-L1-type fusion penalties
   for(j in seq(length(grpindexunique[[1]]))){
    j.which <- which(grpindex[[1]] == j)
    if(fusion == "adja"){
     nrow5x <- length(j.which) - 1
    }else if(fusion == "all"){
     nrow5x <- choose(length(j.which), 2)
    }
    penweights[[5]][[1]][[j]] <- matrix(1, nrow = nrow5x, ncol = K)
   }
  }
  if(hasV){
   penweights[[1]][[2]] <- rep(1, length(grpindexunique[[2]]))
   penweights[[2]][[2]] <- matrix(1, nrow = ncol(dat$y), ncol = length(grpindexunique[[2]]))
  }  
 }else{
 #l2normraw <- function(u){
 # if(length(c(u)) > 1) sqrt(sum(c(u)^2)) #* sqrt(length(c(u))) #* sqrt(1 + var(abs(c(u)))/max(mean(abs(c(u)))^2, 1e-3))
 # else sqrt(sum(c(u)^2))
 #}
  
 varadjust <- function(u){
  1 #if(length(c(u)) > 1) sqrt(1 + var(abs(c(u)))/max(mean(abs(c(u)))^2, 1e-3)) else 1
 } 
 
 penweights <- list()
 penweights[[1]] <- list()
 penweights[[2]] <- list()
 penweights[[3]] <- list()
 penweights[[4]] <- list()
 penweights[[5]] <- list()
 penweights[[1]][[1]] <- numeric(length(grpindexunique[[1]]))
 penweights[[2]][[1]] <- matrix(nrow = nrow(mlfit$coef.stand[[1]]), ncol = length(grpindexunique[[1]]))
 if(fusion != F && ( (!is.numeric(constr) & K > 1) || (is.numeric(constr) & K > 2) )){
  penweights[[3]][[1]] <- numeric(length(grpindexunique[[1]]))
  penweights[[4]][[1]] <- matrix(1, nrow = nrow(R), ncol = length(grpindexunique[[1]]))
 }
 if(any(penindex[[1]] == 18)){
  penweights[[5]][[1]] <- list()
  for(j in seq(length(grpindexunique[[1]]))){
   j.which <- which(grpindex[[1]] == j)
   if(fusion == "adja"){
    nrow5x <- length(j.which) - 1
   }else if(fusion == "all"){
    nrow5x <- choose(length(j.which), 2)
   }
   penweights[[5]][[1]][[j]] <- matrix(1, nrow = nrow5x, ncol = K)
  }
 }
 if(hasV){
  penweights[[1]][[2]] <- numeric(length(grpindexunique[[2]]))
  penweights[[2]][[2]] <- matrix(nrow = nrow(mlfit$coef.stand[[2]]), ncol = length(grpindexunique[[2]]))
 }
 
 for(j in grpindexunique[[1]]){
  j.which <- which(grpindex[[1]] == j)
  penweights[[1]][[1]][j] <- varadjust(mlfit$coef.stand[[1]][,j.which])
  penweights[[2]][[1]][,j] <- t(apply(mlfit$coef.stand[[1]][,j.which, drop=F], 1, varadjust))
 if(fusion != F && ( (!is.numeric(constr) & K > 1) || (is.numeric(constr) & K > 2) )){
   penweights[[3]][[1]][j] <- varadjust(mlfit$coef.stand[[1]][,j.which])
   penweights[[4]][[1]][,j] <- t(apply(R%*%mlfit$coef.stand[[1]][,j.which, drop=F], 1, varadjust))
  }
  if(!is.null(sl.indg)){
   if(all(j.which %in% sl.indg)){
    penweights[[1]][[1]][j] <- varadjust(mlfit$coef.stand[[1]][1,j.which])
   }
  } 
 }
 if(hasV){
  for(j in grpindexunique[[2]]){
   j.which <- which(grpindex[[2]] == j)
   if(all(j.which %in% indg)){
    penweights[[1]][[2]][j-length(grpindexunique[[1]])] <- varadjust(mlfit$coef.stand[[2]][1,j.which])
   }else if(all(j.which %in% indcs)){
    penweights[[1]][[2]][j-length(grpindexunique[[1]])] <- varadjust(mlfit$coef.stand[[2]][,j.which])   
   }else{
    stop("something went wrong with j.which, indg and indcs in the computation of penweights")
   }  
   penweights[[2]][[2]][,j-length(grpindexunique[[1]])] <- t(apply(mlfit$coef.stand[[2]][,j.which, drop=F], 1, varadjust))
  } 
 }

 if(adaptive != F){
  adaexpo <- control@adapt.expo
  control@keepdat <- T
  
  if(adaptive != "ML"){
   if(length(lambda) == 1){
    nrlambda <- control@nrlambda
    lambdabase <- control@lambdabase
    lambdamax <- getlambdamax(dat = dat, offset = offset, weights = weights, 
                              grpindex = grpindex, penindex = penindex, psi = psi,
                              gamma = gamma, model = model, constr = constr, control = control,
                              penweights = penweights, indg = NULL, indcs = NULL, adaptive = FALSE,
                              threshold = threshold, refit = refit, lambdaF = lambdaF, nonneg = nonneg)
    lambdamin <- min(0.05, lambdamax / 100)
    lambdalogvals <- seq(log(lambdamin / lambdabase), log(lambdamax / lambdabase), length.out = nrlambda)
    lambda <- lambdabase * exp(lambdalogvals)
    lambda <- list(lambda)
   } 
  
   if(length(lambdaR) != length(lambda)) lambdaR <- lambda
  
   initialfit <- eval(as.call(c(as.symbol("MRSP.fit"), list(dat = dat, coef.init = coef.init, coef.stand.init = coef.stand.init,
                           coef.pretres.init = coef.pretres.init, offset = offset, weights = weights, 
                           grpindex = grpindex, penindex = penindex, lambda = lambda, lambdaR = lambdaR,
                           lambdaF = lambdaF, gamma = gamma, psi = psi, indg = indg, indcs = indcs, model = model,
                           constr = constr, control = control, fista.control = fista.control, 
                           Proximal.control = Proximal.control, Proximal.args = Proximal.args, 
                           penweights = penweights, mlfit = mlfit, adaptive = FALSE, threshold = FALSE,
                           refit = FALSE, fusion = fusion, nonneg = nonneg), list(...))))
  
   cvinit <- cv(initialfit, k = 5, parallel = control@adaptcv)
   bestm <- which.min(cvinit$mean)
   initialfit <- initialfit[[bestm]]
   initialfit <- list(coef = initialfit@coef,
                      coef.stand = initialfit@coef.stand,
                      coef.pretres = initialfit@coef.pretres)
  }else{initialfit <- mlfit}
    
  lp <- length(penweights[[1]][[1]])
  for(j in grpindexunique[[1]]){
   j.which <- which(grpindex[[1]] == j)
   penweights[[1]][[1]][j] <- penweights[[1]][[1]][j] / max(l2normraw(initialfit$coef.stand[[1]][,j.which])^adaexpo, 1e-10)
   penweights[[2]][[1]][,j] <-  penweights[[2]][[1]][,j] * c(apply(initialfit$coef.stand[[1]][,j.which, drop=F], 1,
                                 function(u){1 / max(l2normraw(u)^adaexpo, 1e-10)}))
   if(fusion != F && ( (!is.numeric(constr) & K > 1) || (is.numeric(constr) & K > 2))){
    penweights[[3]][[1]][j] <- penweights[[3]][[1]][j] / max(l2normraw(R%*%initialfit$coef.stand[[1]][,j.which])^adaexpo, 1e-10)
    penweights[[4]][[1]][,j] <-  penweights[[4]][[1]][,j] / max(rowL2normraw(R%*%initialfit$coef.stand[[1]][,j.which, drop=F])^adaexpo, 1e-10) 
   }
   if(any(penindex[[1]][j.which] == 18)){ ## first we have to setup the difference matrix D for this particular group
    kj <- length(j.which)
    if(kj < 2) stop("a fusion penalty was requested for a too small parameter group")
    if(kj == 2){
     D <- c(-1,1)
    }else if(kj > 2){
     if(fusion == "adja"){
      D <- c(-1, 1, rep(0, kj-2))
      for(i in 1:(kj-2)){
       D <- rbind(D, c(rep(0,i),-1,1,rep(0,kj-2-i)))
      }
      D <- matrix(D, ncol=kj)
     }else if(fusion == "all"){
      D <- cbind(rep(1, kj-1), -diag(kj-1))
      for(i in 1:(kj-1)){
       add <- cbind(matrix(0, nrow=kj-1-i, ncol=i), rep(1, kj-1-i), -diag(kj-1-i))
       D <- rbind(D, add)
      }
     }
    }  
    penweights[[5]][[1]][[j]] <- penweights[[5]][[1]][[j]] * apply(D%*%t(initialfit$coef.stand[[1]][,j.which, drop=F]), c(1,2), 
                                                                   function(u){ 1 / max(abs(u)^adaexpo, 1e-10)})
   } 
   if(!is.null(sl.indg)){
    if(all(j.which %in% sl.indg)){
     penweights[[1]][[1]][j] <- penweights[[1]][[1]][j] / max(l2normraw(initialfit$coef.stand[[1]][1,j.which])^adaexpo, 1e-10)
    }
   }  
  }
  if(hasV){
   for(j in grpindexunique[[2]]){
    j.which <- which(grpindex[[2]] == j)    
    if(all(j.which %in% indg)){
     penweights[[1]][[2]][j - lp] <- penweights[[1]][[2]][j - lp] / max(l2normraw(initialfit$coef.stand[[2]][1,j.which])^adaexpo, 1e-10)
    }else if(all(j.which %in% indcs)){
     penweights[[1]][[2]][j - lp] <- penweights[[1]][[2]][j - lp] / max(l2normraw(initialfit$coef.stand[[2]][,j.which])^adaexpo, 1e-10) 
    }else{
     stop("something went wrong with j.which, indg and indcs in the computation of penweights")
    }      
    penweights[[2]][[2]][,j - lp] <- penweights[[2]][[2]][,j - lp] * c(apply(initialfit$coef.stand[[2]][,j.which, drop=F], 1,
                                 function(u){1 / max(l2normraw(u)^adaexpo, 1e-10)}))
   }
  }
 }
 
 } # end if(missing(mlfit)....
 return(penweights)
} 

## a function for computing the maximal lambda for a given model. this maxlambda
## leads to a model where only the intercepts are left.
## much of the code is the same as in MRSP.fit for a single lambda.
getlambdamax <- function(dat, offset = rep(0, nrow(dat$y)), weights = rep(1, nrow(dat$y)),
                         grpindex, penindex, psi = 1, gamma = 1, fista.control = NULL,
                         model = multinomlogit(), constr = NULL, control = MRSP.control(),
                         penweights, indg = NULL, indcs = NULL, adaptive = FALSE,
                         threshold = FALSE, refit = FALSE, fusion = FALSE, lambdaF = 0, nonneg = FALSE)
{
 if(is.expression(model)) model <- eval(model)
 if(missing(constr) | is.null(constr))  constr <- model@constraint
 if(nonneg == T && constr == "symmetric"){
  constr <- 1
  warning()
  cat(" You cannot combine a nonnegativity constraint on the coefficients with a symmetric side constraint!", "\n", "The first response category is instead used as reference.","\n")
 }
 ## are there any category-specific predictors? they are stored in the 'V'-entry of dat
 hasV <- !is.null(dat$V)
 
 nobs <- nrow(dat$y)
 wnobs <- sum(weights)
 K <- ncol(dat$y)
 P <- ncol(dat$x)
 if(hasV & !is.matrix(dat$V[[1]])){dat$V <- lapply(dat$V, as.matrix)}
 if(hasV){L <- ncol(dat$V[[1]])}
 
 if(hasV){
  if(missing(penindex) | is.null(penindex)){
   if((missing(indg) | is.null(indg)) && (missing(indcs) | is.null(indcs))){
    indg <- integer(0)
    indcs <- seq(1, L)
   }else if((missing(indg) | is.null(indg)) && !(missing(indcs) | is.null(indcs))){
    indg <- seq(1, L)[-which(seq(1, L) %in% indcs)]
   }else if(!(missing(indg) | is.null(indg)) && (missing(indcs) | is.null(indcs))){
    indcs <- seq(1, L)[-which(seq(1, L) %in% indg)]
   }  
  }else{
   indg <- which(penindex[[2]] %in% c(2, 20, 21))
   indcs <- which(penindex[[2]] %in% c(3, 30, 31, 32, 33))
  } 
 }
  

 ## handling the intercept
 which.intercept <- which(apply(dat$x, 2, function(u) diff(range(u)) < .Machine$double.eps ^ 0.5))
 has.intercept <- !(length(which.intercept) == 0)
 if(!has.intercept)
   stop("MRSP currently cannot compute a suitable maximal lambda on its own if the model is specified without an intercept!",'\n',
        "  Please supply argument 'lambdamax' by hand or add an intercept.",'\n')
 if((length(which.intercept) == 1) && which.intercept != 1){
  dat$x[,c(1, which.intercept)] <- dat$x[,c(which.intercept, 1)] 
  penindex[[1]][c(1,which.intercept)] <- penindex[[1]][c(which.intercept,1)]
  grpindex[[1]][c(1,which.intercept)] <- grpindex[[1]][c(which.intercept,1)]
  if(!all(coef[[1]] == 0))
  coef[[1]][,c(1,which.intercept)] <- coef[[1]][,c(which.intercept,1)] 
 }
 if(length(which.intercept) > 1)
   stop("more than one intercept column specified")

 ## are we using an ordinal regression model?
 isordinal <- model@name %in% c("Cumulative Logit Model", "Sequential Logit Model")
 if(isordinal){
  sl.indg <- which(penindex[[1]] %in% c(4, 40, 41))
  sl.indcs <- which(penindex[[1]] %in% c(1, 10, 11, 12, 13, 14))
 }else{
  sl.indg <- NULL
  sl.indcs <- NULL
 } 
 if(isordinal && (constr != "none")) stop("You must use no constraint in ordinal models")

 coef <- list()
 coef[[1]] <- matrix(0, nrow = K, ncol = P)
 ybar <- apply(dat$y, 2, function(u) weighted.mean(u, w = weights))
 if(any(ybar == 0)){
  print(ybar)
  stop("ill-conditioned data: at least one of the response categories was not observed")
 }
 ## if there are ties in the relative frequencies of the categories, break them
 ## note that at least one column of coef[[1]] must contain no ties
 if(length(unique(ybar)) != length(ybar)){
  ybar <- ybar + runif(length(ybar), 0, 1e-2)
  lyb <- sum(ybar)
  ybar <- ybar / lyb
 }
  if(!isordinal){
  if(is.numeric(constr)){
   coef[[1]][-constr,1] <- log(ybar[-constr]/(1-sum(ybar[-constr])))
   coef[[1]][constr,1] <- 0
  }else{
   coef[[1]][,1] <- log(ybar / prod(ybar)^(1/length(ybar)))
   coef[[1]][,1] <- coef[[1]][,1] - mean(coef[[1]][,1])
  }
 }else{
  coef[[1]][,1] <- model@link(mu = matrix(ybar, nrow = 1), constr = constr)
 }

 if(hasV){
  coef[[2]] <- matrix(0, nrow = K, ncol = L)
 }

 if(missing(grpindex) | is.null(grpindex)){
  grpindex <- list()
  grpindex[[1]] <- seq_len(P)
  if(hasV){
   grpindex[[2]] <- seq(from = P+1, to = P+L)
  } 
 }

 if(missing(penindex) | is.null(penindex)){
  penindex <- list()
  penindex[[1]] <- rep(1, P)
  penindex[[1]][1] <- 10                                                               
  if(hasV){
   penindex[[2]][indg] <- rep(2, length(indg))                                  
   penindex[[2]][indcs] <- rep(3, length(indcs))
  }
 } 


 if(!is.list(coef)){coef <- list(coef)}
 if(!is.list(grpindex)){grpindex <- list(grpindex)}
 if(!is.list(penindex)){penindex <- list(penindex)}

 ## set up the difference matrix 'R' for GFL-fusion-type penalties
 if(fusion == F){
  R <- NULL
 }else{
  if(is.numeric(constr)) Kf <- K-1 else Kf <- K         ## Kf = relevant size of the problem for fusion penalties
  #if(Kf < 2) stop("a fusion penalty was requested for a too small parameter group")
  if(Kf <= 2) R <- matrix(c(-1,1), nrow=1)
  if(Kf > 2){
   if(fusion == "adja"){
    R <- c(-1, 1, rep(0, Kf-2))
    for(i in 1:(Kf-2)){
     R <- rbind(R, c(rep(0,i),-1,1,rep(0,Kf-2-i)))
    }
    R <- matrix(R, ncol=Kf)
   }else if(fusion == "all"){
    R <- cbind(rep(1,Kf-1), -diag(Kf-1))
    for(i in 1:(Kf-1)){
     add <- cbind(matrix(0,nrow=Kf-1-i,ncol=i),rep(1,Kf-1-i),-diag(Kf-1-i))
     R <- rbind(R,add)
    }
   }
  }
 }#else stop("option 'fusion' in MRSP.fit must be either FALSE or ``adja`` or ``all``")
 
 ## error checking
 if(!(is.numeric(constr) | constr == "symmetric" | constr == "none"))
   stop("the identifiability constraint is misspecified")
 
 if((weights != rep(1, nrow(dat$y))) && (control@standardize = T))
   warning(paste("the usage of weights for penalized parameter groups together with standardization does not yet produce exact results. maybe try to apply the standardization manually and set control@standardize = F"))
 
 if(!is.matrix(dat$x))
   stop("x has to be a matrix")

 if(any(is.na(dat$x)))
   stop("Missing values in x not allowed!") 
 
 if(!is.numeric(dat$y))
   stop("y must be 'numeric'")
   
 if(!is.matrix(dat$y))
   stop("y must be a matrix. (wide format)")
   
 if(!model@check(dat))
   stop("y values must lie between 0 and 1, i.e. scaled response values are used. if repeated measurements are available, this is accounted for by the use of weights!")
   
 if(nrow(coef[[1]]) != K)
   stop("coef.init has wrong format")
   
 if(ncol(coef[[1]]) != P)
   stop("coef.init has wrong format")
   
 if(is.numeric(constr) && !all(coef[[1]][constr,] == 0))
   stop("coef.init does not comply with the chosen side constraint")  
   
 if(length(weights) != nobs)
   stop("weight vector has wrong length")
   
 if(any(weights <= 0))
   stop("weights must be strictly positive!")  

 if(is.matrix(offset)){
  if(nrow(offset) != nobs)
    stop("nrow(offset) not equal to the number of observations")
  if(ncol(offset) != K)
    stop("ncol(offset) not equal to the number response categoires")
 }
 if(is.vector(offset)){
  if(length(offset) != nobs)
    stop("length(offset) not equal to the number of observations")
 }
    
 if(length(grpindex[[1]]) != P)
   stop("grpindex[[1]] has wrong length")

 if(length(penindex[[1]]) != P)
   stop("penindex[[1]] has wrong length")
    
 if(hasV){
  if(!is.list(dat$V))
    stop("V has to be a list")
    
  if(length(dat$V) != K)
    stop("V has wrong length")
  
  if(nrow(dat$V[[1]]) != nobs)
    stop("the elements of V have wrong size")
    
  if(ncol(dat$V[[1]]) != L)
    stop("the elements of V have wrong size")    
    
  if(Reduce(any, lapply(dat$V, function(u){any(is.na(u))})))
    stop("Missing values in V not allowed!")
  
  if(!Reduce(all, lapply(dat$V, is.numeric)))
    stop("Values in V must be numeric")
    
  if(ncol(coef[[2]]) != L)
    stop("coef.init has wrong format")
    
  if(nrow(coef[[2]]) != K)
    stop("coef.init has wrong format")  
    
  if(length(grpindex[[2]]) != L)
    stop("grpindex[[2]] has wrong length")
    
  if(length(penindex[[2]]) != L)
    stop("penindex[[2]] has wrong length")
  
  if(length(indg) + length(indcs) != L)
    stop("indg and indcs objects have wrong length")  
 }

 
 ## creating objects that give the position, type, length etc. of the different
 ## penalization and groups. most of this is passed to Proximal.args
 
 ## global predictors
 ## if there are GFL penalties (type '14'), use type '11' instead
 penindex[[1]][which(penindex[[1]] == 14)] <- 11
 
 ## structured group lasso penalty with groups across categories
 any.grouppen.x <- any(penindex[[1]] == 1)   
 which.grouppen.x <- which(penindex[[1]] == 1)
 groups.grouppen.x <- unique(grpindex[[1]][which.grouppen.x])
 
 ## unpenalized
 any.notpen.x <- any(penindex[[1]] == 10)
 which.notpen.x <- which(penindex[[1]] == 10)
 groups.notpen.x <- unique(grpindex[[1]][which.notpen.x])
 
 ## sparse group lasso
 any.spgrppen.x <- any(penindex[[1]] == 11)
 which.spgrppen.x <- which(penindex[[1]] == 11)
 groups.spgrppen.x <- unique(grpindex[[1]][which.spgrppen.x])
 
 ## ordinary lasso
 any.lassopen.x <- any(penindex[[1]] == 12)
 which.lassopen.x <- which(penindex[[1]] == 12)
 groups.lassopen.x <- unique(grpindex[[1]][which.lassopen.x])
 
 ## ridge penalty
 any.ridgepen.x <- any(penindex[[1]] == 13)
 which.ridgepen.x <- which(penindex[[1]] == 13)
 groups.ridgepen.x <- unique(grpindex[[1]][which.ridgepen.x])
 
 ## sparse group lasso + group fused lasso
 any.spGFLpen.x <- any(penindex[[1]] == 14)
 which.spGFLpen.x <- which(penindex[[1]] == 14)
 groups.spGFLpen.x <- unique(grpindex[[1]][which.spGFLpen.x])
 
 ## ridge fusion aka discrete P-splines
 any.Psplinepen.x <- any(penindex[[1]] == 15)
 which.Psplinepen.x <- which(penindex[[1]] == 15)
 groups.Psplinepen.x <- unique(grpindex[[1]][which.Psplinepen.x])
 
 ## columnwise ridge fusion
 any.colPsplinepen.x <- any(penindex[[1]] == 16)
 which.colPsplinepen.x <- which(penindex[[1]] == 16)
 groups.colPsplinepen.x <- unique(grpindex[[1]][which.colPsplinepen.x])

 ## fused sparge group lasso (SGL + L1-fusion)
 any.spFGLpen.x <- any(penindex[[1]] == 17)
 which.spFGLpen.x <- which(penindex[[1]] == 17)
 groups.spFGLpen.x <- unique(grpindex[[1]][which.spFGLpen.x])

 ## rowwise L1 fusion
 any.rowFLpen.x <- any(penindex[[1]] == 18)
 which.rowFLpen.x <- which(penindex[[1]] == 18)
 groups.rowFLpen.x <- unique(grpindex[[1]][which.rowFLpen.x])
 
 ## set up the difference matrix 'Omega' for P-spline-type fusion penalties
 if(any.Psplinepen.x || any.rowFLpen.x){
  Omega <- list() #; length(Omega) <- length(groups.Psplinepen.x)
  tD <- list()
  for(j in sort(c(groups.Psplinepen.x, groups.rowFLpen.x))){
   j.which <- which(grpindex[[1]] == j)
   kj <- length(j.which)
   if(kj < 2) stop("a fusion penalty was requested for a too small parameter group")
   if(kj == 2){
    tD[[j]] <- c(-1,1)
    Omega[[j]] <- crossprod(matrix(c(-1,1),nrow=1,ncol=2))
   }else if(kj > 2){
    if(fusion == "adja"){
     D <- c(-1, 1, rep(0, kj-2))
     for(i in 1:(kj-2)){
      D <- rbind(D, c(rep(0,i),-1,1,rep(0,kj-2-i)))
     }
     D <- matrix(D, ncol=kj)
    }else if(fusion == "all"){
     D <- cbind(rep(1, kj-1), -diag(kj-1))
     for(i in 1:(kj-1)){
      add <- cbind(matrix(0, nrow=kj-1-i, ncol=i), rep(1, kj-1-i), -diag(kj-1-i))
      D <- rbind(D, add)
     }
    }
    tD[[j]] <- t(D)
    Omega[[j]] <- crossprod(D)
   }
  }
 }

 if(hasV){
  ## category-specific predictors with global coefficients
  ## ordinary (group-)lasso penalty
  any.globalpen.V <- any(penindex[[2]] == 2)
  which.globalpen.V <- which(penindex[[2]] == 2)
  groups.globalpen.V <- unique(grpindex[[2]][which.globalpen.V])
  
  ## unpenalized
  any.globalunpen.V <- any(penindex[[2]] == 20)
  which.globalunpen.V <- which(penindex[[2]] == 20)
  groups.globalunpen.V <- unique(grpindex[[2]][which.globalunpen.V])
  
  ## ridge
  any.globalridgepen.V <- any(penindex[[2]] == 21)
  which.globalridgepen.V <- which(penindex[[2]] == 21)
  groups.globalridgepen.V <- unique(grpindex[[2]][which.globalridgepen.V])
  
  ## category-specific predictors with category-specific coefficients 
  ## structured group lasso like in "1"
  any.catpen.V <- any(penindex[[2]] == 3)
  which.catpen.V <- which(penindex[[2]] == 3)
  groups.catpen.V <- unique(grpindex[[2]][which.catpen.V])
  
  ## unpenalized
  any.catunpen.V <- any(penindex[[2]] == 30)
  which.catunpen.V <- which(penindex[[2]] == 30)
  groups.catunpen.V <- unique(grpindex[[2]][which.catunpen.V])
  
  ## sparse group lasso
  any.catspgrppen.V <- any(penindex[[2]] == 31)
  which.catspgrppen.V <- which(penindex[[2]] == 31)
  groups.catspgrppen.V <- unique(grpindex[[2]][which.catspgrppen.V])
  
  ## ordinary lasso
  any.catlassopen.V <- any(penindex[[2]] == 32)
  which.catlassopen.V <- which(penindex[[2]] == 32)
  groups.catlassopen.V <- unique(grpindex[[2]][which.catlassopen.V])
  
  ## ridge
  any.catridgepen.V <- any(penindex[[2]] == 33)
  which.catridgepen.V <- which(penindex[[2]] == 33)
  groups.catridgepen.V <- unique(grpindex[[2]][which.catridgepen.V])
 }

 if(isordinal){
  ## ordinal model with global effects
  ## ordinary (group-)lasso penalty
  any.globalpen.o <- any(penindex[[1]] == 4)
  which.globalpen.o <- which(penindex[[1]] == 4)
  groups.globalpen.o <- unique(grpindex[[1]][which.globalpen.o])
  
  ## unpenalized
  any.globalunpen.o <- any(penindex[[1]] == 40)
  which.globalunpen.o <- which(penindex[[1]] == 40)
  groups.globalunpen.o <- unique(grpindex[[1]][which.globalunpen.o])
  
  ## ridge
  any.globalridgepen.o <- any(penindex[[1]] == 41)
  which.globalridgepen.o <- which(penindex[[1]] == 41)
  groups.globalridgepen.o <- unique(grpindex[[1]][which.globalridgepen.o])
 } 
   

 ## degrees of freedom, groups and so on...
 if(isordinal){
  pen.which.x <- sort(c(which.grouppen.x, which.spgrppen.x, which.lassopen.x, which.ridgepen.x,
                        which.spGFLpen.x, which.spFGLpen.x, which.rowFLpen.x,
                        which.globalpen.o, which.globalridgepen.o))
  notpen.which.x <- sort(c(which.notpen.x, which.globalunpen.o))
  Pspline.which.x <- sort(c(which.Psplinepen.x, which.colPsplinepen.x))
 }else{
  pen.which.x <- sort(c(which.grouppen.x, which.spgrppen.x, which.lassopen.x, which.ridgepen.x,
                        which.spGFLpen.x, which.spFGLpen.x, which.rowFLpen.x))
  notpen.which.x <- sort(c(which.notpen.x))
  Pspline.which.x <- sort(c(which.Psplinepen.x, which.colPsplinepen.x))
 }
 nrpen.x <- length(pen.which.x)
 nrnotpen.x <- length(notpen.which.x)
 nrPspline.x <- length(Pspline.which.x)

 grpindex.pen.x <- grpindex[[1]][pen.which.x]
 nrgrp.pen.x <- length(unique(grpindex.pen.x))
 grpindex.notpen.x <- grpindex[[1]][notpen.which.x]
 nrgrp.notpen.x <- length(unique(grpindex.notpen.x))
 grpindex.Pspline.x <- grpindex[[1]][Pspline.which.x]
 nrgrp.Pspline.x <- length(unique(grpindex.Pspline.x))

 dict.pen.x <- sort(unique(grpindex.pen.x))
 pen.tab.x <- table(grpindex.pen.x)[as.character(dict.pen.x)]
  
 if(hasV){
  pen.which.V <- sort(c(which.globalpen.V, which.globalridgepen.V, which.catpen.V, which.catspgrppen.V, which.catlassopen.V, which.catridgepen.V))
  nrpen.V <- length(pen.which.V)
  notpen.which.V <- sort(c(which.globalunpen.V, which.catunpen.V))
  nrnotpen.V <- length(notpen.which.V)

  grpindex.pen.V <- grpindex[[2]][pen.which.V]
  nrgrp.pen.V <- length(unique(grpindex.pen.V))
  grpindex.notpen.V <- grpindex[[2]][notpen.which.V]
  nrgrp.notpen.V <- length(unique(grpindex.notpen.V))

  dict.pen.V <- sort(unique(grpindex.pen.V))
  pen.tab.V <- table(grpindex.pen.V)[as.character(dict.pen.V)]
 }  

 ## remove all penalized predictors from the model
 dat.unpen <- dat
 dat.unpen$x <- dat$x[,penindex[[1]] %in% c(10,15,20,30,40), drop=F]
 penindex.unpen <- penindex
 penindex.unpen[[1]] <- penindex[[1]][penindex[[1]] %in% c(10,15,20,30,40)]
 if(isordinal){
  if(!any(penindex[[1]] == 40)){
   sl.indg.old <- sl.indg
   sl.indg <- NULL
  }else{
   sl.indg.old <- sl.indg
   sl.indg <- which(penindex[[1]] == 40)
  }
  if(!any(penindex[[1]] == 10)){
   sl.indcs.old <- sl.indcs
   sl.indcs <- NULL
  }else{
   sl.indcs.old <- sl.indcs
   sl.indcs <- which(penindex[[1]] == 10)
  }   
 }  
 grpindex.unpen <- grpindex
 grpindex.unpen[[1]] <- grpindex[[1]][penindex[[1]] %in% c(10,15,20,30,40)]
 coef.unpen <- coef
 coef.unpen[[1]] <- coef[[1]][,penindex[[1]] %in% c(10,15,20,30,40), drop=F]
 
 if(hasV){
  if(length(notpen.which.V) > 0){
   dat.unpen$V <- lapply(dat$V, function(u) u[, notpen.which.V, drop=F])
   penindex.unpen[[2]] <- penindex[[2]][penindex[[2]] %in% c(10,20,30)]
   grpindex.unpen[[2]] <- grpindex[[2]][penindex[[2]] %in% c(10,20,30)]
   coef.unpen[[2]] <- coef[[2]][penindex[[2]] %in% c(10,20,30)] 
  }else if(length(notpen.which.V) == 0){
   dat.unpen$V <- NULL
   penindex.unpen[[2]] <- NULL
   grpindex.unpen[[2]] <- NULL
   coef.unpen[[2]] <- NULL
   hasV <- F
  } 
 }
 
 ## signal fistaProximal and penalty that all the penalized groups are gone
 any.grouppen.x = F
 any.spgrppen.x = F
 any.lassopen.x = F
 any.ridgepen.x = F
 any.spGFLpen.x = F
 any.spFGLpen.x = F
 #any.rowFLpen.x = F
 if(hasV){
  any.globalpen.V = F
  any.globalridgepen.V = F
  any.catpen.V = F
  any.catspgrppen.V = F
  any.catlassopen.V = F
  any.catridgepen.V = F
 }
 if(isordinal){
  any.globalpen.o = F
  any.globalridgepen.o = F
 } 


 ## extract control information
 max.iter <- control@max.iter
 rel.tol <- control@rel.tol
 standardize <- control@standardize
 ridgestabil <- control@ridgestabil
 ridgestabilrf <- control@ridgestabilrf
 lambdastabil <- control@lambdastabil
 SPG.iter.max <- control@SPG.iter.max
 SPG.accel <- control@SPG.accel
 SPG.eps <- control@SPG.eps
  
 ## extract the model information
 link <- model@link
 invlink <- model@invlink
 loglik <- model@loglik
 gradient <- model@gradient
 fisher <- model@fisher
 name <- model@name
 
 ## assign the proxy class to the coef object
 class(coef.unpen) <- "MRSP.coef"

 eta <- updateEta(dat=dat.unpen, coef=coef.unpen, offset=offset)
 mu <- invlink(eta)
 
 if(missing(penweights) | is.null(penweights)) stop("penweights missing in getlambdamax")
 
 ## list with arguments for the calls to fistaProximal and the other fistaGenerics
 if(hasV){Proximal.args <- list(
                      any.grouppen.x = any.grouppen.x,
                      any.notpen.x = any.notpen.x,   
                      any.spgrppen.x = any.spgrppen.x,
                      any.lassopen.x = any.lassopen.x,
                      any.ridgepen.x = any.ridgepen.x,
                      any.spGFLpen.x = any.spGFLpen.x,
                      any.Psplinepen.x = any.Psplinepen.x,
                      any.globalpen.V = any.globalpen.V,
                      any.globalunpen.V = any.globalunpen.V,
                      any.globalridgepen.V = any.globalridgepen.V,
                      any.catpen.V = any.catpen.V,
                      any.catunpen.V = any.catunpen.V,
                      any.catspgrppen.V = any.catspgrppen.V,
                      any.catlassopen.V = any.catlassopen.V,
                      any.catridgepen.V = any.catridgepen.V,
                      which.grouppen.x = which.grouppen.x,
                      which.notpen.x = which.notpen.x,   
                      which.spgrppen.x = which.spgrppen.x,
                      which.lassopen.x = which.lassopen.x,
                      which.ridgepen.x = which.ridgepen.x,
                      which.spGFLpen.x = which.spGFLpen.x,
                      which.Psplinepen.x = which.Psplinepen.x,
                      which.globalpen.V = which.globalpen.V,
                      which.globalunpen.V = which.globalunpen.V,
                      which.globalridgepen.V = which.globalridgepen.V,
                      which.catpen.V = which.catpen.V,
                      which.catunpen.V = which.catunpen.V,
                      which.catspgrppen.V = which.catspgrppen.V,
                      which.catlassopen.V = which.catlassopen.V,
                      which.catridgepen.V = which.catridgepen.V,
                      groups.grouppen.x = groups.grouppen.x,
                      groups.notpen.x = groups.notpen.x,
                      groups.spgrppen.x = groups.spgrppen.x,
                      groups.lassopen.x = groups.lassopen.x,
                      groups.ridgepen.x = groups.ridgepen.x,
                      groups.spGFLpen.x = groups.spGFLpen.x,
                      groups.Psplinepen.x = groups.Psplinepen.x,
                      groups.globalpen.V = groups.globalpen.V,
                      groups.globalunpen.V = groups.globalunpen.V,
                      groups.globalridgepen.V = groups.globalridgepen.V,
                      groups.catpen.V = groups.catpen.V,
                      groups.catunpen.V = groups.catunpen.V,
                      groups.catspgrppen.V = groups.catspgrppen.V,
                      groups.catlassopen.V = groups.catlassopen.V,
                      groups.catridgepen.V = groups.catridgepen.V,
                      any.colPsplinepen.x = any.colPsplinepen.x,
                      which.colPsplinepen.x = which.colPsplinepen.x,
                      groups.colPsplinepen.x = groups.colPsplinepen.x,
                      any.spFGLpen.x = any.spFGLpen.x,
                      which.spFGLpen.x = which.spFGLpen.x,
                      groups.spFGLpen.x = groups.spFGLpen.x,
                      any.rowFLpen.x = any.rowFLpen.x,
                      which.rowFLpen.x = which.rowFLpen.x,
                      groups.rowFLpen.x = groups.rowFLpen.x,
                      dict.pen.x = dict.pen.x,
                      dict.pen.V = dict.pen.V,
                      gamma = gamma,
                      psi = psi,
                      hasV = hasV,
                      isordinal = isordinal,
                      ridgestabil = ridgestabil,
                      ridgestabilrf = ridgestabilrf,
                      lambdastabil = lambdastabil,
                      extrastabil = control@extrastabil,
                      sancount = 0,
                      n = nobs,
                      constraint = constr,
                      indg = indg,
                      indcs = indcs,
                      sl.indg = sl.indg,
                      sl.indcs = sl.indcs,
                      R = R,
                      SPG.iter.max = SPG.iter.max,
                      SPG.accel = SPG.accel,
                      SPG.eps = SPG.eps,
                      modelname = model@name,
                      nonneg = nonneg)
 }else{Proximal.args <- list(
                      any.grouppen.x = any.grouppen.x,
                      any.notpen.x = any.notpen.x,   
                      any.spgrppen.x = any.spgrppen.x,
                      any.lassopen.x = any.lassopen.x,
                      any.ridgepen.x = any.ridgepen.x,
                      any.spGFLpen.x = any.spGFLpen.x,
                      any.Psplinepen.x = any.Psplinepen.x,
                      which.grouppen.x = which.grouppen.x,
                      which.notpen.x = which.notpen.x,   
                      which.spgrppen.x = which.spgrppen.x,
                      which.lassopen.x = which.lassopen.x,
                      which.ridgepen.x = which.ridgepen.x,
                      which.spGFLpen.x = which.spGFLpen.x,
                      which.Psplinepen.x = which.Psplinepen.x,
                      groups.grouppen.x = groups.grouppen.x,
                      groups.notpen.x = groups.notpen.x,
                      groups.spgrppen.x = groups.spgrppen.x,
                      groups.lassopen.x = groups.lassopen.x,
                      groups.ridgepen.x = groups.ridgepen.x,
                      groups.spGFLpen.x = groups.spGFLpen.x,
                      groups.Psplinepen.x = groups.Psplinepen.x,
                      any.colPsplinepen.x = any.colPsplinepen.x,
                      which.colPsplinepen.x = which.colPsplinepen.x,
                      groups.colPsplinepen.x = groups.colPsplinepen.x,
                      any.spFGLpen.x = any.spFGLpen.x,
                      which.spFGLpen.x = which.spFGLpen.x,
                      groups.spFGLpen.x = groups.spFGLpen.x,
                      any.rowFLpen.x = any.rowFLpen.x,
                      which.rowFLpen.x = which.rowFLpen.x,
                      groups.rowFLpen.x = groups.rowFLpen.x,
                      dict.pen.x = dict.pen.x,
                      gamma = gamma,
                      psi = psi,
                      hasV = hasV,
                      isordinal = isordinal,
                      ridgestabil = ridgestabil,
                      ridgestabilrf = ridgestabilrf,
                      lambdastabil = lambdastabil,
                      extrastabil = control@extrastabil,
                      sancount = 0,
                      n = nobs,
                      constraint = constr,
                      indg = indg,
                      indcs = indcs,
                      sl.indg = sl.indg,
                      sl.indcs = sl.indcs,
                      R = R,
                      SPG.iter.max = SPG.iter.max,
                      SPG.accel = SPG.accel,
                      SPG.eps = SPG.eps,
                      modelname = model@name,
                      nonneg = nonneg)
 }
 if(isordinal){
  Proximal.args$any.globalpen.o <- any.globalpen.o
  Proximal.args$any.globalunpen.o <- any.globalunpen.o
  Proximal.args$any.globalridgepen.o <- any.globalridgepen.o
  Proximal.args$which.globalpen.o <- which.globalpen.o
  Proximal.args$which.globalunpen.o <- which.globalunpen.o
  Proximal.args$which.globalridgepen.o <- which.globalridgepen.o
  Proximal.args$groups.globalpen.o <- groups.globalpen.o
  Proximal.args$groups.globalunpen.o <- groups.globalunpen.o
  Proximal.args$groups.globalridgepen.o <- groups.globalridgepen.o
 }
 if(any.Psplinepen.x || any.rowFLpen.x){
  Proximal.args$Omega <- Omega
  Proximal.args$tD <- tD
 }

 tuning <- list(lambda = 0,
                gamma  = gamma,
                psi    = psi,
                lambdaR = 0,
                lambdastabil = 1,
                lambdawhatever = 0,
                lambdaF = lambdaF)


 ## create the fista.control argument
 if(missing(fista.control) | is.null(fista.control)){
  fista.control <- fista.control(rel.tol = rel.tol, max.iter = max.iter)
 }else{fista.control@rel.tol <- rel.tol; fista.control@max.iter <- max.iter} 

 ## now the call to fista where the core fitting takes place.
 updated <- fista(coef.init=coef.unpen, dat=dat.unpen, offset=offset, eta.init=eta,
                  mu.init=mu, weights=weights, model=model, control=fista.control,
                  Proximal.args=Proximal.args, #Proximal.control=Proximal.control,
                  grpindex=grpindex.unpen, penindex=penindex.unpen, tuning=tuning,
                  penweights.init=penweights)
                 
 coef0 <- updated$coef
 penweights0 <- updated$penweights
 eta0 <- updated$eta
 mu0 <- updated$mu

 ## a little helper function
 l2norminv <- cmpfun(function(u){sqrt(sum(c(u)^2) / length(c(u)))})
 l2norminv.x <- cmpfun(function(u){
 if(is.matrix(u)){sqrt(sum(c(u)^2) / (length(c(u))-ncol(u)))
  }else{stop("error in l2norminv: non-matrix argument")}}) 
 
 if(constr == "none"){
  l2norminvx <- l2norminv
 }else{
  l2norminvx <- l2norminv.x
 } 

 if(isordinal){
  Proximal.args$sl.indg <- sl.indg.old
  Proximal.args$sl.indcs <- sl.indcs.old
 }

 ## the gradient based on the model where all the penalized parameters are 0
 grad0 <- gradient(dat = dat, mu = mu0, weights = weights, Proximal.args = Proximal.args)

 ## now the norms of the gradient of the penalized groups, based on the null model.
 ## gradnorms is for group lasso, gradnormcoefs for ordinary lasso.
 gradnorms.x <- numeric(nrgrp.pen.x)
 gradnormcoefs.x <- matrix(nrow = ncol(dat$y), ncol = nrgrp.pen.x)
 for(j in unique(grpindex.pen.x)){
  j.which <- which(grpindex[[1]] == j) 
  gradnorms.x[which(unique(grpindex.pen.x) == j)] <- l2norminvx(grad0[[1]][, j.which, drop=F]) / penweights[[1]][[1]][j]
  if(isordinal && !is.null(sl.indg.old)){
   if(all(j.which %in% sl.indg.old)){
    gradnorms.x[which(unique(grpindex.pen.x) == j)] <- l2norminvx(grad0[[1]][1,j.which, drop=F]) / penweights[[1]][[1]][j]
   }
  }    
  for(i in seq(ncol(dat$y))){
   gradnormcoefs.x[i, which(unique(grpindex.pen.x) == j)] <- l2norminv(grad0[[1]][i, j.which, drop=F]) / penweights[[2]][[1]][i,j]
  }   
 }

 if(hasV){
  lp <- length(unique(grpindex[[1]]))
  gradnorms.V <- numeric(nrgrp.pen.V)
  gradnormcoefs.V <- matrix(nrow = ncol(dat$y), ncol = nrgrp.pen.V) 
  for(j in unique(grpindex.pen.V)){
   j.which <- which(grpindex[[2]] == j)
   if(all(j.which %in% indg)){
    gradnorms.V[which(unique(grpindex.pen.V) == j)] <- l2norminv(grad0[[2]][1,j.which]) / penweights[[1]][[2]][j-lp]
   }else if(all(j.which %in% indcs)){
    gradnorms.V[which(unique(grpindex.pen.V) == j)] <- l2norminv(grad0[[2]][,j.which]) / penweights[[1]][[2]][j-lp]
   }else{
    stop("something went wrong with j.which, indg and indcs in getlambdamax")
   }
   for(i in seq(ncol(dat$y))){
    gradnormcoefs.V[i, which(unique(grpindex.pen.V) == j)] <- l2norminv(grad0[[2]][i, j.which, drop=F]) / penweights[[2]][[2]][i,j-lp]  
   }
  } 
 }

 lambdamax <- max(gradnorms.x)

 if(any(penindex[[1]] %in% c(11, 12, 14))) lambdamax <- max(lambdamax, gradnormcoefs.x)
 if(hasV){
  lambdamax <- max(lambdamax, gradnorms.V)
  if(any(penindex[[2]] %in% c(31, 32))) lambdamax <- max(lambdamax, gradnormcoefs.V)
 }

 return(lambdamax)
}


## a function for simulating multinomial data. N = sample size, K = number of
## response categories, coef = coefficients for which to simulate.
#multinom.simdata <- function(coef, nobs, rng, predictor.types, model = multinomlogit(),
#                             weights = rep(1, N), offset = rep(0, N), ...){
# K <- nrow(coef[[1]])
# hasV <- length(coef) > 1
# link <- model@link
# invlink <- model@invlink
# P <- ncol(coef[[1]])
# if(hasV) L <- ncol(coef[[2]])
#
# X <- matrix(nrow=nobs, ncol=1)
# for(j in seq(length(predictor.types))){
#  if(predictor.types[j] == 0){
#   Xj <- eval(as.call(c(as.symbol(rng), c(nobs, list(...)))))
# }
# X[,1] <- 1
# }
# if(length(coef) == 2){
#  L <- length(coef[[2]])
#  V <- list(); length(V) <- K
#  for(k in seq(K)){
#   V[[k]] <- matrix(nrow=nobs,ncol=L)
#   for(j in seq(L)){
#    V[[k]][,j] <- rnorm(nobs, sd=2)
#   }
#  # V[[k]] <- sweep(V[[k]], 2, V[[k]][1, ,drop=T])
#  # V[[k]] <- V[[k]][-1, ,drop=F]
#  }
# }else{V <- NULL}
#
# dat <- list(y = matrix(nrow=N,ncol=K),
#             x = X,
#             V = V)
#
# class(coef) <- "MRSP.coef"
#
# eta <- updateEta(dat=dat, coef=coef, offset=offset, weights=weights)
# mu <- invlink(eta)
# muk <- mu #cbind(mu, 1-rowSums(mu))
# muk[muk < 0] <- 0
#
# y <- matrix(nrow=N, ncol=K)
# for(i in seq(N)){
#  y[i,] <- rmultinom(1, size=weights[i], prob = muk[i,])
# }
# y <- sweep(y, 1, weights, FUN = "/")
# #y <- y[,-K]
# dat$y <- y
#
# return(dat)
#}

## function for simulating coefficients. K = number of categories,
## P = number of global predictors, L = number of category-specific predictors
multinom.simcoef <- function(K, P, L=NULL){
 coef <- list()
 coef[[1]] <- matrix(nrow=K, ncol= P)
 coef[[1]][,1] <- sample(-5:5, size=K, replace=T)
 for(j in 2:P){
  coef[[1]][,j] <- sample(-3:3, size=K, replace=T)
 }
 if(!is.null(L)){
  for(j in 1:L){
   coef[[2]][,j] <- sample(-3:3, size=L, replace=T)
  } 
 }
 coef
}

## helper function for extracting stuff from MRSP.list objects.
## object is an MRSP.list object, what is a character giving the slot to extract
setMethod("extract",
          signature(object = "MRSP.list"),
          function(object, what, ...)
{
 what <- as.character(what)
 if(length(what) != 1)
   stop(" 'extract' cant extract more than one slot at a time")
 if(is.call(object[[length(object)]])) object <- object[-length(object)]        ## this line shouldnt be needed anymore, but I'll leave it for now. In old versions of MRSP,
                                                                                ## the overall call to MRSP.fit was sometimes stored in an additional element of the output list of class MRSP.fit
 out <- lapply(object, function(u) slot(u, what))
 if((!is.matrix(out[[1]]) | !is.list(out[[1]]) | !is.array(out[[1]])) & (length(out[[1]]) == 1))
 out <- Reduce("c", out)   
 return(out)
})

## using the following is ineffective, but its existence makes extract a bit more foolproof:
setMethod("extract",
          signature(object = "MRSP"),
          function(object, what, ...)
{
 what <- as.character(what)
 if(length(what) != 1)
   stop(" 'extract' cannot extract more than one slot at a time")
 out <- slot(object, what)
 return(out)
})


## helper function that creates responses given a set of predictors and a coef 
## object of class MRSP.coef
setMethod("createResponse",
          signature(coef = "MRSP.coef"),
          function(coef, dat, model = multinomlogit(), weights, offset, ...)
{
 isordinal <- model@name %in% c("Cumulative Logit Model", "Sequential Logit Model")
 sane.model <- 0
 while(sane.model < 1){
  if(!is.list(dat)) dat <- list(dat)
  if(!is.list(coef)) coef <- list(coef)
  hasV <- !is.null(dat$V)
  names(dat)[which(sapply(dat, is.matrix))] <- "x"
  if(hasV) names(dat)[which(sapply(dat, is.list))] <- "V"
  nobs <- nrow(dat$x)
  K <- nrow(coef[[1]])
  dat$y <- matrix(nrow=nobs, ncol=K)
  if(missing(weights)) weights <- rep(1, nobs)
  if(missing(offset)) offset <- rep(0, nobs)
  if(hasV & length(coef) < 2) stop("dat and coef don't match")
  if(ncol(coef[[1]]) != ncol(dat$x)) stop("dat and coef don't match")

  invlink <- model@invlink

  eta <- updateEta(dat=dat, coef=coef, offset=offset, weights=weights)
  mu <- invlink(eta)

  if(isordinal){
   mu <- cbind(mu, 1 - rowSums(mu))
   K <- K + 1
  }
  mu[mu < 0] <- 0

  y <- matrix(nrow=nobs, ncol=K)
  for(i in seq(nobs)){
   y[i,] <- rmultinom(1, size = weights[i], prob = mu[i,])
  }
  y <- sweep(y, 1, weights, FUN = "/")
  if(all(colSums(y) >= 2)) sane.model <- 1
  if(isordinal) y <- y[,-K]
 }

 dat$y <- y

 return(dat)
})

## helper function that computes the MSE from an MRSP.sim object:
setMethod("MSE",
          signature(object = "MRSP.sim"),
          function(object, type = "coef", scaled = TRUE, ...)
{
 ncall <- length(object[[1]])
 truecoef <- object[[length(object)]] ## the last entry of an "MRSP.sim" object is the true coef vector!
 object <- object[-length(object)]
 truemus <- lapply(object, function(u){u[[1]]$truemu}) ## the truemu is the same for all calls, therefore
                                                       ## it is always only stored for the first call. 
 
 MSEs <- list(); length(MSEs) <- ncall
 for(i in seq(ncall)){
  if(type == "coef"){
   coefs <- lapply(object, function(u){(u[[i]]$bestfit)@coef})
   MSEs[[i]] <- sapply(coefs, function(u){Reduce('+', lapply(Map('-', u, truecoef), function(v){sum(v^2)})) / ifelse(scaled, sum(sapply(truecoef, function(w){length(c(w))})), 1)})
  }
  if(type == "mu"){
   mus <- lapply(object, function(u){(u[[i]]$bestfit)@mu})
   MSEs[[i]] <- mapply(function(u, v){mean((u - v)^2)}, mus, truemus)
  }
  if(type == "pred.error"){
   prederrs <- sapply(object, function(u){u[[i]]$pred.err})
   MSEs[[i]] <- prederrs #/ ifelse(scaled, length(prederrs), 1)
  }
 }
  
 return(MSEs)
})

## helper function that computes the FPR and FNR from an MRSP.sim object
## first a helper function
FPR.FNR.MRSP.workhorse <- function(object, type)
{
 if(class(object) != "MRSP.sim") stop("FPR.FNR applied to an object of wrong class")
 ncall <- length(object[[1]])
 if(ncall > 1 & length(type) == 1) type = rep(type, ncall)
 truecoef <- object[[length(object)]] ## the last entry of an "MRSP.sim" object is the true coef vector!
 object <- object[-length(object)]

 getzeroes <- function(coef, type){switch(type,
                                  "groups" = lapply(coef, function(u){if(is.matrix(u)){apply(u, 2, function(w) all(w == 0))}else{u == 0}}),
                                  "coefs" = lapply(coef, function(u){u == 0}))
                                 }
                                   
 getnonzeroes <- function(coef, type){switch(type,
                                  "groups" = lapply(coef, function(u){if(is.matrix(u)){apply(u, 2, function(w) any(w != 0))}else{u != 0}}),
                                  "coefs" = lapply(coef, function(u){u != 0}))
                                 }
                                 
 fprs <- vector(length = ncall)
 fnrs <- vector(length = ncall)
 
 for(i in seq(ncall)){
  truezeroes <- getzeroes(truecoef, type[i])
  ntruezeroes <- Reduce('+', lapply(truezeroes, sum))
  
  truenonzeroes <- getnonzeroes(truecoef, type[i])
  ntruenonzeroes <- Reduce('+', lapply(truenonzeroes, sum))
                        
  coefs <- lapply(object, function(u){(u[[i]]$bestfit)@coef})
  fpr <- sapply(lapply(coefs, function(u){getnonzeroes(u, type[i])}), function(v){Reduce('+', lapply(mapply(function(x, y){(x == T) & (y == T)}, v, truezeroes), sum))})
  fpr <- sapply(fpr, function(u) ifelse(ntruezeroes == 0, 0, u / ntruezeroes))
  fprs[i] <- mean(fpr) 
  fnr <- sapply(lapply(coefs, function(u){getzeroes(u, type[i])}), function(v){Reduce('+', lapply(mapply(function(x, y){(x == T) & (y == T)}, v, truenonzeroes), sum))})
  fnr <- sapply(fnr, function(u) ifelse(ntruenonzeroes == 0, 0, u / ntruenonzeroes))
  fnrs[i] <- mean(fnr)
 }
 
 return(rbind(fprs, fnrs))
} 

## now the true method
setMethod("FPR.FNR",
          signature(object = "MRSP.sim"),
          function(object, type = "mixed", ...)
{
 if(type != "mixed"){
  out <- FPR.FNR.MRSP.workhorse(object = object, type = type)
 }else{
  out <- (FPR.FNR.MRSP.workhorse(object = object, type = "groups") + FPR.FNR.MRSP.workhorse(object = object, type = "coefs"))/2
 }
 return(out)
})   

## generic which.min.cv to extract the min of a cv object with a certain tolerance
## arguments:
## x:  a vector with values from which to get the min position.
## reltol, abstol: relative and absolute tolerance.
## left: should the search for tolerable x be started from "the left" or the right?
setMethod("which.min.cv",
          signature(x = "numeric"),
          function(x, reltol = 0.01, abstol = NULL, left = TRUE, ...)
{
 xmin <- min(x)
 besti <- which.min(x)
 xtest <- x[besti]
 bestabs <- besti
 bestrel <- besti
 
 if(left){
  if(!is.null(abstol)){
   for(i in seq(1, besti)){
    xtest <- x[i]
    if(xtest - xmin < abstol){bestabs <- i; break}
   }
  }else{bestabs <- 1}
  for(i in seq(1, besti)){
   xtest <- x[i]
   if((xtest - xmin)/xmin < reltol){bestrel <- i; break}
  }
  bestm <- max(bestabs, bestrel)
 }else{
  if(!is.null(abstol)){
   for(i in seq(length(x), besti)){
    xtest <- x[i]
    if(xtest - xmin < abstol){bestabs <- i; break}
   }
  }else{bestabs <- length(x)}
  for(i in seq(length(x), besti)){
   xtest <- x[i]
   if((xtest - xmin)/xmin < reltol){bestrel <- i; break}
  }
  bestm <- min(bestabs, bestrel)
 }   

 if(!(bestm %in% seq_along(x))) stop("something went wrong in which.min.cv")  
 return(bestm)
})  


## which.max.cv
setMethod("which.max.cv",
          signature(x = "numeric"),
          function(x, reltol = 0.01, abstol = 4, left = TRUE, ...)
{
 which.min.cv(-x, reltol = reltol, abstol = abstol, left = left, ...)
})       