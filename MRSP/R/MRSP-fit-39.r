################################################################################
#####    An R Package for fitting (M)ultinomial (L)ogit models with a      #####
#####    (S)tructured (P)enalization for global and category-specific      #####
#####    predictors.                                                       #####
################################################################################
#####    Author: Wolfgang Pößnecker                                        #####
#####    Last modified: 14.10.2014, 16:21                                  #####
################################################################################

## the main fitting function:
## -------------------------------------------------------------------------- ##
## Arguments: 
   ## dat: a list that contains the data. must have entries "y", "x" and "V".
   ##      "y" is a matrix with entries between 0 and 1 giving the response. in 
   ##      the case of repeated observations, the entries in y are the relative 
   ##      frequency of the response categories for a particular combination of
   ##      covariates. the amount of repeated measurements is then captured by 
   ##      weights. size nobs x K for K response categories.
   ##      "x" is a matrix with predictors that dont depend on the response 
   ##      category but whose coefficients are category-specific. size nobs x P.
   ##      "V" is a list (!) of length K whose entires are nobs x L matrices 
   ##      giving the category-specific predictors. 
   ## coef.init: a list of initial coefficients. the first entry is a K x P 
   ##      matrix corresponding to x, the second entry is a K x L matrix
   ##      belonging to V.
   ## coef.stand.init: same as coef.init, but corresponding to the standardized
   ##      predictors.
   ## offset: a vector of length nobs with offset values.
   ## weights: a vector of observation weights.
   ## grpindex: a list of two vectors that indicate which columns of the design
   ##      matrix form a group that has to be penalized jointly, e.g. the 
   ##      different dummies of a categorical predictor. the first element is 
   ##      the grouping vector for x, the second one for V. those columns with 
   ##      the same number belong to one group. the numbers must begin with 1!
   ##      MRSP always uses a group lasso for such parameter groups. if you dont
   ##      want this behaviour, simly dont specify grpindex.
   ## penindex: a list of two vectors, same size as grpindex, whose numbers 
   ##      indicate how each predictors shall be penalized. elements that share 
   ##      the same number must obviously also share the same number in penindex
   ##      possible values:
   ##      1:  global predictor ("x") whose coefficients shall be penalized with  
   ##          a group lasso penalty with grouping "across" categories.
   ##      10: global predictor, unpenalized.
   ##      11: global predictor, sparse group lasso.
   ##      12: global predictor, ordinary lasso.
   ##      13: global predictor, ridge penalty. does not support penweights.
   ##      14: global predictor, sparse group lasso + group fused lasso.
   ##      15: global predictor, "rowwise" ridge fusion on adjacent parameters.
   ##          (i.e. a discrete P-spline penalty.) only works for parameter
   ##          groups with more than one element.
   ##      16: global predictor, "columnwise" ridge fusion on adjacent para-
   ##          meters. only works if K > 1.
   ##      17: global predictor, sparge group lasso + L1 fused lasso.
   ##      18: global predictor, "rowwise" L1 fused lasso.
   ##      2:  category-specific predictor whose coefficients are penalized with
   ##          the ordinary (group-)lasso. (depending on grpindex.)
   ##      20: category-specific, unpenalized.
   ##      21: category-specific, ridge penalty.
   ##      3:  category-speficic predictor with category-specific coefficients
   ##          that are penalized by a group lasso like in "1". this variant is
   ##          not implemented yet!
   ##      30: category-specific with category-specific coefs, unpenalized.
   ##      31: category-specific with category-specific coefficients and sparse
   ##          group lasso penalty.
   ##      32: cat-cat-specific, with ordinary lasso.
   ##      33: cat-cat, with ridge. does not support penweights.
   ##      34-37: to come in the future, like the previous four versions, but
   ##             includes automatic distinction between category-specific and
   ##             global parameterizations for category-specific preditors.
   ##             in other words: these 4 will be able to shrink the "3"-series
   ##             to the "2"-series.
   ##      4:  ordinal, global effect, penalized. (cf. the "2"-series).
   ##      40: ordinal, global effect, unpenalized.
   ##      41: ordinal, global effect, ridge penalty. 
   ## tuning: a list of tuning parameters for the model. the first slot contains
   ##      the tuning parameters for the global predictors in "x", the second 
   ##      slot that for the category-specific predictors with global parameters
   ##      and the third slot those for the "3"-series from above. if a sparse 
   ##      group lasso is used in the "1"- or the "3"-series, the corresponding 
   ##      slot in tuning must contain 2 parameters instead of just one.
   ## model: an object of class MRSP.model.
   ## constr: either "symmetric" or "none", which mean that a symmetric side 
   ##         constraint or no side constraint are used. (if constr = "none",
   ##         mean-centering is used for coefficients of unpenalized predictors)
   ##         alternatively, constr can be numeric; then the corresponding 
   ##         category of y is the reference category.
   ## control: an object of class MRSP.control.
   ## fista.control: a list with control information for fista.
   ## Proximal.control: a list with control information for fistaProximal.
   ## Proximal.args: a list with information for the call to fistaProximal during
   ##      the fista algorithm. if missing, the required information will be 
   ##      computed within MRSP.fit.
   ## penweights: a list with 2 elements, one for the group lasso type penalties
   ##      and one for the ordinary lasso type penalties. each of these elements
   ##      gives the weight that is assigned to the various predictors. this is 
   ##      used for adaptive lasso estimation.
   ## adaptive: should an adaptive penalty be used? set to F to not use adaptive
   ##      penalization. set to T to use the same penalty for both initial and
   ##      the adaptively weighted final fit. set to "ML" to use the ML-fit as 
   ##      as the initial estimator. 
   ## refit: should a refit be performed?
   ## fusion: if fusion="adja", fusion-type penalties apply to all differences
   ##      between adjacent coefficients. if fusion="all", it applies to all
   ##      pairwise differences. fusion = FALSE means no fusion is considered.
   ## nonneg: if TRUE, a nonnegativity constraint is incorporated on all coef-
   ##      ficients that are subject to a grouped or ordinary lasso penalty.
## -------------------------------------------------------------------------- ##

## the main fitting method which fits the model for one conrete lambda:
setMethod("MRSP.fit",
          signature(lambda = "numeric"),
function(dat, coef.init, coef.stand.init, coef.pretres.init, offset = rep(0, nrow(dat$y)),
                     weights = rep(1, nrow(dat$y)), grpindex = NULL, penindex = NULL,
                     lambda, lambdaR = lambda, lambdaF = lambda, gamma = 1, psi = 1,
                     indg = NULL, indcs = NULL,
                     model = multinomlogit(), constr = NULL, control = MRSP.control(),
                     fista.control = NULL, Proximal.control = NULL,
                     Proximal.args = NULL, penweights = NULL, mlfit = NULL,
                     adaptive = FALSE, threshold = FALSE, refit = FALSE, fusion = FALSE, nonneg = FALSE,
                     ...)
{
 mycall <- match.call()
 ## if any objects that are not formals of MRSP.fit are passed via ... , we need
 ## to explicitly assign the corresponding objects:
 dotlist <- list(...)
 for(i in seq_len(length(dotlist))) assign(names(dotlist)[i], dotlist[[i]])
 if(is.expression(model)) model <- eval(model)

 
 ## replace all entries of mycall with their actual value. this is useful, for
 ## example, for bootstrapping or predict methods (and so on...)
 if(control@expandcall){
  excl <- c(1, which(names(mycall) %in% c("dat", "model")))
  if(is.null(mycall$model)){
   mycall$model <- expression(multinomlogit())
   excl <- c(excl, length(mycall))
  }else{
   mycall[[which(names(mycall) == "model")]] <- switch(model@name,
    "Multinomial Logit Model" = expression(multinomlogit()),
    "Sequential Logit Model" = expression(sequentiallogit()),
    "CUB Binomial Logit Model" = expression(CUBbinomiallogit()),
    "OLS Regression" = expression(OLSreg())
   )
  }
  mycall[-excl] <- mget(names(mycall[-excl]), envir = as.environment(-1))
 }
 
 ## if more than one lambda was supplied, fix this by going back to
 ## MRSP.fit(..., signature(lambda=list)):
 if(length(lambda) > 1){mycall$lambda <- list(lambda); out <- eval(mycall); return(out)}

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

 ## if the first response category is chosen as reference, this can cause 
 ## problems on many occasions. therefore, we rearrange stuff in this case.
 ## really? probably not! commented out for now!
 #if(constr == 1){
 # dat$y[,c(1, K)] <- dat$y[, c(K, 1)]
 # colnames(dat$y)[c(1, K)] <- colnames(dat$y)[c(K, 1)]
 # dat$V[c(1, K)] <- dat$V[c(K, 1)]
 # if(!(missing(mlfit) | is.null(mlfit)) | !(missing(coef.init) | is.null(coef.init)) |
 #    !(missing(coef.stand.init) | is.null(coef.stand.init)) | 
 #    !(missing(coef.pretres.init) | is.null(coef.pretres.init)) |
 #    !(missing(penweights) | is.null(penweights))  ){
 #     stop("something went completely wrong in the passing on of arguments for the case of 'constr == 1'. please use a response category different from the first one as reference.")
 # }  
 #} 

 if(!(missing(grpindex) | is.null(grpindex))){
  if(!is.list(grpindex)){grpindex <- list(grpindex)}
 }
 if(!(missing(penindex) | is.null(penindex))){
  if(!is.list(penindex)){penindex <- list(penindex)}
 }

 if((missing(mlfit) | is.null(mlfit)) & control@doMLfit == T){
  mlfit <- getmlfit(dat = dat, coef.init = coef.init, coef.stand.init = coef.stand.init,
                              offset = offset, weights = weights, grpindex = grpindex,
                              penindex = penindex, lambda = lambda, lambdaR = lambdaR,
                              lambdaF = lambdaF, gamma = gamma,
                              psi = psi, indg = indg, indcs = indcs, model = model,
                              control = control, fista.control = fista.control, constr = constr,
                              Proximal.control = Proximal.control, Proximal.args = Proximal.args,
                              adaptive = adaptive, threshold = threshold, refit = refit,
                              nonneg = nonneg, ...)
 }                             

 if(missing(penweights) | is.null(penweights)){
  penweights <- getpenweights(dat = dat, coef.init = coef.init, coef.stand.init = coef.stand.init,
                              offset = offset, weights = weights, grpindex = grpindex,
                              penindex = penindex, lambda = lambda, lambdaR = lambdaR,
                              lambdaF = lambdaF, gamma = gamma,
                              psi = psi, indg = indg, indcs = indcs, model = model, 
                              control = control, fista.control = fista.control, constr = constr,
                              Proximal.control = Proximal.control, Proximal.args = Proximal.args,
                              adaptive = adaptive, threshold = threshold, refit = refit,
                              fusion = fusion, mlfit = mlfit, nonneg = nonneg, ...)
 }
 

 if(is.list(lambda) | length(lambda) != 1 | lambda < 0)
   stop("something went wrong with lambda")
   
 tuning <- list(lambda = lambda,
                gamma  = gamma,
                psi    = psi,
                lambdaR = lambdaR,
                lambdaF = lambdaF)


 ## handling the intercept
 which.intercept <- which(apply(dat$x, 2, function(u) diff(range(u)) < .Machine$double.eps ^ 0.5))
 has.intercept <- !(length(which.intercept) == 0)
 if((length(which.intercept) == 1) && which.intercept != 1){
  dat$x[,c(1, which.intercept)] <- dat$x[,c(which.intercept, 1)] 
  penindex[[1]][c(1,which.intercept)] <- penindex[[1]][c(which.intercept,1)]
  grpindex[[1]][c(1,which.intercept)] <- grpindex[[1]][c(which.intercept,1)]
  if(!all(coef[[1]] == 0))
  coef[[1]][,c(1,which.intercept)] <- coef[[1]][,c(which.intercept,1)]
  which.intercept <- 1
 }
 if(length(which.intercept) > 1)
   stop("more than one intercept column specified")

 ## are we using an ordinal regression model?
 isordinal <- model@name %in% c("Cumulative Logit Model", "Sequential Logit Model")
 if(isordinal){
  if(!is.list(penindex)){
   if(is.null(penindex)) stop("penindex must be supplied in ordinal regression")
   penindex <- list(penindex)
  }
  sl.indg <- which(penindex[[1]] %in% c(4, 40, 41))
  sl.indcs <- which(penindex[[1]] %in% c(1, 10, 11, 12, 13, 14))
 }else{
  sl.indg <- NULL
  sl.indcs <- NULL
 } 
 if(isordinal && (constr != "none")) stop("You must use no constraint in ordinal models")

 ## now fill in the potentially missing arguments:
 if(missing(coef.init) | is.null(coef.init)){
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
  if(K == 1){
   coef[[1]][,1] <- weighted.mean(dat$y, w = weights)
  }
 }else{coef <- coef.init}

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

 if(missing(coef.stand.init) | is.null(coef.stand.init)){
  coef.stand <- coef
 }else{coef.stand <- coef.stand.init}
 if(missing(coef.pretres.init) | is.null(coef.pretres.init)){
  coef.pretres <- coef.stand
 }else{coef.pretres <- coef.pretres.init}
 if(!is.list(coef.pretres)){coef.pretres <- list(coef.pretres)} 
 if(!is.list(coef.stand)){coef.stand <- list(coef.stand)}
 if(!is.list(coef)){coef <- list(coef)}

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
   
 if(!all(weights == rep(1, nrow(dat$y))) && (control@standardize == T))
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
   stop("y values outside the appropriate range for the specified model class!")
   
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
    
 if(length(tuning) != 5)
   stop("for the given data, the tuning object needs to be a list of length 4")
     
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
 
 ## if spGFL is used, prepare some objects for it:
 if(any.spGFLpen.x && !isordinal) stop("GFL penalty is currently only supported for ordinal models")

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

 ## extract control information
 max.iter <- control@max.iter
 rel.tol <- control@rel.tol
 standardize <- control@standardize
 keepdat <- control@keepdat
 ridgestabil <- control@ridgestabil
 ridgestabilrf <- control@ridgestabilrf
 lambdastabil <- control@lambdastabil
 trace <- control@trace
 cut.to.rel <- control@cut.to.rel
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

 #####################################################################################
 ## orthonormalization of the original observations. should be done in the wrap-around
 ## call to MRSP.fit, which applies MRSP.fit to the sequence of lambda-
 ## values. 
 ##### note: a correct backtransformation from the coefficients belonging to
 ## standardized predictors to the corresponding coefficients belonging to the
 ## original predictors is only possible if an intercept column is present!
 
 warnold <- options("warn")
 if(warnold < 2)
  options(warn = 1)
 if(!has.intercept & standardize){
  warning()
  cat(" Standardization of predictors was requested while no intercept column could be found.",'\n',
  "For this case, backtransformation of coefficients to the original scale is not possible, so that output for the standardized predictors is returned instead.",'\n',
  "Be careful!",'\n')
 }
 if(warnold < 2)
  options(warn = as.numeric(warnold))
 
 ## orthonormalization of observations
 x.original <- dat$x
 if(hasV) V.original <- dat$V
    
 if(standardize){
  sweights <- weights / sum(weights)
  x.centered <- dat$x
  mean.x <- apply(dat$x[,,drop=F], 2, function(u) weighted.mean(u, weights))
  if(has.intercept){
   x.centered[,-1] <- sweep(dat$x[,-1, drop=F], 2, mean.x[-1])
  }else{x.centered <- sweep(dat$x[,, drop=F], 2, mean.x)} 

  xw <- x.centered
  scale.pen.x <- list()
  scale.notpen.x <- NULL
  ## if there are no unpenalized groups, skip the following:

  if(length(notpen.which.x) > 0){
   if(has.intercept){
    scale.notpen.x <- sqrt(drop(sweights %*% (x.centered[, notpen.which.x[-1]]^2)))
    xw[,notpen.which.x[-1]] <- scale(x.centered[, notpen.which.x[-1]], FALSE, scale.notpen.x)
   }else{
    scale.notpen.x <- sqrt(drop(sweights %*% (x.centered[, notpen.which.x]^2)))
    xw[,notpen.which.x] <- scale(x.centered[, notpen.which.x], FALSE, scale.notpen.x)     
   }
  }
  
  for(j in dict.pen.x){
   j.which <- which(grpindex[[1]] == j)
   decomp.x <- qr(sqrt(weights) * x.centered[, j.which])
   if(decomp.x$rank < length(j.which)) ## warn if block doesnt have full rank
     stop("Block belonging to columns ", paste(j.which), "does not have full rank!")
   scale.pen.x[[j]] <- qr.R(decomp.x) * 1/sqrt(wnobs)
   xw[, j.which] <- qr.Q(decomp.x) * sqrt(wnobs)
  }
  
  dat$x <- xw  
  
  if(hasV){
   eweights <- rep(weights, times=K)
   seweights <- rep(sweights, times=K)
   ## extract V into a matrix of size (nobs * K) x L:
   Vs <- do.call(rbind, dat$V)
   mean.V <- apply(Vs, 2, function(u) weighted.mean(u, eweights))
   V.centered <- sweep(Vs, 2, mean.V)

   Vw <- V.centered
   scale.pen.V <- list()
   scale.notpen.V <- NULL
   if(length(notpen.which.V) > 0){
    scale.notpen.V <- sqrt(drop(seweights %*% (V.centered[, notpen.which.V]^2)))
    Vw[, notpen.which.V] <- scale(V.centered[, notpen.which.V], FALSE, scale.notpen.V)
   }
   
   for(j in dict.pen.V){
    j.which <- which(grpindex[[2]] == j)
    decomp.V <- qr(sqrt(eweights) * V.centered[, j.which])
    if(decomp.V$rank < length(j.which))
      stop("Block belonging to columns ", paste(j.which), "in V does not have full rank!")
    scale.pen.V[[j]] <- qr.R(decomp.V) * 1/sqrt(wnobs * K)
    Vw[, j.which] <- qr.Q(decomp.V) * sqrt(wnobs * K)
   }
   
   for(i in seq_len(K)){dat$V[[i]] <- Vw[seq(nobs*i-nobs+1,nobs*i),, drop=F]}                                      
  }
  
  if(threshold != F){
   coef <- coef.pretres
  }else{  
   coef <- coef.stand
  }
 }     

 ## assign the proxy class to the coef object
 class(coef) <- "MRSP.coef"

 eta <- updateEta(dat=dat, coef=coef, offset=offset)
 mu <- invlink(eta)
                    
 ## prepare the tuning object for the call to fista. note that only true "lambda"
 ## parameters may be left in tuning when fista calls it. gamma is for spgrplasso
 ## psi is for global vs cat-specific
 ## tuning[[1]] is for the group lasso penalty,
 ## tuning[[2]] for the ordinary lasso
 ## tuning[[3]] is for the category-specific predictors with global effects or
 ## for cat-spec predictors with cat-spec effects and group lasso penalty
 ## tuning[[6]] is for the category-specific predictors with cat-spec and 
 ## ordinary lasso penalty
 ## tuning[[4]] is for ridge penalties
 ## tuning[[5]] is for a ridge penalty that may be used for all parameters which
 ## otherwise would be unpenalized (including the intercept), which can be tried
 ## in order to stabilize the model if divergence of parameters occurs.
 ## tuning[[7]] is for the GFL penalty and, generally speaking, fusion penalties.
 tuning.old <- tuning
 lambda <- tuning$lambda
 gamma <- tuning$gamma
 psi <- tuning$psi
 lambdaR <- tuning$lambdaR
 lambdaF <- tuning$lambdaF
 tuning <- list()
 tuning[[1]] <- lambda *  psi
 tuning[[2]] <- lambda * gamma *psi
 tuning[[4]] <- lambdaR
 if(hasV){
  tuning[[3]] <- lambda * (2 - psi)
  tuning[[6]] <- lambda * (2 - psi) * gamma
 }
 tuning[[7]] <- lambdaF
 lmod <- K*P
 if(hasV) lmod <- lmod + L
 if(nobs <= lmod) lambdastabil <- 1*lambdastabil   
 tuning[[5]] <- 4*sqrt(2*log(lmod)/nobs) * lambdastabil

 ## list with arguments for the calls to fistaProximal and the other fistaGenerics
 if(missing(Proximal.args) | is.null(Proximal.args)){
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
 }
 
 logl.old <- loglik(y=dat$y, mu=mu, weights=weights)
 pen.old <- penalty(coef=coef, tuning=tuning, penweights=penweights,
                    penindex=penindex, grpindex=grpindex, Proximal.args=Proximal.args)
 fn.val.old <- -logl.old + pen.old
 coef.old <- coef

 ## create the fista.control argument
 if(missing(fista.control) | is.null(fista.control)){
  fista.control <- fista.control(rel.tol = rel.tol, max.iter = max.iter)
 }else{fista.control@rel.tol <- rel.tol; fista.control@max.iter <- max.iter} 

 ## now the call to fista where the core fitting takes place.
 updated <- fista(coef.init=coef, dat=dat, offset=offset, eta.init=eta,
                  mu.init=mu, weights=weights, model=model, control=fista.control,
                  Proximal.args=Proximal.args, Proximal.control=Proximal.control,
                  grpindex=grpindex, penindex=penindex, tuning=tuning,
                  penweights.init=penweights)
                 
 coef <- updated$coef
 penweights <- updated$penweights
 eta <- updated$eta
 mu <- updated$mu
 logl <- updated$loglikval
 logl.approx <- 0
 pen <- updated$penval
 iter.count <- updated$iter.count
 best.iter <- updated$best.iter
 d.fn.updated <- updated$d.fn
 fn.val <- updated$fn.val
 ridgestabil <- updated$ridgestabil
 control@ridgestabilrf <- updated$ridgestabilrf
 
 ## some checks that nothing has gone completely wrong
 if(is.na(pen) | is.na(logl) | is.na(logl.approx) | is.na(fn.val) | is.na(d.fn.updated))
   warning(paste("pen, logl, logl.approx, fn.val, d.fn or several of them were NA"))
 if(class(coef) != "MRSP.coef")
   stop("the coef object returned by fista did not retain the correct class")
 
 ## convergence checks
 d.fn <- (fn.val.old - fn.val) / (0.001 + abs(fn.val))

 if(sign(d.fn) != sign(d.fn.updated)){
  if(trace){
   warning(paste("the value of the penalized loglikelihood has worsened"))
  }
 } 
 if(abs(d.fn.updated) > sqrt(rel.tol)){
  if(trace){
   print(paste("d.fn.updated:"))
   print(d.fn.updated)
   print(paste("rel.tol:"))
   print(rel.tol)
   warning(paste("the fista algorithm might have converged earlier than wanted"))
  }
 }  
 d.par <- max(mapply(function(x,y){max(abs(x-y)/(0.01+abs(y)))}, coef.old, coef))
 
 ## now prepare for the output
 #############################
 ## cut values of the standardized coefficients which are smaller than
 ## control@cut.to.rel. if it's set to 0, there is no cutting.
 ## note that this has to be done before the retransformation. otherwise, the 
 ## coefs of predictors with large variance might get cut although they are
 ## highly significant!
 coef <- lapply(coef, function(u){u[which(abs(u) <= cut.to.rel)] <- 0; return(u)})

 ## save the coefficients before thresholding so that the next model in the list
 ## of lambdas gets computed as fast as possible.
 coef.pretres <- coef
 
 ## perform thresholding
 if(is.function(control@numerictresh) & threshold == F){
  threshold <- control@numerictresh(nobs)
  if(refit == T | refit == "L") threshold <- threshold / 2
 } 
 if(threshold == T){ 
  if(hasV) dimmod <- P + L else dimmod <- P
  threshold <- 0.5*sqrt(2*log(dimmod)/nrow(dat$y))
 } 
 if(is.numeric(threshold) & length(threshold) == 1){
  for(j in dict.pen.x){
   j.which <- which(grpindex[[1]] == j)
   if(penindex[[1]][j.which[1]] %in% c(1, 13)){
    if(constr == "none"){
     thresholdval <- l2norminv(coef[[1]][, j.which])
    }else{
     thresholdval <- l2norminv.x(coef[[1]][, j.which])
    }       
    if(thresholdval < threshold) coef[[1]][, j.which] <- 0
   }else if(penindex[[1]][j.which[1]] %in% c(12)){
    if(constr == "none"){
     thresholdval <- apply(coef[[1]][, j.which, drop=F], 1, l2norminv)
    }else{
     thresholdval <- apply(coef[[1]][, j.which, drop=F], 1, l2norminv.x)
    } 
    whichnull <- which(thresholdval < threshold)
    coef[[1]][whichnull, j.which] <- 0
   }else if(penindex[[1]][j.which[1]] %in% c(11)){
    if(constr == "none"){
     thresholdval1 <- l2norminv(coef[[1]][, j.which])
     thresholdval2 <- apply(coef[[1]][, j.which, drop=F], 1, l2norminv)
    }else{
     thresholdval1 <- l2norminv.x(coef[[1]][, j.which])
     thresholdval2 <- apply(coef[[1]][, j.which, drop=F], 1, l2norminv.x)
    }                        
    if(thresholdval1 < threshold) coef[[1]][, j.which] <- 0
    whichnull <- which(thresholdval2 < threshold)
    coef[[1]][whichnull, j.which] <- 0
   }else if(penindex[[1]][j.which[1]] %in% c(14)){
    if((constr == "none" | constr == "symmetric")){
     thresholdval1 <- l2norminv(coef[[1]][, j.which])
     thresholdval2 <- apply(coef[[1]][, j.which, drop=F], 1, l2norminv)
     thresholdval3 <- l2norminv(R%*%coef[[1]][, j.which])
    }else{ stop("GFL penalty can only be used in models with no or a symmetric side constraint for now") }
    if(thresholdval1 < threshold) coef[[1]][, j.which] <- 0
    if(thresholdval3 < threshold){
     coef[[1]][, j.which] <- matrix(rep(colMeans(coef[[1]][, j.which]), each=K), ncol=length(j.which))
    }
    whichnull <- which(thresholdval2 < threshold)
    coef[[1]][whichnull, j.which] <- 0
   }else if(penindex[[1]][j.which[1]] %in% c(17)){
    if(constr == "none" | constr == "symmetric"){
     thresholdval1 <- l2norminv(coef[[1]][, j.which])
     thresholdval2 <- apply(coef[[1]][, j.which, drop=F], 1, l2norminv)
    }else{ stop("Sparse group lasso with L1-fusion can only be used in models without side constraint for now") }
    if(thresholdval1 < threshold) coef[[1]][, j.which] <- 0
    clusterj <- sapply(seq(nrow(coef[[1]][, j.which])), list)
    for(s in seq(nrow(coef[[1]][, j.which]))){
     for(l in seq(j, nrow(coef[[1]][, j.which]))){
      if(l2norminv(coef[[1]][s, j.which] - coef[[1]][l, j.which]) < threshold){
       s.ind <- which(clusterj %in% coef[[1]][s, j.which])
       l.ind <- which(clusterj %in% coef[[1]][l, j.which])
       clusterj[[s.ind]] <- c(clusterj[[s.ind]], clusterj[[l.ind]])
       clusterj[[l.ind]] <- NULL
       coef[[1]][clusterj[[s.ind]], j.which] <- matrix(colMeans(coef[[1]][clusterj[[s.ind]], j.which]), nrow=length(clusterj[[s.ind]]), ncol=length(j.which), byrow=T)
      }
     }
    }
    whichnull <- which(thresholdval2 < threshold)
    coef[[1]][whichnull, j.which] <- 0
   }else if(penindex[[1]][j.which[1]] %in% c(18)){
    if(constr == "none"){
     for(i in seq(nrow(coef[[1]][, j.which]))){
      clusteri <- sapply(seq(ncol(coef[[1]][, j.which])), list)                 ## first every coefficient forms its own cluster
      for(s in seq(ncol(coef[[1]][, j.which])-1)){
       for(l in seq(s+1, ncol(coef[[1]][, j.which]))){
        if(abs(coef[[1]][i, j.which][,s] - coef[[1]][i, j.which][,l]) < threshold){
         s.ind <- which(sapply(clusteri, function(x) s %in% x))   #which(coef[[1]][i, j.which] %in% coef[[1]][i, j.which][,s])
         l.ind <- which(sapply(clusteri, function(x) l %in% x))   #which(coef[[1]][i, j.which] %in% coef[[1]][i, j.which][,l])
         clusteri[[s.ind]] <- unique(c(clusteri[[s.ind]], clusteri[[l.ind]]))
         if(s.ind != l.ind){
          clusteri[[l.ind]] <- NULL
          if(s.ind > l.ind) s.ind <- s.ind - 1
         }
         coef[[1]][i, j.which][, clusteri[[s.ind]]] <- mean(coef[[1]][i, j.which][, clusteri[[s.ind]]])
        }
       }
      }
     }
    }else{
     stop("Rowwise L1-fusion can only be used in models without side constraint for now")
    }
   }
   if(penindex[[1]][j.which[1]] %in% c(4, 41)){
    thresholdval <- l2norminv(coef[[1]][1, j.which])                                        
    if(thresholdval < threshold) coef[[1]][, j.which] <- 0
   }      
  }
  if(hasV){
   for(j in dict.pen.V){
    j.which <- which(grpindex[[2]] == j) 
    if(penindex[[2]][j.which[1]] %in% c(2, 21)){
     thresholdval <- l2norminv(coef[[2]][1, j.which])                                        
     if(thresholdval < threshold) coef[[2]][, j.which] <- 0
    }else if(penindex[[2]][j.which[1]] %in% c(3, 33)){
     thresholdval <- l2norminv(coef[[2]][, j.which])
     if(thresholdval < threshold) coef[[2]][, j.which] <- 0
    }else if(penindex[[2]][j.which[1]] %in% c(32)){
     thresholdval <- apply(coef[[2]][, j.which, drop=F], 1, l2norminv)
     whichnull <- which(thresholdval < threshold)
     coef[[2]][whichnull, j.which] <- 0
    }else if(penindex[[2]][j.which[1]] %in% c(31)){
     thresholdval1 <- l2norminv(coef[[2]][, j.which])
     thresholdval2 <- apply(coef[[2]][, j.which, drop=F], 1, l2norminv)
     if(thresholdval1 < threshold) coef[[2]][, j.which] <- 0
     whichnull <- which(thresholdval2 < threshold)
     coef[[2]][whichnull, j.which] <- 0
    }                            
   }
  }
 }


 ##### degrees of freedom
 ## note that the numeric thresholding performed in the previous lines means
 ## that the symmetric side constraint, if used, might not be fully satisfied
 ## for the thresholded coefs. therefore, subtracting from the edf in order to 
 ## account for the side constraint might in rare instances make the following
 ## edf formulas return negative values - thus, we use a max(edf, 0)-construct.
 ## for the general edf formula, see Tutz, Pößnecker, Uhlmann (2014) or
 ## Yuan & Lin (2006). 
 df <- 0
 
 if(control@computeDF){
  if(any.notpen.x){
   df <- df + max(length(c(coef[[1]][,which.notpen.x])) - ifelse(constr == "none", 0, length(which.notpen.x)), 0)
  }

  if(any.ridgepen.x){
   #add <- sum(diag(dat$x%*%solve(crossprod(dat$x) + lambdaR * diag(ncol(dat$x)))%*%t(dat$x)))   ## <- that only works for linear models..
   add <- max(length((coef[[1]][,which.ridgepen.x])) - ifelse(constr == "none", 0, length(which.ridgepen.x)), 0) / (1 + lambdaR)
   df <- df + add #max(length((coef[[1]][,which.ridgepen.x])) - ifelse(constr == "none", 0, length(which.ridgepen.x)), 0)
  }
  if(any.grouppen.x){
   for(j in groups.grouppen.x){
    j.which <- which(grpindex[[1]] == j)
    coefj <- coef[[1]][,j.which]
    df <- df + max(ifelse(l2normraw(coefj) > 0, 1 +
    (length(c(coefj)) - 1 - ifelse(constr == "none", 0, ncol(coefj))) *
    l2normraw(coefj) / l2normraw(mlfit$coef.stand[[1]][,j.which]), 0), 0)
   }
  }
 
  if(any.spgrppen.x){
   for(j in groups.spgrppen.x){
    j.which <- which(grpindex[[1]] == j)
    coefj <- coef[[1]][,j.which]
    df <- df + max((0.5 * ifelse(l2normraw(coefj) > 0, 1 +
    (length(c(coefj)) - 1 - ifelse(constr == "none", 0, ncol(coefj))) *
    l2normraw(coefj) / l2normraw(mlfit$coef.stand[[1]][,j.which]), 0)   +
    0.5 * gamma * sum(sapply(1:nrow(coefj), function(i){ifelse(l2normraw(coefj[i,]) > 0, 1 +
    ifelse(length(c(coefj[i,])) == 1, 0, (length(c(coefj[i,])) - 1) * l2normraw(coefj[i,]) / l2normraw(mlfit$coef.stand[[1]][i,j.which])), 0)}))
    - ifelse(constr == "symmetric", ncol(coefj), 0)) / (0.5 + 0.5 * gamma), 0)
   }
  }
 
  if(any.lassopen.x){
   for(j in groups.lassopen.x){
    j.which <- which(grpindex[[1]] == j)
    coefj <- coef[[1]][,j.which]
    df <- df + max(sum(sapply(1:nrow(coefj), function(i){ifelse(l2normraw(coefj[i,]) > 0, 1 +
    ifelse(length(c(coefj[i,])) == 1, 0, (length(c(coefj[i,])) - 1) * l2normraw(coefj[i,]) / l2normraw(mlfit$coef.stand[[1]][i,j.which])), 0)}))
    - ifelse(constr == "symmetric", ncol(coefj), 0), 0)
   }
  }

  if(any.spGFLpen.x){
   for(j in groups.spGFLpen.x){
    j.which <- which(grpindex[[1]] == j)
    coefj <- coef[[1]][,j.which]
    #if(all(apply(coefj, 2, function(u) all(diff(range(u)) < .Machine$double.eps ^ 0.5)))){
    # coefj <- coef[[1]][1,j.which]
    # add <- ifelse(l2normraw(coefj) > 0, 1 + ifelse(length(c(coefj)) == 1, 0, (length(c(coefj)) - 1) * l2normraw(coefj) / l2normraw(mlfit$coef.stand[[1]][1,j.which])), 0)
    #}else{
     diffnormratio <- ifelse(lambdaF > 0, (ifelse(l2normraw(R%*%coefj) > 0, 1, 0) + (K-2) * ifelse(l2normraw(coefj) > 0, l2normraw(mlfit$coef.stand[[1]][,j.which]) / l2normraw(coefj), 1) *
                                           l2normraw(R%*%coefj) / l2normraw(R%*%mlfit$coef.stand[[1]][,j.which])) / (K-1), 1)

     lassoadd <-  0.5 * gamma * sum(sapply(1:nrow(coefj), function(i){ifelse(l2normraw(coefj[i,]) > 0, 1 +
      ifelse(length(c(coefj[i,])) == 1, 0, (length(c(coefj[i,])) - 1) * l2normraw(coefj[i,]) / l2normraw(mlfit$coef.stand[[1]][i,j.which])), 0)}))

     grplassoadd <- 0.5 * ifelse(l2normraw(coefj) > 0, 1 +
      (length(c(coefj)) - 1 - ifelse(constr == "none", 0, ncol(coefj)))*
      l2normraw(coefj) / l2normraw(mlfit$coef.stand[[1]][,j.which])*diffnormratio, 0)
    
     add <- max(( grplassoadd + 1/K * lassoadd + (K-1)/K * diffnormratio * lassoadd
      - ifelse(constr == "symmetric", ncol(coefj), 0)) / (0.5 + 0.5 * gamma), 0)
    #}
    df <- df + add
   }
  }

  if(any.Psplinepen.x){
   for(j in groups.Psplinepen.x){
    j.which <- which(grpindex[[1]] == j)
    coefj <- coef[[1]][,j.which]
    add <- 0
    for(i in seq(K)){
     add <- add + ifelse(constr != i, 1 + (nrow(coefj) - 1) * crossprod(t(tD[[j]])%*%coefj[i,,drop=T]) / crossprod(t(tD[[j]])%*%mlfit$coef.stand[[1]][i,j.which,drop=T]), 0)
    }
    df <- df + add
   }
  }
  
  if(any.colPsplinepen.x){
   for(j in groups.colPsplinepen.x){
    j.which <- which(grpindex[[1]] == j)
    coefj <- coef[[1]][,j.which]
    df <- df + ncol(coefj) + ncol(coefj)*(nrow(coefj) - 1 - ifelse(constr == "none", 0, 1)) * crossprod(R%*%coefj)/crossprod(R%*%mlfit$coef.stand[[1]][,j.which])
   }
  }
  
  if(any.spFGLpen.x){
   for(j in groups.spFGLpen.x){
    j.which <- which(grpindex[[1]] == j)
    coefj <- coef[[1]][,j.which]
    diffnormratio <- ifelse(lambdaF > 0, (ifelse(l2normraw(R%*%coefj) > 0, 1, 0) + (K-2) * ifelse(l2normraw(coefj) > 0, l2normraw(mlfit$coef.stand[[1]][,j.which]) / l2normraw(coefj), 1) *
                                          l2normraw(R%*%coefj) / l2normraw(R%*%mlfit$coef.stand[[1]][,j.which])) / (K-1), 1)

    coefju <- unique(coefj)
    clustertable <- matrix(c(seq(nrow(coefj)), rep(0, nrow(coefj))), ncol=2)
    for(i in seq(nrow(coefj))){
     for(l in seq(nrow(coefju))){
      if(all(coefj[i,] == coefju[l,])) clustertable[i,2] <- l
     }
    }
      
    lassoadd <-  0.5 * gamma * sum(sapply(1:nrow(coefj), function(i){ifelse(l2normraw(coefj[i,]) > 0, 1 +
     ifelse(length(c(coefj[i,])) == 1, 0, (length(c(coefj[i,])) - 1) * l2normraw(coefj[i,]) / l2normraw(mlfit$coef.stand[[1]][i,j.which])), 0) / sum(clustertable[,2] == clustertable[i,2])}))

    grplassoadd <- 0.5 * ifelse(l2normraw(coefj) > 0, 1 +
     (length(c(coefj)) - 1 - ifelse(constr == "none", 0, ncol(coefj)))*
     l2normraw(coefj) / l2normraw(mlfit$coef.stand[[1]][,j.which])*diffnormratio, 0)

    add <- max(( grplassoadd + lassoadd
     - ifelse(constr == "symmetric", ncol(coefj), 0)) / (0.5 + 0.5 * gamma), 0)
    df <- df + add
   }
  }

  if(any.rowFLpen.x){
   for(j in groups.rowFLpen.x){
    j.which <- which(grpindex[[1]] == j)
    coefj <- coef[[1]][,j.which]
    df <- df + sum(apply(coefj, 1, function(u){length(unique(u[u!=0]))}))
   }
  }
 
  if(hasV){
   if(any.globalunpen.V){
    df <- df + length(c(coef[[2]][1,which.globalunpen.V]))
   }

   if(any.globalridgepen.V){
    df <- df + length(c(coef[[2]][1,which.globalridgepen.V]))
   }
    
   if(any.globalpen.V){
    for(j in groups.globalpen.V){
     j.which <- which(grpindex[[2]] == j)
     coefj <- coef[[2]][1,j.which]
     df <- df + ifelse(l2normraw(coefj) > 0, 1 +
     ifelse(length(c(coefj)) == 1, 0, (length(c(coefj)) - 1) * l2normraw(coefj) / l2normraw(mlfit$coef.stand[[2]][1,j.which])), 0)
    }
   }
  
   if(any.catpen.V){
    for(j in groups.catpen.V){
     j.which <- which(grpindex[[2]] == j)
     coefj <- coef[[2]][,j.which]
     df <- df + ifelse(l2normraw(coefj) > 0, 1 +
     (length(c(coefj)) - 1) * l2normraw(coefj) / l2normraw(mlfit$coef.stand[[2]][,j.which]), 0)
    }
   }

   if(any.catunpen.V){
    df <- df + sum(coef[[2]][,which.catunpen.V] != 0)
   }

   if(any.catridgepen.V){
    df <- df + sum(coef[[2]][,which.catridgepen.V] != 0)
   }
  
   if(any.catspgrppen.V){
    for(j in groups.catspgrppen.V){
     j.which <- which(grpindex[[2]] == j)
     coefj <- coef[[2]][,j.which]
     df <- df + (0.5 * ifelse(l2normraw(coefj) > 0, 1 +
     (length(c(coefj)) - 1) *
     l2normraw(coefj) / l2normraw(mlfit$coef.stand[[2]][,j.which]), 0)   +
     0.5 * gamma * sum(sapply(1:nrow(coefj), function(i){ifelse(l2normraw(coefj[i,]) > 0, 1 +
     ifelse(length(c(coefj[i,])) == 1, 0, (length(c(coefj[i,])) - 1) * l2normraw(coefj[i,]) / l2normraw(mlfit$coef.stand[[2]][i,j.which])), 0)}))
     ) / (0.5 + 0.5 * gamma)
    }
   }
  
   if(any.catlassopen.V){
    for(j in groups.catlassopen.V){
     j.which <- which(grpindex[[2]] == j)
     coefj <- coef[[2]][,j.which]
     df <- df + sum(sapply(1:nrow(coefj), function(i){ifelse(l2normraw(coefj[i,]) > 0, 1 +
     ifelse(length(c(coefj[i,])) == 1, 0, (length(c(coefj[i,])) - 1) * l2normraw(coefj[i,]) / l2normraw(mlfit$coef.stand[[2]][i,j.which])), 0)}))
    }
   }
  } # end 'hasV'
 
  if(isordinal){
   if(any.globalunpen.o){
    df <- df + length(c(coef[[1]][1,which.globalunpen.o]))
   }

   if(any.globalridgepen.o){
    df <- df + length(c(coef[[1]][1,which.globalridgepen.o]))
   }
    
   if(any.globalpen.o){
    for(j in groups.globalpen.o){
     j.which <- which(grpindex[[1]] == j)
     coefj <- coef[[1]][1,j.which]
     df <- df + ifelse(l2normraw(coefj) > 0, 1 +
    ifelse(length(c(coefj)) == 1, 0, (length(c(coefj)) - 1) * l2normraw(coefj) / l2normraw(mlfit$coef.stand[[1]][1,j.which])), 0)
    }
   }
  } # end 'isordinal'
 }else{
  df <- sum(v.allequal(Reduce("c", coef), 0))
 }                                               # end if(control@computeDF)

 ## now retransform the coefficients back to the original scale if 
 ## standardization and centering were used.
 coef.stand <- coef

 if(control@backtransf & has.intercept){
  if(length(notpen.which.x) > 1 && (which.intercept %in% notpen.which.x)){
   coef[[1]][, notpen.which.x[-1]] <- scale(coef[[1]][, notpen.which.x[-1]], FALSE, scale.notpen.x)
  }else if(length(notpen.which.x) > 0 && !(which.intercept %in% notpen.which.x)){
   coef[[1]][, notpen.which.x] <- scale(coef[[1]][, notpen.which.x], FALSE, scale.notpen.x)
  }


  for(j in dict.pen.x){
   j.which <- which(grpindex[[1]] == j)
   if(any(coef[[1]][,j.which] != 0)){
    coef[[1]][, j.which] <- t(apply(coef[[1]][, j.which, drop=F], 1, function(u) solve(scale.pen.x[[j]], u)))
   }
  }
  
  if(hasV){
   if(length(notpen.which.V) > 0){
    coef[[2]][,notpen.which.V] <- scale(coef[[2]][,notpen.which.V], FALSE, scale.notpen.V)
   }
  
   for(j in dict.pen.V){
    j.which <- which(grpindex[[2]] == j)
    if(any(coef[[2]][,j.which] != 0)){                                    ## fixme: for cs-g variables with more than one column, the retransformation might break up the "equal coefficients in all rows of coef[[2]]"-constraint
     coef[[2]][,j.which] <- t(apply(coef[[2]][,j.which, drop=F], 1, function(u) solve(scale.pen.V[[j]], u)))
    }
   }  
  }   
  
  coefsum.x <- rowSums(sweep(coef[[1]][,-1, drop=F], 2, mean.x[-1], FUN="*"))
  coef[[1]][,1] <- coef[[1]][,1] - coefsum.x
  if(hasV){
   coefsum.V <- rowSums(sweep(coef[[2]], 2, mean.V, FUN="*"))                                            
   coef[[1]][,1] <- coef[[1]][,1] - coefsum.V
  }
  ## sometimes, there are spill-over effects on the intercept for the reference
  ## category, making it unequal to zero.
  if(is.numeric(constr)){
   if(abs(coef[[1]][constr,1]) > 1e-4){
    if(penindex[[1]][1] %in% c(10, 13)){
     coef[[1]][,1] <- coef[[1]][,1] - rep(coef[[1]][constr,1], length(coef[[1]][,1]))
    }else{
     coef[[1]][constr,1] <- 0
    }
   }else{
    coef[[1]][constr,1] <- 0
   }    
  }
 } 
 ## ensure the symmetric side constraint. (note for future self and curious 
 ## readers of this source code: the following should not work in theory, but
 ## extensive tests have shown that it does work in practice. dont ask me why..)
 if(constr == "symmetric"){
  coef[[1]] <- apply(coef[[1]], 2, function(u){u[u!=0] <- u[u!=0] - mean(u[u!=0]); u})
 } 

 ## now, we have to check the signs. sometimes the qr decomposition turns around
 ## the signs of estimated coefficient groups. the backtransformed coefficients
 ## belonging to the original observations always have the correct sign, but 
 ## coef.stand and coef.pretres sometimes dont. therefore:
 if(standardize | control@backtransf){
  coef.stand <- Map(function(x,y){x*(sign(x)*sign(y))}, coef.stand, coef)
  coef.pretres <- Map(function(x,y){x*(sign(x)*sign(y))}, coef.pretres, coef)
 }       

 ## prepare the rest of the output
 if(standardize){
  x.stand <- dat$x
  dat$x <- x.original
  if(hasV){
   V.stand <- dat$V
   dat$V <- V.original
  } 
 }
 
 ## list of active parameter groups/predictors
 guessed.active <- list()
 guessed.active[[1]] <- numeric(length(unique(grpindex[[1]])))
 guessed.active.coef <- list()
 guessed.active.coef[[1]] <- matrix(nrow = K, ncol = length(unique(grpindex[[1]])))
 guessed.active.groupdiff <- list()
 guessed.active.groupdiff[[1]] <- numeric(length(unique(grpindex[[1]])))
 guessed.active.diff <- list()
 if(fusion != F){
  guessed.active.diff[[1]] <- matrix(nrow = nrow(R), ncol = P)
 }
 if(hasV){
  guessed.active[[2]] <- numeric(length(unique(grpindex[[2]])))
  guessed.active.coef[[2]] <- matrix(nrow = K, ncol = length(unique(grpindex[[2]])))
 }

 for(j in unique(grpindex[[1]])){
  j.which <- which(grpindex[[1]] == j)
  guessed.active[[1]][j] <- any(coef[[1]][, j.which] != 0)
  guessed.active.coef[[1]][,j] <- apply(coef[[1]][, j.which, drop=F], 1, function(u) any(u != 0))
  if(fusion != F && ( (!is.numeric(constr) & K > 1) || (is.numeric(constr) & K > 2) )){
   guessed.active.groupdiff[[1]][j] <- any(R%*%coef[[1]][, j.which] != 0)
   guessed.active.diff[[1]][,j.which] <- apply(R%*%coef[[1]][, j.which, drop=F], 1, function(u) any(u != 0))
  }
 }
 if(hasV){
  lp <- length(unique(grpindex[[1]]))
  for(j in unique(grpindex[[2]])){
   j.which <- which(grpindex[[2]] == j)
   guessed.active[[2]][j - lp] <- any(coef[[2]][, j.which] != 0)
   guessed.active.coef[[2]][,j - lp] <- apply(coef[[2]][, j.which, drop=F], 1, function(u) any(u != 0))
  }
 }   

 ## residual matrix or whatever. still to come
 res <- NULL
 
 AICfit <- -2 * logl + 2 * df
 if(all(round(weights) == weights)){tnobs <- wnobs}else{tnobs <- nobs}
 BICfit <- -2 * logl + log(tnobs) * df
 if(name %in% c("Sequential Logit Model")){
  brier <- Brier(y = cbind(dat$y, 1-rowSums(dat$y)), mu = cbind(mu, 1-rowSums(mu)), weights = weights)
 }else{
  brier <- Brier(y = dat$y, mu = mu, weights = weights)
 }
 
 ## fisher matrix of the model
 if(control@fisher){
  fish <- fisher(dat=dat, mu=mu, weights=weights, Proximal.args=Proximal.args)
 }else{
  fish <- NULL
 }
  
 ## names of the data:
 colnames(coef[[1]]) <- colnames(dat$x)
 rownames(coef[[1]]) <- colnames(dat$y)
 if(hasV){
  colnames(coef[[2]]) <- colnames(dat$V[[1]])
  rownames(coef[[2]]) <- colnames(dat$y)
 } 
 colnames(coef.stand[[1]]) <- colnames(dat$x)
 rownames(coef.stand[[1]]) <- colnames(dat$y)
 if(hasV){
  colnames(coef.stand[[2]]) <- colnames(dat$V[[1]])
  rownames(coef.stand[[2]]) <- colnames(dat$y)
 } 
 colnames(coef.pretres[[1]]) <- colnames(dat$x)
 rownames(coef.pretres[[1]]) <- colnames(dat$y)
 if(hasV){
  colnames(coef.pretres[[2]]) <- colnames(dat$V[[1]])
  rownames(coef.pretres[[2]]) <- colnames(dat$y)
 }
 
 ## dont attach the data to the object if keepdat is FALSE.
 y <- dat$y
 if(!keepdat){
  dat <- NULL
  y <- NULL
 }

 if(!keepdat | !standardize){
  x.stand <- NULL
  x.original <- NULL
 }
 if(!keepdat | !standardize | !hasV){ 
  V.stand <- NULL
  V.original <- NULL
 }

 ## set the class of the coef objects to MRSP.coef
 class(coef) <- "MRSP.coef"
 class(coef.stand) <- "MRSP.coef"
 class(coef.pretres) <- "MRSP.coef"

 ## list with all arguments present in the current environment. this is needed
 ## to compute the refits for a list of lambda values more efficiently.
 if(refit == "L" | control@keeparglist == T){
  arglist <- list()
  helplist <- ls()[-which(ls() == "arglist" | ls() == "ybar")]
  if(control@noarglistdat == T){
   helplist <- helplist[-which(helplist %in% c("dat", "eta", "mu", "V.original", "V.stand", "Vs", "Vw",
                                               "V.centered", "decomp.V", "scale.pen.V", "scale.notpen.V",
                                               "x.centered", "x.original", "x.stand", "xw", "decomp.x",
                                               "scale.pen.x", "scale.notpen.x", "y", "Proximal.args",
                                               "mycall", "fish", "fisher", "eweights", "seweights",
                                               "sweights", "updated"))]
  }
  arglist <- mget(helplist, envir = as.environment(-1))
  names(arglist) <- helplist
 }else{arglist <- NULL}

 ## the final output, an object of class "MRSP"
 out <- new("MRSP",
             coef = coef,
             coef.stand = coef.stand,
             coef.pretres = coef.pretres,
             dat = dat,
             x.original = x.original,
             x.stand = x.stand,
             V.original = V.original,
             V.stand = V.stand,
             y = y,
             weights = weights,
             penindex = penindex,
             grpindex = grpindex,
             penweights = penweights,
             guessed.active = guessed.active,
             guessed.active.coef = guessed.active.coef,
             guessed.active.groupdiff = guessed.active.groupdiff,
             guessed.active.diff = guessed.active.diff,
             df = df,
             tuning = tuning.old,
             lambda = lambda,
             lambdaR = lambdaR,
             lambdaF = lambdaF,
             fusion = fusion,
             gamma = gamma,
             psi = psi,
             eta = eta,
             mu = mu,
             indg = indg,
             indcs = indcs,
             offset = offset,
             residuals = res,
             loglik = logl,
             penalty = pen,
             fn.val = fn.val,
             model = model,
             constr = constr,
             control = control,
             iter.count = iter.count,
             best.iter = best.iter,
             ridgestabil = ridgestabil,
             name = name,
             AIC = AICfit,
             BIC = BICfit,
             Brier = brier,
             threshold = threshold,
             refit = refit,
             fisher = fish,
             arglist = arglist,
             call = mycall)
             
 if(refit == T) out <- refit(out)
 return(out)
})


## this applies MRSP.fit to a sequence of lambda-values:
setMethod("MRSP.fit",
          signature(lambda = "list"),
function(dat, coef.init, coef.stand.init, coef.pretres.init, offset = rep(0, nrow(dat$y)),
                     weights = rep(1, nrow(dat$y)), grpindex = NULL, penindex = NULL,
                     lambda, lambdaR = lambda, lambdaF = lambda, gamma = 1, psi = 1,
                     indg = NULL, indcs = NULL,
                     model = multinomlogit(), constr = NULL, control = MRSP.control(),
                     fista.control = NULL, Proximal.control = NULL, Proximal.args = NULL,
                     penweights = NULL, mlfit = NULL, adaptive = FALSE, threshold = FALSE,
                     refit = FALSE, fusion = FALSE, nonneg = FALSE, ...)
{
 mycall <- match.call()
 ## if any objects that are not formals of MRSP.fit are passed via ... , we need
 ## to explicitly assign the corresponding objects:
 dotlist <- list(...)
 for(i in seq_len(length(dotlist))) assign(names(dotlist)[i], dotlist[[i]])
 if(is.expression(model)) model <- eval(model)
 
 if(control@expandcall){
  excl <- c(1, which(names(mycall) %in% c("dat", "model")))
  if(is.null(mycall$model)){
   mycall$model <- expression(multinomlogit())
   excl <- c(excl, length(mycall))
  }else{
   mycall[[which(names(mycall) == "model")]] <- switch(model@name,
    "Multinomial Logit Model" = expression(multinomlogit()),
    "Sequential Logit Model" = expression(sequentiallogit()),
    "CUB Binomial Logit Model" = expression(CUBbinomiallogit()),
    "OLS Regression" = expression(OLSreg())
   )
  }
  mycall[-excl] <- mget(names(mycall[-excl]), envir = as.environment(-1))
 }

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
 if((length(which.intercept) == 1) && which.intercept != 1){
  dat$x[,c(1, which.intercept)] <- dat$x[,c(which.intercept, 1)]
  penindex[[1]][c(1,which.intercept)] <- penindex[[1]][c(which.intercept,1)]
  grpindex[[1]][c(1,which.intercept)] <- grpindex[[1]][c(which.intercept,1)]
  if(!all(coef[[1]] == 0))
  coef[[1]][,c(1,which.intercept)] <- coef[[1]][,c(which.intercept,1)]
  which.intercept <- 1
 }
 if(length(which.intercept) > 1)
   stop("more than one intercept column specified")

 if(missing(penindex) | is.null(penindex)){
  penindex <- list()
  penindex[[1]] <- rep(1, P)
  penindex[[1]][1] <- 10
  if(hasV){
   penindex[[2]][indg] <- rep(2, length(indg))
   penindex[[2]][indcs] <- rep(3, length(indcs))
  }
 }
 
 #K <- ncol(dat$y)
 ## if the first response category is chosen as reference, this can cause 
 ## problems on many occasions. therefore, we rearrange stuff in this case.
 ## really? probably not! commented out for now!
 #if(constr == 1){
 # dat$y[,c(1, K)] <- dat$y[, c(K, 1)]
 # colnames(dat$y)[c(1, K)] <- colnames(dat$y)[c(K, 1)]
 # dat$V[c(1, K)] <- dat$V[c(K, 1)]
 # if(!(missing(mlfit) | is.null(mlfit)) | !(missing(coef.init) | is.null(coef.init)) |
 #    !(missing(coef.stand.init) | is.null(coef.stand.init)) | 
 #    !(missing(coef.pretres.init) | is.null(coef.pretres.init)) |
 #    !(missing(penweights) | is.null(penweights))  ){
 #     stop("something went completely wrong in the passing on of arguments for the case of 'constr == 1'. please use a response category different from the first one as reference.")
 # }  
 #}
  
 if(!(missing(grpindex) | is.null(grpindex))){
  if(!is.list(grpindex)){grpindex <- list(grpindex)}
 }
 if(!(missing(penindex) | is.null(penindex))){
  if(!is.list(penindex)){penindex <- list(penindex)}
 }


 if(!is.list(lambda)) stop("the supplied lambda was not a list although it should have been")
 if(length(lambda) == 1){lambda <- lambda[[1]]}else{lambda <- unlist(lambda)}
 if(!is.numeric(lambda)) lambda <- as.numeric(lambda)

 if(is.list(lambdaR)){
  if(length(lambdaR) == 1 ){lambdaR <- lambdaR[[1]]}else{lambdaR <- unlist(lambdaR)}
 }
 if(is.numeric(lambdaR) & (length(lambdaR) == 1)){
  lambdaR <- rep(lambdaR, length(lambda))
 }
 if(length(lambdaR) != length(lambda)) stop("lambdaR not of appropriate length")
 if(!is.numeric(lambdaR)) lambdaR <- as.numeric(lambdaR)

 lambdaR <- sort(lambdaR, decreasing = T)
 lambda <- sort(lambda, decreasing = T)
 nrlambda <- length(lambda)

 if(is.list(lambdaF)){
  if(length(lambdaF) == 1 ){lambdaF <- lambdaF[[1]]}else{lambdaF <- unlist(lambdaF)}
 }
 if(is.numeric(lambdaF) & (length(lambdaF) == 1)){
  lambdaF <- rep(lambdaF, length(lambda))
 }
 if(length(lambdaF) != length(lambda)) stop("lambdaF not of appropriate length")
 if(!is.numeric(lambdaF)) lambdaF <- as.numeric(lambdaF)

 lambdaF <- sort(lambdaF, decreasing = T)

 if((missing(mlfit) | is.null(mlfit)) & control@doMLfit == T){
  mlfit <- getmlfit(dat = dat, coef.init = coef.init, coef.stand.init = coef.stand.init,
                              offset = offset, weights = weights, grpindex = grpindex,
                              penindex = penindex, lambda = lambda, lambdaR = lambdaR,
                              lambdaF = lambdaF, gamma = gamma,
                              psi = psi, indg = indg, indcs = indcs, model = model, constr = constr,
                              control = control, fista.control = fista.control,
                              Proximal.control = Proximal.control, Proximal.args = Proximal.args,
                              adaptive = adaptive, threshold = threshold, refit = refit,
                              nonneg = nonneg, ...)
  control@doMLfit <- F
 }

 if(missing(penweights) | is.null(penweights)){
  penweights <- getpenweights(dat = dat, coef.init = coef.init, coef.stand.init = coef.stand.init,
                              coef.pretres.init = coef.pretres.init, offset = offset, weights = weights,
                              grpindex = grpindex, penindex = penindex, lambda = lambda, lambdaR = lambdaR,
                              lambdaF = lambdaF, gamma = gamma, psi = psi, indg = indg, indcs = indcs,
                              model = model, constr = constr, control = control, fista.control = fista.control,
                              Proximal.control = Proximal.control, Proximal.args = Proximal.args,
                              adaptive = adaptive, threshold = threshold, refit = refit, mlfit = mlfit,
                              fusion = fusion, nonneg = nonneg, ...)
 }

 refitold <- refit
 if(refit == T) refit <- "L"
 
 fitlist <- list()

 fitlist[[1]] <- MRSP.fit(dat = dat, coef.init = coef.init, coef.stand.init = coef.stand.init,
                          coef.pretres.init = coef.pretres.init, offset = offset, weights = weights,
                          grpindex = grpindex, penindex = penindex, lambda = lambda[1], lambdaR = lambdaR[1],
                          lambdaF = lambdaF[1], gamma = gamma, psi = psi, indg = indg, indcs = indcs,
                          model = model, constr = constr, control = control, 
                          fista.control = fista.control, Proximal.control = Proximal.control,
                          Proximal.args = Proximal.args, penweights = penweights, mlfit = mlfit,
                          adaptive = adaptive, threshold = threshold, refit = refit, fusion = fusion,
                          nonneg = nonneg, ...)
                          
 if(nrlambda > 1){
  control@keepdat <- F
  control@noarglistdat <- T
  for(pos in seq(2, nrlambda)){
   if((fitlist[[pos-1]]@control)@ridgestabilrf == T) control@ridgestabilrf <- T  
   
   fitlist[[pos]] <- MRSP.fit(dat = dat, coef.init = fitlist[[pos-1]]@coef,
                              coef.stand.init = fitlist[[pos-1]]@coef.stand,
                              coef.pretres.init = fitlist[[pos-1]]@coef.pretres,
                              offset = offset, weights = weights, grpindex = grpindex,
                              penindex = penindex, lambda = lambda[pos], lambdaR = lambdaR[pos],
                              lambdaF = lambdaF[pos], gamma = gamma, psi = psi, indg = indg,
                              indcs = indcs, model = model,
                              constr = constr, control = control, fista.control = fista.control,
                              Proximal.control = Proximal.control, Proximal.args = Proximal.args, 
                              penweights = penweights, mlfit = mlfit, adaptive = adaptive, 
                              threshold = threshold, refit = refit, fusion = fusion,
                              nonneg = nonneg, ...)
                              
  }
 }
 
 if(nrlambda == 1){fit <- fitlist[[1]]}
 if(nrlambda > 1){
  fit <- fitlist
  #fit[[nrlambda + 1]] <- mycall
  class(fit) <- "MRSP.list"
  fit <- structure(fit, call = mycall, dat = dat)
 }
 
 if(refitold == T){
  argl1 <- fit[[1]]@arglist
  if(!control@keeparglist) fit[[1]]@arglist <- NULL
  argl1$coef.init <- argl1$coef
  argl1$coef.stand.init <- argl1$coef.stand
  argl1$coef.pretres.init <- argl1$coef.pretres
  
  fit[[1]] <- refit(fit[[1]], argl1)
  rm(argl1)
  
  if(length(fit) > 1){
   for(pos in seq(2, length(fit)-1)){
    argl <- fit[[pos]]@arglist
    if(!control@keeparglist) fit[[pos]]@arglist <- NULL
    argl$dat <- fit[[1]]@dat
    argl$coef <- argl$coef.init <- fit[[pos-1]]@coef
    argl$coef.stand <- argl$coef.stand.init <- fit[[pos-1]]@coef.stand
    argl$coef.pretres <- argl$coef.pretres.init <- fit[[pos-1]]@coef.pretres
  
    fit[[pos]] <- refit(fit[[pos]], argl)
    rm(argl)
   }
  }
 }
 
 if(!control@keeparglist){
  for(pos in seq_len(length(fit)-1)){
   fit[[pos]]@arglist <- NULL
  }
 }
 
 fit <- asS4(fit)   
 return(fit)
})


## standard method for missing lambda
setMethod("MRSP.fit",
          signature(lambda = "missing"),
function(dat, coef.init, coef.stand.init, coef.pretres.init, offset = rep(0, nrow(dat$y)),
                     weights = rep(1, nrow(dat$y)), grpindex = NULL, penindex = NULL,
                     lambda, lambdaR = NULL, lambdaF = 0, gamma = 1, psi = 1, indg = NULL, indcs = NULL,
                     model = multinomlogit(), constr = NULL, control = MRSP.control(),
                     fista.control = NULL, Proximal.control = NULL, Proximal.args = NULL,
                     penweights = NULL, mlfit = NULL, adaptive = FALSE, threshold = FALSE,
                     refit = FALSE, fusion = FALSE, nonneg = FALSE, ...)
{
 mycall <- match.call()
 ## if any objects that are not formals of MRSP.fit are passed via ... , we need
 ## to explicitly assign the corresponding objects:
 dotlist <- list(...)
 for(i in seq_len(length(dotlist))) assign(names(dotlist)[i], dotlist[[i]])
 if(is.expression(model)) model <- eval(model)
 
 if(control@expandcall){
  excl <- c(1, which(names(mycall) %in% c("dat", "model")))
  if(is.null(mycall$model)){
   mycall$model <- expression(multinomlogit())
   excl <- c(excl, length(mycall))
  }else{
   mycall[[which(names(mycall) == "model")]] <- switch(model@name,
    "Multinomial Logit Model" = expression(multinomlogit()),
    "Sequential Logit Model" = expression(sequentiallogit()),
    "CUB Binomial Logit Model" = expression(CUBbinomiallogit()),
    "OLS Regression" = expression(OLSreg())
   )
  }
  mycall[-excl] <- mget(names(mycall[-excl]), envir = as.environment(-1))
 }

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
 if((length(which.intercept) == 1) && which.intercept != 1){
  dat$x[,c(1, which.intercept)] <- dat$x[,c(which.intercept, 1)]
  penindex[[1]][c(1,which.intercept)] <- penindex[[1]][c(which.intercept,1)]
  grpindex[[1]][c(1,which.intercept)] <- grpindex[[1]][c(which.intercept,1)]
  if(!all(coef[[1]] == 0))
  coef[[1]][,c(1,which.intercept)] <- coef[[1]][,c(which.intercept,1)]
  which.intercept <- 1
 }
 if(length(which.intercept) > 1)
   stop("more than one intercept column specified")

 if(missing(penindex) | is.null(penindex)){
  penindex <- list()
  penindex[[1]] <- rep(1, P)
  penindex[[1]][1] <- 10
  if(hasV){
   penindex[[2]][indg] <- rep(2, length(indg))
   penindex[[2]][indcs] <- rep(3, length(indcs))
  }
 }
 
 #K <- ncol(dat$y)
 ## if the first response category is chosen as reference, this can cause 
 ## problems on many occasions. therefore, we rearrange stuff in this case.
 ## really? probably not! commented out for now!
 #if(constr == 1){
 # dat$y[,c(1, K)] <- dat$y[, c(K, 1)]
 # colnames(dat$y)[c(1, K)] <- colnames(dat$y)[c(K, 1)]
 # dat$V[c(1, K)] <- dat$V[c(K, 1)]
 # if(!(missing(mlfit) | is.null(mlfit)) | !(missing(coef.init) | is.null(coef.init)) |
 #   !(missing(coef.stand.init) | is.null(coef.stand.init)) | 
 #   !(missing(coef.pretres.init) | is.null(coef.pretres.init)) |
 #   !(missing(penweights) | is.null(penweights))  ){
 #    stop("something went completely wrong in the passing on of arguments for the case of 'constr == 1'. please use a response category different from the first one as reference.")
 #}  
 #} 

 if(!(missing(grpindex) | is.null(grpindex))){
  if(!is.list(grpindex)){grpindex <- list(grpindex)}
 }
 if(!(missing(penindex) | is.null(penindex))){
  if(!is.list(penindex)){penindex <- list(penindex)}
 }

 if((missing(mlfit) | is.null(mlfit)) & control@doMLfit == T){
  mlfit <- getmlfit(dat = dat, coef.init = coef.init, coef.stand.init = coef.stand.init,
                              offset = offset, weights = weights, grpindex = grpindex,
                              penindex = penindex, lambda = 123, lambdaR = lambdaR,
                              lambdaF = lambdaF, gamma = gamma,
                              psi = psi, indg = indg, indcs = indcs, model = model, constr = constr,
                              control = control, fista.control = fista.control,
                              Proximal.control = Proximal.control, Proximal.args = Proximal.args,
                              adaptive = adaptive, threshold = threshold, refit = refit,
                              nonneg = nonneg, ...)
  control@doMLfit <- F
 }

 if(!missing(lambda))
   stop("something went wrong with the passing on of lambda in the calls to MRSP.fit")
 
 if(missing(penweights) | is.null(penweights)){
  penweights <- getpenweights(dat = dat, coef.init = coef.init, coef.stand.init = coef.stand.init,
                              coef.pretres.init = coef.pretres.init, offset = offset, weights = weights,
                              grpindex = grpindex, penindex = penindex, lambda = list(123), lambdaR = lambdaR,
                              lambdaF = lambdaF, gamma = gamma, psi = psi, indg = indg, indcs = indcs,
                              model = model, constr = constr, control = control, fista.control = fista.control,
                              Proximal.control = Proximal.control, Proximal.args = Proximal.args,
                              adaptive = adaptive, threshold = threshold, refit = refit, mlfit = mlfit,
                              fusion = fusion, nonneg = nonneg, ...)
 }
  
 nrlambda <- control@nrlambda
 lambdabase <- control@lambdabase
 lambdamax <- control@lambdamax
 if(is.null(lambdamax)){
  lambdamax <- getlambdamax(dat = dat, offset = offset, weights = weights,
                            grpindex = grpindex, penindex = penindex, psi = psi,
                            gamma = gamma, model = model, constr = constr,
                            control = control, penweights = penweights,
                            indg = indg, indcs = indcs, adaptive = adaptive,
                            threshold = threshold, refit = refit, fusion = fusion,
                            lambdaF = lambdaF, nonneg = nonneg)
 }
 lambdamin <- min(control@lambdamin, lambdamax / 100)
 lambdalogvals <- seq(log(lambdamin / lambdabase), log(lambdamax / lambdabase), length.out = nrlambda)
 lambda <- lambdabase * exp(lambdalogvals)
 if(missing(lambdaR) | is.null(lambdaR)) lambdaR <- lambda ## attention, lambdaR may not be a list yet here!
 if(missing(lambdaF) | is.null(lambdaF)) lambdaF <- lambda
                                         
 fit <- MRSP.fit(dat = dat, coef.init = coef.init, coef.stand.init = coef.stand.init,
                          coef.pretres.init = coef.pretres.init, offset = offset, weights = weights,
                          grpindex = grpindex, penindex = penindex, lambda = list(lambda),
                          lambdaR = list(lambdaR), lambdaF = list(lambdaF), gamma = gamma, psi = psi,
                          indg = indg, indcs = indcs, model = model, constr = constr, control = control,
                          fista.control = fista.control, Proximal.control = Proximal.control,
                          Proximal.args = Proximal.args, penweights = penweights, adaptive = adaptive,
                          threshold = threshold, refit = refit, mlfit = mlfit, fusion = fusion,
                          nonneg = nonneg, ...)

 if(length(fit) > 1){
  fit <- structure(fit, call = mycall, dat = dat)
 }
 fit <- asS4(fit)
 return(fit)
})                               
 
