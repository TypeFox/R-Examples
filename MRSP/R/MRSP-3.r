################################################################################
#####    An R Package for fitting (M)ultinomial (L)ogit models with a      #####
#####    (S)tructured (P)enalization for global and category-specific      #####
#####    predictors.                                                       #####
################################################################################
#####    Author: Wolfgang Pößnecker                                        #####
#####    Last modified: 04.11.2014, 17:21                                  #####
################################################################################

## function 'MRSP' is actually just a convenience function that facilitates the
## fitting of the models by calls to the actual workhose, MRSP.fit
## -------------------------------------------------------------------------- ##
## Mandatory Arguments: 
   ## y: the response. either a vector that consists of integers from 1 to k, 
   ##    where k = max(y), or a character string giving the name of the response
   ##    variable to be found in 'data'. 
   ## x: a RHS formula or the design matrix for the mixture components/binomial
   ##    models.
   ## data: a data frame containing the variables specified in y, x and possibly
   ##    z. factor variables with more than one dummy will automatically get
   ##    penalized jointly if you leave the construction of argument grpindex to
   ##    MRSP
   ## lambda: the tuning parameter for the lasso penalty that is applied to 
   ##    regression coefficients within the mixture components. it must be a 
   ##    single numeric number.
## Optional Arguments: 
   ## z: a RHS formula or the design matrix for the concomitant model, i.e. the
   ##    multinomial logit model that tries to explain the probabilities/mixing
   ##    proportions of the overall model with the variables found in z. if left
   ##    unspecified, constant mixing probabilities are used.
   ## c.lambda: lambda for the lasso penalty on the regression parameters of the
   ##    concomitant model. 
   ## penindex: a vector of length ncol(x) that specifies how each predictor in
   ##    x shall be penalized. the number 10 means the corresponding predictor 
   ##    is not penalized, 12 stands for a lasso penalty, 13 for ridge. do not
   ##    use 1, 11 or 14 with MRSP, the corresponding penalties wouldnt make
   ##    sense in the CUB context. if left unspecified, all true covariates are
   ##    subject to a lasso penalty and the intercept is left unpenalized.
   ## grpindex: a vector of length ncol(x) that contains integer numbers in 
   ##    ascending order, starting from 1. columns of x that share the same 
   ##    number in grpindex are treated as a group whose parameters are pena-
   ##    lized jointly via a group lasso penalty. it is best to leave this argu-
   ##    ment unspecified. for factor variables, the corresponding dummies 
   ##    should always be treated as such a group. this is achieved by supplying
   ##    the factor variable as a column of type 'factor' within x. MRSP
   ##    then appropriately constructs dummy variables and adjust penindex and
   ##    grpindex.
   ## weights: a vector of observation weights.
   ## offset: an offset vector.
   ## lambdaR: the lambda for ridge penalties. use this is you want to use a mix
   ##    of lasso and ridge penalties and want them to use different tuning 
   ##    parameters. if ridge is the only penalty you want to use, you can just
   ##    use the lambda argument instead.
   ## adaptive: should adaptive weights be used for the lasso penalties on the
   ##    component model parameters? set 'adaptive = "ML" ' if you want to use 
   ##    it. for both conceptual and computational reasons, I discourage it at 
   ##    the moment. never ever set 'adaptive=T', you wouldnt like the outcome.
   ## refit: if true, refitting is used. like adaptive weights, this can improve
   ##    the variable selection performed by the lasso, but it is very time con-
   ##    suming in the mixture context that we have here.
   ## penweights: an object that can be used to supply different, user-specified
   ##    weights for the penalty of different predictors by hand. ask Wolfgang
   ##    Pößnecker for details if you want to use this.
   ## lambdaF, fusion, mlfit: included for compatibility reasons, but you should 
   ##    never use them.
   ## concomitant: name of a genmix driver function for the concomitant model.
   ##    do not change this. 
   ## c.penindex, c.grpindex, c.offset and so on: same arguments as above, but 
   ##    this time referring to the concomitant model and thus the z variables.
   ## model.coef.init, concomitant.coef.init: you can start MRSP at specific
   ##    values of the model.coef and concomitant.coef using these arguments.
   ## postweight.init, pi.init: two matrices of dimension length(y) x M that 
   ##    contain initial values of the prior and posterior weights of the model.
   ## standardize: if true, all predictor variables are centered and standard-
   ##    ized by MRSP before the penalized fit is computed internally via
   ##    MRSP.fit. it is strongly recommended to use this standardization.
   ## control: the control object for the mixture model/the EM algorithm facili-
   ##    tated by genmix. 
   ## model.control, concomitant.control: control objects of class MRSP.control
   ##    for the component and concomitant models, respectively. DO NOT CHANGE!
   ## checks.default: a check function for genmix. you might want to include 
   ##    some additional checks into it.
   ## initialseed: used as a little trick to make results independent from 
   ##    random initializations across different versions of MRSP. see the
   ##    comments above 'genmix' for details.
   ## psi: a factor in [0,2] that balances the penalization on global and cate-
   ##    gory-specific predictors. global ones are weighted with psi, catspec.
   ##    ones are weighted with (2-psi). the default of psi=1 thus means that
   ##    both predictor types are weighted equally.
## -------------------------------------------------------------------------- ##
setMethod("MRSP",
          character(),      #signature(lambda = "ANY"),
cmpfun(function(formula, data, class.names=NULL, model=multinomlogit(), constr=NULL, offset=NULL, weights=NULL,
                    penweights=NULL, standardize=TRUE, nrlambda=50, lambdamin=0.01, lambdamax=NULL,
                    control=MRSP.control(nrlambda=nrlambda, lambdamin=lambdamin, lambdamax=lambdamax),
                    penalty=TRUE, group.classes=TRUE, group.dummies=TRUE, sparse.groups=FALSE, adaptive=FALSE,
                    threshold=FALSE, refit=FALSE, lambda, lambdaR=0, lambdaF=0, gamma=1, psi=1,
                    fusion=FALSE, nonneg=FALSE, y=NULL, X=NULL, Z=NULL, penindex=NULL, grpindex=NULL,
                    mlfit=NULL, perform.fit=TRUE, ...)
{
 mycall <- match.call()
 ## if any objects that are not formals of MRSP or MRSP.fit are passed via ... , we need
 ## to explicitly assign the corresponding objects:
 dotlist <- list(...)
 for(i in seq_len(length(dotlist))) assign(names(dotlist)[i], dotlist[[i]])
 if(is.expression(model)) model <- eval(model)

 ## replace all entries of mycall with their actual value. this is useful, for
 ## example, for bootstrapping or predict methods (and so on...)
 if(control@expandcall){
  excl <- c(1, which(names(mycall) %in% c("data", "model")))
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
 if(model@name %in% c("Sequential Logit Model", "Cumulative Logit Model")){
  if(constr != "none"){
   constr <- "none"
   warning(paste("Identifiability constraints are not needed in ordinal regression. Argument 'constr' was changed to value 'none'!"))
  }
 }
 ## in multinomial models, a symmetric side constraint means that the parameters
 ## of a covariate for the different response categories must sum to zero. Thus,
 ## it cannot be combined with a nonnegativity constraint on the parameters.
 if(nonneg == T && constr == "symmetric"){
  constr <- 1
  warning()
  cat(" You cannot combine a nonnegativity constraint on the coefficients with a symmetric side constraint!", "\n", "The first response category is instead used as reference.","\n")
 }

 ## check that we have enough info about the model specification to be able to work with it:
 if(missing(formula) && any(c(missing(y), missing(X), missing(grpindex), missing(penindex))))
   stop("if you do not specify a formula in MRSP, you must correctly supply arguments y, X, grpindex, penindex and (possibly) Z!")

 if(!missing(formula)){
  ## make sure the specificed formula has the required 'Formula' class:
  if(!(is(formula, "formula") | is(formula, "Formula"))) formula <- as.formula(formula)
  formula <- Formula(formula)

  ## in the following, create, from the data,  the response matrix y, a design
  ## matrix called X that contains global/individual-specific covariates and,
  ## possibly, a list called Z that contains category-/alternative-specific
  ## covariates. due to case distinctions, it becomes a bit messy....
  dataobjects <- modelmatrix.MRSP(formula, data, model@name)      ## creates all the data objects that we need
 
  formula <- dataobjects$formula                   # an object of class 'Formula', i.e. with rhs and lhs attributes
  if(missing(y) | is.null(y)) y <- dataobjects$y   # response matrix, it is a 0-1 matrix of dimension nobs x K, with rowSums(y) = rep(1,nobs)
  if(missing(X) | is.null(X)) X <- dataobjects$X   # overall X-matrix
  hasZ <- !(missing(Z) | is.null(Z)) | !is.null(dataobjects$Z)  ## indicator whether there are any category-specific predictors
  if(hasZ){
   if(missing(Z) | is.null(Z)) Z <- dataobjects$Z  # data object with category-specific predictors, stored in a length K list of nobs x (L1+L2) matrices
  }
  len.form <- dataobjects$len.form                 # length of the rhs formula
  #which.intercept <- dataobjects$which.intercept  # index of the intercept column
  #has.intercept <- dataobjects$has.intercept      # indicator for the presence of an intercept
  if(!is.null(dataobjects$x)){
   x <- dataobjects$x                              # same as X for multinomial or cumulative models, the global predictors with global effects in sequential models
   x.vars <- dataobjects$x.vars                    # name of all the actual variables used in x. disregards, for example, multiple dummies for one and the same categorical predictor
  }
  if(!is.null(dataobjects$x2)){
   x2 <- dataobjects$x2                            # in sequential models, these are the global predictors with category-specific effects
   x2.vars <- dataobjects$x2.vars
  }
  if(!is.null(dataobjects$z)){
   z <- dataobjects$z                              # category-specific predictors with global effects, stored in matrix form
   z.vars <- dataobjects$z.vars
  }
  if(!is.null(dataobjects$z2)){
   z2 <- dataobjects$z2                            # category-specific predictors with category-specific effects. stored in matrix form
   z2.vars <- dataobjects$z2.vars
  }
 
  ## construct the names of all actual variables used in X:
  X.vars <- NULL                                                                 # implicitly, this assumed that there always is at least one global, i.e. 'X-type' variable
  if(exists("x.vars",inherits=F)) X.vars <- c(X.vars, x.vars)
  if(exists("x2.vars",inherits=F)) X.vars <- c(X.vars, x2.vars)
  ## now the same for Z:
  if(hasZ){
   Z.vars <- NULL
   if(exists("z.vars",inherits=F)) Z.vars <- c(Z.vars, z.vars)
   if(exists("z2.vars",inherits=F)) Z.vars <- c(Z.vars, z2.vars)
  }
  
  ## check that formulas do not contain 'I(...)' since MRSP cant't handle it:
  #evalh(Iformulacheck)
  if(exists("x",inherits=F)){
   for(j in seq_along(x.vars)){
    if(any(sapply(strsplit(colnames(x), x.vars[j]), function(u) u[1] == "I(")))
      stop("MRSP currently does not support the use of variable transformations ",
           "via 'I(...)'-constructs in the x-formula")
   }
  }
  if(exists("x2",inherits=F)){
   for(j in seq_along(x2.vars)){
    if(any(sapply(strsplit(colnames(x2), x2.vars[j]), function(u) u[1] == "I(")))
      stop("MRSP currently does not support the use of variable transformations ",
           "via 'I(...)'-constructs in the x-formula")
   }
  }
  if(exists("z",inherits=F)){
   for(j in seq_along(z.vars)){
    if(any(sapply(strsplit(colnames(z), z.vars[j]), function(u) u[1] == "I(")))
      stop("MRSP currently does not support the use of variable transformations ",
           "via 'I(...)'-constructs in the z-formula")
   }
  }
  if(exists("z2",inherits=F)){
   for(j in seq_along(z2.vars)){
    if(any(sapply(strsplit(colnames(z2), z2.vars[j]), function(u) u[1] == "I(")))
      stop("MRSP currently does not support the use of variable transformations ",
           "via 'I(...)'-constructs in the z-formula")
   }
  }
 }

 ## number of classes/repsonse categories/alternatives:
 K <- ncol(y)
 nobs <- nrow(y)
 if(missing(class.names) | is.null(class.names)){
  class.names <- paste0("Class ", seq(1, K))
  if(model@name %in% c("Sequential Logit Model")) class.names[K+1] <- paste0("Class ",K+1)
 }
 if(length(class.names) == K){
  colnames(y) <- class.names
 }else if(length(class.names) == K+1){                                          ## this is the case for ordinal regression since we there always have K-1 coefs per predictor. (unlike multinom models with, e.g. symmetric constraint)
  #colnames(y) <- class.names[-length(class.names)]
  #class.names[K+1] <- paste0("Class ", K+1)
  colnames(y) <- paste(class.names[1:K],class.names[2:(K+1)], sep="|")
 }
 
 ## get the weights for the observations (not to be confused with 'penweights'):
 if(missing(weights) | is.null(weights))  weights <- rep(1, nobs)
 ## get the offset for observations...
 if(missing(offset) | is.null(offset)) offset <- rep(0, nobs)
 
 #### now comes the section in which we prepare the penindex and grpindex for
 #### for the global predictors, i.e. the 'X-part'. first some intercept handling:
 which.intercept <- which(apply(X, 2, function(u) diff(range(u)) < .Machine$double.eps ^ 0.5))
 if(length(which.intercept) > 1) stop("too many intercept columns supplied; or multicollinearity is present in the dataset")
 has.intercept <- length(which.intercept == 1)
 if(has.intercept && which.intercept != 1) X[,c(1,which.intercept)] <- X[,c(which.intercept,1)]
 if(has.intercept) colnames(X)[1] <- "(Intercept)"
 #if(!has.intercept){
 # x <- cbind(1, x)
 # colnames(x)[1] <- "(Intercept)"
 # has.intercept <- T
 #}


 ####                         penindex section                                  #####
 ## create the penindex object, which specifies how different predictors are
 ## penalized by MRSP.
 if(missing(penindex) | is.null(penindex)){
  penindex <- list()
  penindex[[1]] <- numeric(ncol(X))
  if(model@name == "Multinomial Logit Model"){
   xind <- charmatch(colnames(x), colnames(X))          # index of global variables to be equipped with category-specific coefs
   #penindex <- evalh(getpenindex.x.catspec)
   if(has.intercept){
    xind <- xind[-1]
    penindex[[1]][1] <- 10          # penindex = 10 means unpenalized
   }
   if(penalty==F){
    penindex[[1]][xind] <- 10
   }else if(penalty==T){
    if(group.classes){
     if(sparse.groups){
      penindex[[1]][xind] <- 11     # penindex = 11 means sparse CATS lasso
     }else{
      penindex[[1]][xind] <- 1      # penindex = 1 means CATS lasso without within-group sparsity
     }
    }else{
     penindex[[1]][xind] <- 12      # penindex = 12 means that an ordinary lasso or (depending on group.dummies)
                                    # groups-across-dummies--grouplasso is used, but no grouping over repsonse categories
    }
   }else if(penalty=="Ridge" | penalty=="ridge"){
    penindex[[1]][xind] <- 13       # penindex = 13 means a ridge penalty
   }
  }else if(model@name %in% c("Cumulative Logit Model", "Sequential Logit Model")){
   if(has.intercept)   penindex[[1]][[1]] <- 10
   xglobind <- charmatch(setdiff(colnames(x),"(Intercept)"), colnames(X))         # index of global variables to be equipped with global coefs
   #penindex <- evalh(getpenindex.x.global)
   if(penalty==F){
    penindex[[1]][xglobind] <- 40       # penindex = 40 means global coefficient, unpenalized
   }else if(penalty==T){
    penindex[[1]][xglobind] <- 4        # penindex = 4 means global coefficient, penalized by a lasso or grouplasso-type penalty (depending on group.dummies)
   }else if(penalty=="Ridge" | penalty=="ridge"){
    penindex[[1]][xglobind] <- 41       # penindex = 41 means global coefficient with a ridge penalty
   }
   if(model@name %in% c("Sequential Logit Model") && exists("x2",inherits=F)){
    xind <- charmatch(setdiff(colnames(x2),"(Intercept)"), colnames(X))
    #penindex <- evalh(getpenindex.x.catspec)
    if(penalty==F){
     penindex[[1]][xind] <- 10
    }else if(penalty==T){
     if(group.classes){
      if(sparse.groups){
       penindex[[1]][xind] <- 11     # penindex = 11 means sparse CATS lasso
      }else{
       penindex[[1]][xind] <- 1      # penindex = 1 means CATS lasso without within-group sparsity
      }
     }else{
      penindex[[1]][xind] <- 12      # penindex = 12 means that an ordinary lasso or (depending on group.dummies)
                                     # groups-across-dummies--grouplasso is used, but no grouping over repsonse categories
     }
    }
   }
  }
  if(hasZ){
   penindex[[2]] <- numeric(ncol(Z[[1]]))
   zglobind <- charmatch(colnames(z), colnames(Z[[1]]))
   #penindex <- evalh(getpenindex.z.global)
   if(penalty==F){
    penindex[[2]][zglobind] <- 20       # penindex = 20 means category-specific coefficient, unpenalized
   }else if(penalty==T){
    penindex[[2]][zglobind] <- 2        # penindex = 2 means category-specific coefficient, penalized by a lasso or grouplasso-type penalty (depending on group.dummies)
   }else if(penalty=="Ridge" | penalty=="ridge"){
    penindex[[2]][zglobind] <- 21       # penindex = 21 means category-specific coefficient with a ridge penalty
   }
   if(exists("z2",inherits=F)){
    zind <- charmatch(colnames(z2), colnames(Z[[1]]))
    #penindex <- evalh(getpenindex.z.catspec)
    if(penalty==F){
     penindex[[2]][zind] <- 30       # penindex = 30 means category-specific coef, unpenalized
    }else if(penalty==T){
     if(group.classes){
      if(sparse.groups){
       penindex[[2]][zind] <- 31     # penindex = 31 means sparse CATS lasso
      }else{
       penindex[[2]][zind] <- 3      # penindex = 3 means CATS lasso without within-group sparsity
      }
     }else{
      penindex[[2]][zind] <- 32      # penindex = 32 means that an ordinary lasso or (depending on group.dummies)
                                     # groups-across-dummies--grouplasso is used, but no grouping over repsonse categories
     }
    }else if(penalty=="Ridge" | penalty=="ridge"){
     penindex[[2]][zind] <- 33       # penindex = 33 means a ridge penalty
    }
   }
  }
 }
 ## check if the user-supplied penindex has the correct length
 if(!(length(penindex[[1]]) == ncol(X)))
   stop("the length of the specified penindex does not match the column number of the X-model matrix")
 if(hasZ && !(length(penindex[[2]]) == ncol(Z[[1]])))                                                   # this line does not require penindex[[2]] or Z to exist if hasZ = F!
   stop("the length of the specified penindex does not match the column number of the Z-model matrix")
 
 
 ####                         grpindex section                                  #####
 ## now, create the grpindex object that specifies which predictors are penalized
 ## as a joint group. grpindex is a vector whose length is the same as the number
 ## of columns in the design matrix. it contains integer values in ascending order.
 ## each number denotes one parameter group that is penalized jointly. columns that
 ## share the same number form one group. an example:
 ## grpindex = c(1,2,3,3,4,4,4) means that variables 1 and 2 form their own group.
 ## variables 3 and 4  as well as variables 5, 6 and 7 form groups.
 ## the following code creates grpindex and adjusts x and penindex for factor or
 ## ordered predictors.
 if(missing(grpindex) | is.null(grpindex)){
  grpindex <- list()
  grpindex[[1]] <- numeric()
  if(has.intercept | ncol(X) == 1) grpindex[[1]][1] <- 1       # the intercept should always form its own group

  if(ncol(X) > 1){
   ## create a matrix that indicates which actual variables are involved in the
   ## creation of which columns of the design matrix X:
   varind <- matrix(F, nrow=length(X.vars), ncol=ncol(X))
   for(i in seq(nrow(varind))){
    for(j in seq(ncol(varind))){
     ivarname <- X.vars[i]
     ivar <- data[,which(colnames(data) == ivarname)]
     if(is.factor(ivar)){
      ilevels <- levels(ivar)[-1]
      for(l in seq(ilevels)){
       ilevelname <- paste0(ivarname, ilevels[l])
       if(ilevelname == colnames(X)[j] |
          length(grep(paste0(ilevelname, ":"), colnames(X)[j])) > 0 |
          length(grep(paste0(":", ilevelname), colnames(X)[j])) > 0)
       { varind[i,j] <- T }
      }
     }else{
      if(ivarname == colnames(X)[j] |
         length(grep(paste0(ivarname, ":"), colnames(X)[j])) > 0 |
         length(grep(paste0(":", ivarname), colnames(X)[j])) > 0)
      { varind[i,j] <- T }
     }
    }
   }
   
   for(j in seq(2,ncol(X))){
    if(!all(varind[,j] == varind[,j-1])){
     grpindex[[1]][j] <- grpindex[[1]][j-1] + 1
    }else{
     grpindex[[1]][j] <- grpindex[[1]][j-1]
    }
   }
  }
 }

 ## now prepare the grpindex for Z-variables.
 if(hasZ){
  if(length(grpindex) < 2){
   grpindex[[2]] <- numeric()
   grpindex[[2]][1] <- max(grpindex[[1]]) + 1
   
   if(ncol(Z[[1]]) > 1){
    ## create a matrix that indicates which actual variables are involved in the
    ## creation of which columns of the matrices in Z:
    varind <- matrix(F, nrow=length(Z.vars), ncol=ncol(Z[[1]]))
    for(i in seq(nrow(varind))){
     for(j in seq(ncol(varind))){
      ivarname <- Z.vars[i]
      ivar <- data[,which(colnames(data) == ivarname)]
      if(is.factor(ivar)){
       ilevels <- levels(ivar)[-1]
       for(l in seq(ilevels)){
        ilevelname <- paste0(ivarname, ilevels[l])
        if(ilevelname == colnames(Z[[1]])[j] |
           length(grep(paste0(ilevelname, ":"), colnames(Z[[1]])[j])) > 0 |
           length(grep(paste0(":", ilevelname), colnames(Z[[1]])[j])) > 0)
        { varind[i,j] <- T }
       }
      }else{
       if(ivarname == colnames(Z[[1]])[j] |
          length(grep(paste0(ivarname, ":"), colnames(Z[[1]])[j])) > 0 |
          length(grep(paste0(":", ivarname), colnames(Z[[1]])[j])) > 0)
       { varind[i,j] <- T }
      }
     }
    }

    for(j in seq(2,ncol(Z[[1]]))){
     if(!all(varind[,j] == varind[,j-1])){
      grpindex[[2]][j] <- grpindex[[2]][j-1] + 1
     }else{
      grpindex[[2]][j] <- grpindex[[2]][j-1]
     }
    }
   }
  }
 }
 ####  end grpindex section  ####

 ## for technical reasons, penindex, penweights and grpindex must be lists in MRSP.fit
 if(!is.list(penindex)) penindex <- list(penindex)
 if(!is.list(grpindex)) grpindex <- list(grpindex)
 if(!is.null(penweights) && !is.list(penweights)) penweights <- list(penweights)
 
 ######################################################################################
 ##### code for the standardization of predictors. do not change the following    #####
 ##### code section unless you REALLY know what you're doing!                     #####
 ######################################################################################
 pen.which.x <- which(!(penindex[[1]] %in% c(10, 40)))
 notpen.which.x <- which(penindex[[1]] %in% c(10, 40))
 grpindex.pen.x <- grpindex[[1]][pen.which.x]
 dict.pen.x <- sort(unique(grpindex.pen.x))
 wnobs <- sum(weights)      # typically, wnobs = nobs

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
 if(standardize){
  x.original <- X
  sweights <- weights / sum(weights)
  x.centered <- X
  mean.x <- apply(X[,,drop=F], 2, function(u) weighted.mean(u, weights))
  if(has.intercept){
   x.centered[,-1] <- sweep(X[,-1, drop=F], 2, mean.x[-1])
  }else{x.centered <- sweep(X[,, drop=F], 2, mean.x)}

  xw <- x.centered
  scale.pen.x <- list()
  scale.notpen.x <- NULL
  ## the following lines are necessary to handle the backtransformation of
  ## standardized coefficients during the ml fit via MRSP in the final step of
  ## this script.
  scale.notpen.ML <- NULL
  if(has.intercept){
   scale.notpen.ML <- sqrt(drop(sweights %*% (x.centered[, -1]^2)))
  }else{
   scale.notpen.ML <- sqrt(drop(sweights %*% (x.centered^2)))
  }
  
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

  X <- xw

  if(hasZ){
   Z.original <- Z
   pen.which.Z <- which(!(penindex[[2]] %in% c(20, 30)))
   notpen.which.Z <- which(penindex[[2]] %in% c(20, 30))
   grpindex.pen.Z <- grpindex[[2]][pen.which.Z]
   dict.pen.Z <- sort(unique(grpindex.pen.Z))
   
   eweights <- rep(weights, times=K)
   seweights <- rep(sweights, times=K)
   ## extract Z into a matrix of size (nobs * K) x L:
   Zs <- do.call(rbind, Z)
   mean.Z <- apply(Zs, 2, function(u) weighted.mean(u, eweights))
   Z.centered <- sweep(Zs, 2, mean.Z)

   Zw <- Z.centered
   scale.pen.Z <- list()
   scale.notpen.Z <- NULL
   ## again some stuff that is necessary for backtransformations inside of getmlfit:
   scale.notpen.ZML <- sqrt(drop(seweights %*% (Z.centered^2)))
   
   if(length(notpen.which.Z) > 0){
    scale.notpen.Z <- sqrt(drop(seweights %*% (Z.centered[, notpen.which.Z]^2)))
    Zw[, notpen.which.Z] <- scale(Z.centered[, notpen.which.Z], FALSE, scale.notpen.Z)
   }

   for(j in dict.pen.Z){
    j.which <- which(grpindex[[2]] == j)
    decomp.Z <- qr(sqrt(eweights) * Z.centered[, j.which])
    if(decomp.Z$rank < length(j.which))
      stop("Block belonging to columns ", paste(j.which), "in Z does not have full rank!")
    scale.pen.Z[[j]] <- qr.R(decomp.Z) * 1/sqrt(wnobs * K)
    Zw[, j.which] <- qr.Q(decomp.Z) * sqrt(wnobs * K)
   }

   for(i in seq_len(K)){Z[[i]] <- Zw[seq(nobs*i-nobs+1,nobs*i),, drop=F]}
  }
 }
 ##### end of standardization #######################################################

 
 ## create a list of arguments to pass on to MRSP.fit
 arglist <- list(adaptive = adaptive,
                 refit = refit,
                 mlfit = mlfit,
                 penweights = penweights,
                 threshold = threshold,
                 penindex = penindex,
                 grpindex = grpindex,
                 fusion = fusion,
                 #lambda = lambda,                                              # this is handled separately later on due to case distinctions...
                 lambdaR = lambdaR,
                 lambdaF = lambdaF,
                 control = control,
                 model = model,
                 psi = psi,
                 constr = constr,
                 offset = offset,
                 weights = weights,
                 gamma = gamma,
                 nonneg = nonneg)
                 
 if(standardize){                                                               # those arguments are need for the backtransformation of standardized coefs to the original scale
  arglist$notpen.which.x <- notpen.which.x
  arglist$scale.notpen.x <- scale.notpen.x
  arglist$scale.notpen.ML <- scale.notpen.ML
  arglist$dict.pen.x <- dict.pen.x
  arglist$scale.pen.x <- scale.pen.x
  arglist$mean.x <- mean.x
  #arglist$control@standardize <- F
  #arglist$control@backtransf <- T
  slot(arglist$control, "standardize") <- F
  slot(arglist$control, "backtransf") <- T
  if(hasZ){                                                                     ## the following arguments are denoted with the letter 'V' in MRSP.fit. I should change that someday...
   arglist$notpen.which.V <- notpen.which.Z
   arglist$scale.notpen.V <- scale.notpen.Z
   arglist$scale.notpen.VML <- scale.notpen.ZML
   arglist$dict.pen.V <- dict.pen.Z
   arglist$scale.pen.V <- scale.pen.Z
   arglist$mean.V <- mean.Z
  }
 }

 ## do not evaluate the call entries during the ordinary fitting to save memory
 #arglist$control@expandcall <- F

 ## the data object for the call to MRSP.fit:
 dat <- list(y=y,x=X)
 if(hasZ) dat$V <- Z          ## it is a relict that the category-specific covariates are stored in an object called 'V' in MRSP.fit. too tedious to change that to 'Z'...

 ## prepare the call to MRSP.fit:
 fitcall  <- as.call(c(as.symbol("MRSP.fit"), list(dat=quote(dat)), arglist, dotlist))
 
 if(!missing(lambda)){
  if(is.numeric(lambda)){
   if(length(lambda) > 1){
    lambda <- sort(lambda, decreasing=T)
    lambda <- list(lambda)
   }
  }else stop("lambda must either be left out or supplied as a numeric vector or scalar.")
  fitcall$lambda <- lambda
 }else if(penalty == F){                                                        # this makes sure we fit only one model if no penalization if requested,
  fitcall$lambda <- list(12345567890)                                           # even if lambda is missing. (which would usually result in a grid of lambda values...)
 }


 #### now its finally time to compute stuff:
 if(perform.fit){
  fit <- eval(fitcall)                                                          ########### this is where the action happens!
  if(class(fit) == "MRSP.list"){
   fit <- lapply(fit, function(u){structure(u, topcall = mycall)})              ## assign the topcall that created the objects to every slot
   class(fit) <- "MRSP.list"
  }
  out <- structure(fit, topcall = mycall, call = fitcall, dat = dat)            ## topcall = call of the toplevel function, in this case MRSP. fitcall = call that is actually used for model fitting. dat = data object in the form required by MRSP.
  out <- asS4(out)                                                              ## this is necessary to signal generic functions that this object is not just a list...
  #if(exists("formula",inherits=F)){                  ## <- not needed since the formula is part of 'topcall' if it exists.
  # out <- structure(out, formula = formula)
  # if(class(out) == "MRSP.list"){
  #  out <- lapply(out, function(u){structure(u, formula = formula)})
  # }
  #}
 }else{
  out <- list(topcall = mycall, call = fitcall, dat=dat)                                            ## we dont need to store the topcall if we dont actually compute stuff...
  if(exists("formula",inherits=F)) out$formula <- formula
 }

 return(out)
}))






 
 
 
 
 
 
 
 
