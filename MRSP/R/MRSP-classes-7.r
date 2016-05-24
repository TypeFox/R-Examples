################################################################################
#####    Setting up the main classes of MRSP.                              #####
################################################################################
#####    Author: Wolfgang Pößnecker                                        #####
#####    Last modified: 11.03.2014, 20:26                                  #####
################################################################################

## control class
setClass("MRSP.control",
         representation = representation(
           max.iter     = "numeric",
           rel.tol      = "numeric",
           fisher       = "logical",
           standardize  = "logical",
           keepdat      = "logical",
           nrlambda     = "numeric",
           lambdamin    = "numeric",
           lambdamax    = "ANY",
           lambdabase   = "numeric",
           ridgestabil  = "logical",
           ridgestabilrf= "logical",
           lambdastabil = "numeric",
           extrastabil  = "ANY",
           adapt.expo   = "numeric",
           adaptcv      = "logical",
           doMLfit      = "logical",
           computeDF    = "logical",
           backtransf   = "logical",
           trace        = "logical",
           cut.to.rel   = "numeric",
           expandcall   = "logical",
           keeparglist  = "logical",
           noarglistdat = "logical",
           numerictresh = "function",
           SPG.iter.max = "numeric",
           SPG.accel    = "numeric",
           SPG.eps      = "numeric"),

         prototype = list(
           max.iter     = 2000,
           rel.tol      = 1e-6,
           fisher       = FALSE,
           standardize  = TRUE,
           keepdat      = TRUE,
           nrlambda     = 50,
           lambdamin    = 0.01,
           lambdamax    = NULL,
           lambdabase   = exp(1),
           ridgestabil  = FALSE,
           ridgestabilrf= FALSE,
           lambdastabil = 0.02,
           extrastabil  = NULL,
           adapt.expo   = 1,
           adaptcv      = FALSE,
           doMLfit      = TRUE,
           computeDF    = TRUE,
           backtransf   = TRUE,
           trace        = FALSE,
           cut.to.rel   = 0,
           expandcall   = TRUE,
           keeparglist  = FALSE,
           noarglistdat = FALSE,
           numerictresh = function(x){1/x^0.4},
           SPG.iter.max = 200,
           SPG.accel    = seq(0,200)/(seq(0,200) + 3),
           SPG.eps      = 1e-4),

         validity = function(object){

           if(ceiling(object@max.iter) != floor(object@max.iter) |
              object@max.iter <= 0)
             return("max.iter must be a positive integer!")

           if(object@rel.tol <= 0)
             return("rel.tol has to be positive")

           return(TRUE)
         }
)

MRSP.control <- function(max.iter = 2000, rel.tol = 1e-6, fisher = FALSE,
                         standardize = TRUE, keepdat = TRUE, nrlambda = 50,
                         lambdamin = 0.01, lambdamax = NULL,
                         lambdabase = exp(1), ridgestabil = FALSE, ridgestabilrf = FALSE,
                         lambdastabil = 0.02, extrastabil = NULL, adapt.expo = 1, adaptcv = FALSE,
                         doMLfit = TRUE, computeDF = doMLfit, backtransf = standardize,
                         trace = FALSE, cut.to.rel = 0, expandcall = TRUE,
                         keeparglist = FALSE, noarglistdat = FALSE,
                         numerictresh = function(x){min(max(1/x^0.4, 0.01), 0.20913)}, SPG.iter.max = 200,
                         SPG.accel = seq(0,SPG.iter.max)/(seq(0,SPG.iter.max) + 3),
                         SPG.eps = 1e-4){

  ## Purpose: Options for the MRSP functions
  ## ---------------------------------------------------------------------------
  ## Arguments:
  ## max.iter:    maximal number of iterations of the FISTA algorithm.
  ## rel.tol:     convergence relative tolerance; the smaller the more precise.
  ## standardize: should the predictors be centered and standardized?
  ## backtransf:  should the coefficients be backtransformed to the original
  ##              scale of the predictors? defaults to 'standardize', i.e. yes
  ##              if standardization was used and no otherwise.
  ## keepdat:  should the data (x, V, etc.) be stored in the MRSP-object that is
  ##           currently being computed? mostly used to prevent storing the data
  ##           several times if a whole list of MRSP-objects is computed.
  ## nrlambda: if no sequence of lambda values is specified, such a sequence is
  ##           constructed automatically - how many such lambdas should be used?
  ## lambdabase: the lambda values that are chosen automatically if necessary
  ##           are equidistant on the log-scale. lambdabase gives the base of
  ##           this logarithm.
  ## ridgestabil: should a little ridgepenalty of magnitude lambdastabil be 
  ##           applied to all unpenalized parameter groups ? a remark:
  ##           lambdastabil only gives the baseline strenght of the ridge 
  ##           penalty which is then multiplied with sqrt(log(p)/n), a factor
  ##           which is known to be optimal for the tuning parameter in the
  ##           lasso context, see Bühlmann & van de Geer (2011, Book) Lemma 6.2.
  ##           this choice worked very well in our simulations, but note that it
  ##           still is completely heuristic.
  ## extrastabil: a list of the same size and format as penindex with logical
  ##           entries. columns of the model matrices where extrastabil is TRUE
  ##           are subject to an untuned ridge penalty of magnitude lambdastabil.
  ##           this can be used as a cheap trick to make models identifiable
  ##           that would not be without it. for example the intercept column
  ##           and the time-specific intercepts in competing risks models.
  ## trace:    should various warnings and diagnostic informations be printed?
  ## cut.to.rel: the resulting standardized (!) coefficients whose abs value is 
  ##           smaller than cut.to.rel are set to 0. set it to 0 for no cutting.
  ## expandcall: if true, the arguments in the calls to MRSP.fit are evaluated,
  ##           with the exception of dat. this allows to bootstrap a single MRSP
  ##           object or to manually refit it. set to false if you don't need
  ##           this functionality and want to save some memory.
  ## adaptcv:  if adaptive = T, aka the initial estimator used for the adaptive
  ##           weights is a penalized fit, then should this cv be computed in
  ##           parallel?
  ## numerictresh: a function that computes the small threshold for setting coefs
  ##           or coefficient groups to zero in order to avoid numerical inade-
  ##           quacies which might prevent coefs that the algorithm undoubtedly
  ##           tries to shrink to zero from remaining at small, non-zero values.
  ##           its argument is the sample size. it should be chosen so that the
  ##           numerical threshold goes to zero for n -> infinity.
  ## SPG.iter.max: if something inside MRSP has to be computed using the
  ##           Smoothing Proximal Gradient method of Chen et al (2012), this
  ##           is its maximal number of iterations.
  ## SPG.accel: sequence of acceleration factors for SPG. do not change this.
  ## SPG.eps: absolute (!) tolerance for SPG.
  ## ---------------------------------------------------------------------------

  RET <- new("MRSP.control",
             max.iter     = max.iter,
             rel.tol      = rel.tol,
             fisher       = fisher,
             standardize  = standardize,
             keepdat      = keepdat,
             nrlambda     = nrlambda,
             lambdamin    = lambdamin,
             lambdamax    = lambdamax,
             lambdabase   = lambdabase,
             ridgestabil  = ridgestabil,
             ridgestabilrf= ridgestabilrf,
             lambdastabil = lambdastabil,
             extrastabil  = extrastabil,
             adapt.expo   = adapt.expo,
             adaptcv      = adaptcv,
             doMLfit      = doMLfit,
             computeDF    = computeDF,
             backtransf   = backtransf,
             trace        = trace,
             cut.to.rel   = cut.to.rel,
             expandcall   = expandcall,
             keeparglist  = keeparglist,
             noarglistdat = noarglistdat,
             numerictresh = numerictresh,
             SPG.iter.max = SPG.iter.max,
             SPG.accel    = SPG.accel,
             SPG.eps      = SPG.eps)
  RET
}


## a class for MRSP model objects that implement the loglikelihood, score
## function and so on.
setClass("MRSP.model", representation=representation(
     invlink  = "function",
     link     = "function",
     logl     = "function",
     loglik   = "function",
     gradient = "function",
     fisher   = "function",
     constraint = "ANY",
     check    = "function",
     sancheck = "function",
     name     = "character",
     comment  = "character"
))


## now the main class, most slots should be self-explanatory.
setClass("MRSP",
         representation = representation(
         coef = "list",
         coef.stand = "list",
         coef.pretres = "list",
         dat = "ANY",
         x.original = "ANY",
         x.stand = "ANY",
         V.original = "ANY",
         V.stand = "ANY",
         y = "ANY",
         weights = "vector",
         penindex = "list",
         grpindex = "list",
         penweights = "list",
         guessed.active = "list",
         guessed.active.coef = "list",
         guessed.active.groupdiff = "list",
         guessed.active.diff = "list",
         df = "numeric",
         tuning = "ANY",
         lambda = "numeric",
         lambdaR = "numeric",
         lambdaF = "numeric",
         fusion = "ANY",
         gamma = "numeric",
         psi = "numeric",
         eta = "ANY",
         mu = "ANY",
         offset = "vector",
         residuals = "ANY",
         mlfit = "ANY",
         AIC = "numeric",
         BIC = "numeric",
         Brier = "numeric",
         threshold = "ANY",
         refit = "ANY",
         indg = "ANY",
         indcs = "ANY",
         model = "MRSP.model",
         constr = "ANY",
         control = "MRSP.control",
         fn.val = "numeric",
         loglik = "numeric",
         penalty = "numeric",
         iter.count = "numeric",
         best.iter = "numeric",
         ridgestabil = "ANY",
         name = "character",
         fisher = "ANY",
         arglist = "ANY",
         call = "ANY")
)

setClass("MRSP.list")#,
        # representation = representation("list")
#)


setClass("MRSP.sim",
         representation = representation("list")
)

## proxy class for the calls to the FISTA-generics:
setClass("MRSP.coef",
         representation = representation("list")
)

setClass("MRSP.bootstrap",
         representation = representation("list")
)

setGeneric("se",
           function(object, dat, B = 50, bootinds = NULL, method = "bootstrap", what = "coef",
                    parallel = TRUE, cores = detectCores(), ...)
           standardGeneric("se"))
           
setGeneric("pval",
           function(object, SE, what = "coef", ...)
           standardGeneric("pval"))
           
setClass("MRSP.se",
         representation = representation("list")
)

## generic "MRSP.fit"
setGeneric("MRSP.fit",
           function(dat, coef.init=NULL, coef.stand.init=NULL, coef.pretres.init=NULL, offset=NULL,
                    weights=NULL, grpindex=NULL, penindex=NULL, lambda, lambdaR=lambda, lambdaF=lambda,
                    gamma=1, psi=1, indg=NULL, indcs=NULL, model=NULL, constr=NULL, control=NULL,
                    fista.control=NULL, Proximal.control=NULL, Proximal.args=NULL, penweights=NULL,
                    mlfit=NULL, adaptive=FALSE, threshold=FALSE, refit=FALSE, fusion=FALSE, nonneg=FALSE, ...)
           standardGeneric("MRSP.fit"))
           
## generic "MRSP"
setGeneric("MRSP",
           function(formula, data, class.names=NULL, model=multinomlogit(), constr=NULL, offset=NULL, weights=NULL, penweights=NULL,
                    standardize=TRUE, nrlambda=50, lambdamin=0.01, lambdamax=NULL,  control=NULL, penalty=TRUE,
                    group.classes=TRUE, group.dummies=TRUE, sparse.groups=FALSE, adaptive=FALSE, threshold=FALSE, refit=FALSE,
                    lambda, lambdaR=0, lambdaF=0, gamma=1, psi=1, fusion=FALSE, nonneg=FALSE, y=NULL, X=NULL,
                    Z=NULL, penindex=NULL, grpindex=NULL, mlfit=NULL, perform.fit=TRUE, ...)                   # perform.fit = F means that the prepared call to MRSP.fit is returned, but not evaluated!
           standardGeneric("MRSP"))

## generic "cv"
setGeneric("cv",
           function(object, k = 10, type = "deviance", cvinds = NULL, parallel = TRUE, cores = detectCores(),
                    adaptcv = FALSE, initialseed = sample(1000, size=1), ...)
           standardGeneric("cv"))
      
## a class for MRSP crossvalidation objects
setClass("MRSP.cv",
         representation = representation("list")
)           

## generic "simulation"
setGeneric("simulation",
           function(coef, calls, nrep = 50, nobs = 100, nnewdata = 3*nobs, generator, type = "deviance",
                    k = 10, refit = FALSE, keepfit = FALSE, keepcv = TRUE, parallel = TRUE, cores = detectCores(),
                    initialseed = sample(1000, size=1), adaptcv = FALSE, cvreltol = 0.01, theta = 0.5, ...)
           standardGeneric("simulation"))
           
## generic "createResponse"
setGeneric("createResponse", 
           function(coef, dat, model = multinomlogit(), weights, offset, theta = 0.5, ...)
           standardGeneric("createResponse"))

## generic "MSE"
setGeneric("MSE",
          function(object, type = "coef", scaled = TRUE, ...)
          standardGeneric("MSE"))
          
## generic "FPR.FNR"
setGeneric("FPR.FNR",
           function(object, type = "mixed", ...)
           standardGeneric("FPR.FNR"))
           
## generic "refit"
setGeneric("refit",
           function(object, arglist, ...)
           standardGeneric("refit")) 
           
## generic "which.min.cv"
setGeneric("which.min.cv",
           function(x, reltol = 0.01, abstol = NULL, left = TRUE, ...)
           standardGeneric("which.min.cv"))

## generic "which.max.cv"
setGeneric("which.max.cv",
           function(x, reltol = 0.01, abstol = NULL, left = TRUE, ...)
           standardGeneric("which.max.cv"))        

## generic "bootstrap"
setGeneric("bootstrap",
           function(object, dat, what = NULL, B = 50, bootinds = NULL, parallel = TRUE, cores = detectCores(), ...)
           standardGeneric("bootstrap"))
           
## generic "predict"
setGeneric("predict",
           function(object, ...)
           standardGeneric("predict"))
           
## generic "select"
setGeneric("select",
           function(object, ...)
           standardGeneric("select"))
           
## make sure 'getCall' is S4-generic
setGeneric("getCall",
           function(x, ...)
           standardGeneric("getCall"))
           
