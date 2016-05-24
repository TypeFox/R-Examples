###################################################################################
##                               Strategy.R                                      ##
###################################################################################

###################################################################################
##' @include global.R
NULL
###################################################################################

###################################################################################
##' Create an instance of [\code{\linkS4class{Strategy}}] class 
##'
##' This class will contain all the parameters needed by the estimation algorithms.
##'
##' There are different ways to initialize an algorithm :
##'
##'    \describe{
##'
##'        \item{random}{Initialization from a random position is a standard way to
##'        initialize an algorithm. This random initial position is obtained by
##'        choosing at random centers in the data set. This simple strategy is
##'        repeated \eqn{5} times (the user can choose the number of times) from
##'        different random positions and the position that maximises the
##'        likelihood is selected.}
##'
##'        \item{smallEM}{A maximum of \eqn{50} iterations of the EM algorithm according to the process : \eqn{n_i} numbers of iterations
##'        of EM are done (with random initialization) until the \code{smallEM} stop criterion value has been reached. 
##'        This action is repeated until the sum of \eqn{n_i}
##'
##'        reaches \eqn{50} iterations (or if in one action \eqn{50} iterations are reached before the stop criterion value).\\
##'        It appears that repeating runs of EM is generally profitable since using a single run
##'        of EM can often lead to suboptimal solutions.}
##'
##'        \item{CEM}{\eqn{10} repetitions of \eqn{50} iterations of the CEM algorithm are done.
##'        One advantage of initializing an algorithm with CEM lies in the fact
##'        that CEM converges generally in a small number of iterations. Thus,
##'        without consuming a large amount of CPU times, several runs of CEM are
##'        performed. Then EM is run with the best solution among the \eqn{10} repetitions.}
##'
##'        \item{SEMMax}{A run of \eqn{500} iterations of SEM. The idea is that an SEM sequence is
##'        expected to enter rapidly in the neighbourhood of the global maximum
##'        of the likelihood function.}
##'
##'    }
##'
##' Defining the algorithms used in the strategy, the stopping rule and when to stop.
##'    \itemize{
##'        \item Algorithms :
##'           \describe{
##'               \item{EM}{Expectation Maximisation}
##'               \item{CEM}{Classification EM}
##'               \item{SEM}{Stochastic EM}
##'           }
##'        \item Stopping rules for the algorithm :
##'           \describe{
##'               \item{nbIterationInAlgo}{Sets the maximum number of iterations}
##'               \item{epsilonInAlgo}{Sets relative increase of the log-likelihood criterion}
##'           }
##'        \item Default values are \eqn{200} \code{nbIterationInAlgo} of \code{EM} with an \code{epsilonInAlgo} value of \eqn{10-3}.
##'    }
##' 
##' @param algo list of character string with the estimation algorithm.  Possible values: "EM", "SEM", "CEM", c("EM","SEM"). Default value is "EM".
##' @param nbTry integer defining the number of tries. nbTry must be a positive integer. Option available only if \code{init} is "random" or "smallEM" or "CEM" or "SEMMax". Default value: 1.
##' @param initMethod a character string with the method of initialization of the algorithm specified in the \code{algo} argument. Possible values: "random", "smallEM", "CEM", "SEMMax". Default value: "smallEM".
##' @param nbTryInInit integer defining number of tries in \code{initMethod} algorithm. nbTryInInit must be a positive integer. Option available only if \code{init} is "smallEM" or "CEM". Default value: 50.
##' @param nbIterationInInit integer defining the number of "EM" or "SEM" iterations in \code{initMethod}. nbIterationInInit must be a positive integer. Only available if \code{initMethod} is "smallEM" or "SEMMax". Default values: 5 if \code{initMethod} is "smallEM" and 100 if \code{initMethod} is "SEMMax".
##' @param nbIterationInAlgo list of integers defining the number of iterations if you want to use nbIteration as rule to stop the algorithm(s). Default value: 200. 
##' @param epsilonInInit real defining the epsilon value in the initialization step. Only available if \code{initMethod} is "smallEM". Default value: 0.001.
##' @param epsilonInAlgo list of reals defining the epsilon value for the algorithm. Warning: epsilonInAlgo doesn't have any sens if \code{algo} is SEM, so it needs to be set as NaN in that case. Default value: 0.001.
##' @param seed a positive integer defining the seed of the random number generator. Setting a particular seed allows the user to (re)-generate a particular serie of random numbers. NULL or negative value for a random seed.
##'
##' @examples
##'    mixmodStrategy()
##'    mixmodStrategy(algo="CEM",initMethod="random",nbTry=10,epsilonInInit=0.00001)
##'    mixmodStrategy(algo=c("SEM","EM"), nbIterationInAlgo=c(200,100), epsilonInAlgo=c(NA,0.000001))
##'
##' @references 
##'   R. Lebret, S. Iovleff, F. Langrognet, C. Biernacki, G. Celeux, G. Govaert (2015), "Rmixmod: The R Package of the Model-Based Unsupervised, Supervised, and Semi-Supervised Classification Mixmod Library", Journal of Statistical Software, 67(6), 1-29, doi:10.18637/jss.v067.i06
##'   Biernacki, C., Celeux, G., Govaert, G., 2003. "Choosing starting values for the EM algorithm for getting the highest likelihood in multivariate gaussian mixture models". Computational Statistics and Data Analysis 41, 561-575.
##'
##' @return a [\code{\linkS4class{Strategy}}] object
##' @author Remi Lebret and Serge Iovleff and Florent Langrognet, with contributions from C. Biernacki and G. Celeux and G. Govaert \email{contact@@mixmod.org}
##' @export
##'
mixmodStrategy <- function( algo="EM", nbTry=1, initMethod="smallEM", nbTryInInit=50, nbIterationInInit=5, nbIterationInAlgo=200, epsilonInInit=0.001, epsilonInAlgo=0.001, seed=NULL ){
  # create a new class Strategy
  new("Strategy", algo=algo, nbTry=nbTry, initMethod=initMethod, nbTryInInit=nbTryInInit, nbIterationInInit=nbIterationInInit, nbIterationInAlgo=nbIterationInAlgo, epsilonInInit=epsilonInInit, epsilonInAlgo=epsilonInAlgo, seed=seed)
}
###################################################################################


###################################################################################
##' Constructor of [\code{\linkS4class{Strategy}}] class
##'
##' This class defines the Mixmod strategies.
##'
##' \describe{
##'   \item{algo}{list of character string with the estimation algorithm.  Possible values: "EM", "SEM", "CEM", c("EM","SEM"). Default value is "EM".}
##'   \item{nbTry}{integer defining the number of tries. Default value: 1.}
##'   \item{initMethod}{a character string with the method of initialization of the algorithm specified in the \code{algo} argument. Possible values: "random", "smallEM", "CEM", "SEMMax". Default value: "smallEM".}
##'   \item{nbTryInInit}{integer defining number of tries in \code{initMethod} algorithm. Default value: 50.}
##'   \item{nbIterationInInit}{integer defining the number of "EM" or "SEM" iterations in \code{initMethod}. Default values: 5 if \code{initMethod} is "smallEM" and 100 if \code{initMethod} is "SEMMax".}
##'   \item{nbIterationInAlgo}{list of integers defining the number of iterations if user want to use nbIteration as rule to stop the algorithm(s). Default value: 200.} 
##'   \item{epsilonInInit}{real defining the epsilon value in the initialization step. Only available if \code{initMethod} is "smallEM". Default value: 0.001.}
##'   \item{epsilonInAlgo}{list of reals defining the epsilon value for the algorithm. Warning: epsilonInAlgo doesn't have any sens if \code{algo} is SEM, so it needs to be set as NaN in that case. Default value: 0.001.}
##'   \item{seed}{integer defining the seed of the random number generator. Setting a particular seed allows the user to (re)-generate a particular serie of random numbers. Default value is NULL, i.e. a random seed.}
##' }
##'
##' @examples
##'   new("Strategy")
##'   new("Strategy", algo="SEM", initMethod="SEMMax")
##'
##'   getSlots("Strategy")
##'
##' @name Strategy-class
##' @rdname Strategy-class
##' @exportClass Strategy
##'
setClass(
    Class="Strategy",
    representation=representation(
        algo = "character",
        nbTry = "numeric",
        initMethod = "character",
        nbTryInInit = "numeric",
        nbIterationInInit = "numeric",
        nbIterationInAlgo = "numeric",
        epsilonInInit = "numeric",
        epsilonInAlgo = "numeric",
        seed = "numeric"
    ),
    prototype=prototype(
        algo = "EM",
        nbTry = 1,
        initMethod = "smallEM",
        nbTryInInit = 50,
        nbIterationInInit = 5,
        nbIterationInAlgo = 200,
        epsilonInInit = 0.001,
        epsilonInAlgo = 0.001,
        seed = -1
    ),
    # validity function
    validity=function(object){
      # for algo
      if ( sum(object@algo %in% c("EM","SEM","CEM")) != length(object@algo) ){
        stop("At least one algorithm is not valid. See ?mixmodAlgo for the list of available algorithms.")
      }
      # for 'initMethod'
      if( (object@initMethod != "smallEM") & (object@initMethod != "random") & (object@initMethod != "CEM") & (object@initMethod != "SEMMax") ){
        stop("initMethod name is not valid.")
      }
      # for 'nbTry'
      if (!is.wholenumber(object@nbTry)){
        stop("nbTry must be an integer.")
      }
      if (object@nbTry < 1){
        stop("nbTry must be positive.")
      }
      # for 'nbTryInInit'
      if ( (object@initMethod == "smallEM") | (object@initMethod == "CEM") ){
        if (!is.wholenumber(object@nbTryInInit)){
          stop("nbTryInInit must be an integer.")
        }
        if (object@nbTryInInit < 1){
          stop("nbTryInInit must be positive.")
        }
      }
      # for 'epsilonInInit'
      if ( (object@initMethod == "smallEM") ){
        if (!is.double(object@epsilonInInit)){
          stop("epsilonInInit must be a real. Default value will be used.")
        }
        if ( (object@epsilonInInit > 1) | (object@epsilonInInit < 0) ){
          stop("epsilonInInit must be less than one and positive.")
        }
      } 
      # for 'nbIterationInInit'
      if ( (object@initMethod == "smallEM") | (object@initMethod == "SEMMax") ){
        if (!is.wholenumber(object@nbIterationInInit)){
          stop("nbIterationInInit must be an integer.")
        }
        if(object@nbIterationInInit < 1){
          stop("nbIterationInInit must be positive.")
        }
      }
      # for 'nbIterationInAlgo'
      if (sum(!is.wholenumber(object@nbIterationInAlgo))){
        stop("nbIterationInAlgo must be an integer.")
      }
      if ( min(object@nbIterationInAlgo,na.rm=TRUE) < 1){
        stop("nbIterationInAlgo must be positive.")
      }
      
      # for 'epsilonInAlgo'
      if ( sum(!is.nan(object@epsilonInAlgo[which(object@algo=="SEM")])) ){
        stop("epsilonInAlgo must be NaN for the SEM algorithm.")
      }
      if (!is.double(object@epsilonInAlgo[which(object@algo!="SEM")])){
        stop("epsilonInAlgo must be a real.")
      }
      if ( sum(!is.nan(object@epsilonInAlgo)) ){
        if( (max(object@epsilonInAlgo,na.rm=TRUE) > 1) | (min(object@epsilonInAlgo,na.rm=TRUE) < 0) ) {
          stop("epsilonInAlgo must be less than one and positive.")
        }
      }
      return(TRUE)
    }
)
###################################################################################


###################################################################################
##' Create an instance of the [\code{\linkS4class{Strategy}}] class using new/initialize.
##' 
##' Initialization method. Used internally in the `Rmixmod' package.
##' 
##' @seealso \code{\link{initialize}}
##'
##' @keywords internal
##'
##' @rdname initialize-methods
##'
setMethod(
  f="initialize",
  signature=c("Strategy"),
  definition=function(.Object,algo,nbTry,initMethod,nbTryInInit,nbIterationInInit,nbIterationInAlgo,epsilonInInit,epsilonInAlgo,seed
){
    if(!missing(algo)){
      if(length(algo)>1){
        # for epsilon in Algo
        if( missing(epsilonInAlgo) ){ 
          .Object@epsilonInAlgo<-rep(0.001,length(algo))
        }
        else if(length(epsilonInAlgo)!=length(algo)){
          .Object@epsilonInAlgo<-epsilonInAlgo[1:length(algo)]
          .Object@epsilonInAlgo[is.na(.Object@epsilonInAlgo)]<-0.001
        }
        else{
          .Object@epsilonInAlgo<-epsilonInAlgo
        }
        # for nbIteration in Algo
        if( missing(nbIterationInAlgo) ){ 
          .Object@nbIterationInAlgo<-rep(200,length(algo))
        }
        else if(length(nbIterationInAlgo)!=length(algo)){
          .Object@nbIterationInAlgo<-nbIterationInAlgo[1:length(algo)]
          .Object@nbIterationInAlgo[is.na(.Object@nbIterationInAlgo)]<-200
        }
        else{
          .Object@nbIterationInAlgo<-nbIterationInAlgo
        }
        # check whether SEM algo is in list to set epsilon as NaN
        if(sum(algo=="SEM")){
          .Object@epsilonInAlgo[which(algo=="SEM")]<-NaN
        }
      }
      else{
        if(algo=="SEM"){ .Object@epsilonInAlgo<-NaN }
        else if( missing(epsilonInAlgo) ){ .Object@epsilonInAlgo<-0.001 }
        else{.Object@epsilonInAlgo<-epsilonInAlgo[1]}
        if(missing(nbIterationInAlgo)){ .Object@nbIterationInAlgo<-200 }
        else{.Object@nbIterationInAlgo<-nbIterationInAlgo[1]}
      }
      .Object@algo<-algo
    }
    else{
      .Object@algo<-"EM"
      if(missing(epsilonInAlgo)){ .Object@epsilonInAlgo<-0.001 }
      else{.Object@epsilonInAlgo<-epsilonInAlgo[1]}
      if(missing(nbIterationInAlgo)){ .Object@nbIterationInAlgo<-200 }
      else{.Object@nbIterationInAlgo<-nbIterationInAlgo[1]}
    }
    
    if(!missing(seed)){ 
      if(is.null(seed)){
        .Object@seed<-(-1)
      }else{
        .Object@seed<-seed 
      }
    }
    else{.Object@seed<-(-1)}

    if(!missing(nbTry)){ .Object@nbTry<-nbTry }
    else{.Object@nbTry<-1}

    if(!missing(nbTryInInit)){  .Object@nbTryInInit<-nbTryInInit }
    else{.Object@nbTryInInit<-50}

    if(!missing(nbIterationInInit)){  .Object@nbIterationInInit<-nbIterationInInit }
    else{.Object@nbIterationInInit<-5}

    if(!missing(epsilonInInit)){  .Object@epsilonInInit<-epsilonInInit }
    else{.Object@epsilonInInit<-0.001}
  
    if(!missing(initMethod)){ 
      .Object@initMethod<-initMethod 
      # change default value if necessary
      if ( (.Object@initMethod == "SEMMax") & (.Object@nbIterationInInit == 5) ){
        .Object@nbIterationInInit <- 100
      }
    }
    else{.Object@initMethod<-"smallEM"}
  
    validObject(.Object)        
    return(.Object)
  }
)
###################################################################################


###################################################################################
##' @rdname print-methods
##' @aliases print print,Strategy-method
##'
setMethod(
  f="print",
  signature=c("Strategy"),
  function(x,...){
    cat("****************************************\n")
    cat("*** MIXMOD Strategy:\n")
    cat("* algorithm            = ", x@algo, "\n")
    cat("* number of tries      = ", x@nbTry, "\n")
    cat("* number of iterations = ", x@nbIterationInAlgo, "\n")
    cat("* epsilon              = ", x@epsilonInAlgo, "\n")
    cat("*** Initialization strategy:\n")
    cat("* algorithm            = ", x@initMethod, "\n")
    cat("* number of tries      = ", x@nbTryInInit, "\n")
    cat("* number of iterations = ", x@nbIterationInInit, "\n")
    cat("* epsilon              = ", x@epsilonInInit, "\n")
    cat("* seed                 = ", ifelse(x@seed<0,"NULL",x@seed), "\n")
    cat("****************************************\n")
  }
)
###################################################################################


###################################################################################
##' @rdname show-methods
##' @aliases show show,Strategy-method
##'
setMethod(
  f="show",
  signature=c("Strategy"),
  function(object){
    cat("****************************************\n")
    cat("*** MIXMOD Strategy:\n")
    cat("* algorithm            = ", object@algo, "\n")
    cat("* number of tries      = ", object@nbTry, "\n")
    cat("* number of iterations = ", object@nbIterationInAlgo, "\n")
    cat("* epsilon              = ", object@epsilonInAlgo, "\n")
    cat("*** Initialization strategy:\n")
    cat("* algorithm            = ", object@initMethod, "\n")
    cat("* number of tries      = ", object@nbTryInInit, "\n")
    cat("* number of iterations = ", object@nbIterationInInit, "\n")
    cat("* epsilon              = ", object@epsilonInInit, "\n")
    cat("* seed                 = ", ifelse(object@seed<0,"NULL",object@seed), "\n")
    cat("****************************************\n")
  }
)
###################################################################################


###################################################################################
##' @rdname extract-methods
##' @aliases [,Strategy-method
##'
setMethod(
  f="[", 
  signature(x = "Strategy"),
  definition=function(x,i,j,drop){
    if ( missing(j) ){
      switch(EXPR=i,
        "algo"={return(x@algo)},
        "nbTry"={return(x@nbTry)},
        "nbIterationInAlgo"={return(x@nbIterationInAlgo)},
        "epsilonInAlgo"={return(x@epsilonInAlgo)},
        "initMethod"={return(x@initMethod)},
        "nbTryInInit"={return(x@nbTryInInit)},
        "nbIterationInInit"={return(x@nbIterationInInit)},
        "epsilonInInit"={return(x@epsilonInInit)},
        "seed"={return(x@seed)},
        stop("This attribute doesn't exist !")
      )
    }else{
      stop("This attribute is not a list !")
    }
  }
)
###################################################################################



###################################################################################
##' @name [
##' @rdname extract-methods
##' @aliases [<-,Strategy-method
##'
setReplaceMethod(
  f="[", 
  signature(x = "Strategy"), 
  definition=function(x,i,j,value){
    if ( missing(j) ){
      switch(EXPR=i,
        "algo"={x@algo<-value},
        "nbTry"={x@nbTry<-value},
        "nbIterationInAlgo"={x@nbIterationInAlgo<-value},
        "epsilonInAlgo"={x@epsilonInAlgo<-value},
        "initMethod"={x@initMethod<-value},
        "nbTryInInit"={x@nbTryInInit<-value},
        "nbIterationInInit"={x@nbIterationInInit<-value},
        "epsilonInInit"={x@epsilonInInit<-value},
        "seed"={x@seed<-value},
        stop("This attribute doesn't exist !")
      )
    }else{
      stop("This attribute is not a list !")
    }
    validObject(x)
    return(x)
  }
)
###################################################################################

