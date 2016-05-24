#' \code{Predict} ANFIS' network output
#'
#' Forward Pass to predict the ANFIS' output
#'
#' @param object ANFIS class object.
#' @param x numeric matrix [patterns x inputs] of input patterns.
#'
#' @return matrix with the output values
#'
#' @include Anfis-plot.R
#' @exportMethod predict
#' @docType methods
#' @name predict
#' @rdname ANFIS-predict
#' @aliases predict,ANFIS-method
#' @seealso \code{\link{ANFIS-class}}
#' @note see full example in \code{\link{ANFIS-class}}
#' @family ANFIS
#' @author Cristobal Fresno \email{cfresno@@bdmg.com.ar}, Andrea S. Llera 
#'  \email{ALlera@@leloir.org.ar} and Elmer A. Fernandez 
#'  \email{efernandez@@bdmg.com.ar}
setMethod(f="predict", signature="ANFIS", definition=function(object, x){
  ##Check for valid x
  stopifnot(is.matrix(x) & length(x)!=0 & ncol(x)==length(object@premises))
  ##Initialization of variables
  w <- NULL ## matrix R x nrow(X) with tita_2 output in each column
  wSum <- NULL ## vector with sum values
  wNormalized <- NULL ## For output generation
    
  A <- NULL ## for LSE consequents estimation
      
  ##Pattern Foward Half Pass  	 
  invisible(lapply(1:nrow(x),function(pattern){
    ##Obtain tita_1 output
    tita_1 <- sapply(seq(along=object@premises), function(input,pattern){
      unlist(lapply(object@premises[[input]], function(MF){
        evaluateMF(MF,x[pattern,input])}))
    }, pattern)
    ##If not homogenious MFs
    if(class(tita_1)=="list"){
      aux <- matrix(NA,nrow=max(sapply(tita_1,length)), ncol=ncol(x))
      invisible(sapply(seq(along=tita_1), function(input){
        aux[1:length(tita_1[[input]]),input] <<- tita_1[[input]]
        return(NULL)
      }))
      tita_1 <- aux
    }
    ##MF=1 for each input
    if(!is.matrix(tita_1)){
      tita_1 <- t(tita_1)
    } 
#   if(ncol(x)==1){tita_1 <- t(tita_1)}  ##Input=1

    ##Obtain tita_2 output a.k.a. w
    w<<-cbind(w,apply(object@rules,1,function(rule,tita_1){
      prod(sapply(seq(along=rule), function(input,tita_1){
        tita_1[rule[input],input]},tita_1))
    },tita_1))
    
    ## Obtain tita_3 a.k.a. w normalization 
    wSum<<-c(wSum, sum(w[, pattern]))
    wNormalized<<-cbind(wNormalized, w[,pattern]/wSum[pattern])

    ## Obtain tita_4 
    ## Generete the LSE system for the consequents
    # AX=B with each of the patterns 
    # (w1 x)p11 + (w1 y)q11 + (w1)r11 + ..... + 
    # + (w16 x)p1_16 + (w16 y)q1_16 + (w16)r1_16 = J1 for output 1
    # A is (Pattern)x(Rules*(input+1)); 
    # X is (Rules*(input+1))x(outputs) and B is (Patterns)x(outputs)
    A<<-rbind(A, matrix(matrix(c(x[pattern, ], 1), ncol=1) %*% 
      wNormalized[, pattern], nrow=1, byrow=FALSE))

    return(NULL)
  }))

  ##Rest of the Forward Pass
  ## Obtain tita_5 output estimation
  return(A%*%object@consequents)
})
