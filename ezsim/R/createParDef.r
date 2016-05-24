#' createParDef creates a new parameterDef object from a list of scalar parameters and a list of other parameters.
#' parameterDef is a short hand of "parameter definition". It defines parameters used by the \code{\link{dgp}} which is the most important part of a simulation. For each simulation,There is a particular set of parameter. parameterDef allow us to define several parameters for different simulation at once. There are two types of parameter in parameterDef, scalar parameters and other parameters.
#' Scalar parameters must be a scalar. Any vectors or matrix is regarded as a sequence of scalar parameters. For example, n=seq(10,50,10), first simulation  takes n=10, second simulation takes n=20 and so on.
#' Other parameters can be anything and it is banker over the scalar parameters.
#' For example, we would like to know how would the sample size affect the variance of the sample mean of normally distributed variable. We can set n=seq(10,50,10), mean=1 and sd=2.  (see example)
#' @name createParDef
#' @aliases parameterDef createParDef
#' @title Create a parameterDef Object.
#' @param selection A list of scalar parameters
#' @param banker A list of other parameters
#' @return A parameterDef object
#' @author TszKin Julian Chan \email{ctszkin@@gmail.com}
#' @export
#' @seealso \code{\link{setBanker.parameterDef}},\code{\link{setSelection.parameterDef}},\code{\link{evalFunctionOnParameterDef}},\code{\link{generate.parameterDef}}
#' @examples       
#' par_def1<-createParDef(selection=list(mean=1,sd=2,n=seq(10,50,10)))
#' 
#' par_def2<-createParDef()
#' setSelection(par_def2,mean=1,sd=2,n=seq(10,50,10))
#' 
#' identical(par_def1,par_def2)
#' 
#' evalFunctionOnParameterDef(par_def1, function() rnorm(n,mean,sd) )  # 10 random number
#' evalFunctionOnParameterDef(par_def1, function() rnorm(n,mean,sd), index=3)  # 30 random number
#' 
#' generate(par_def1)
#' 
#' # More than 1 selection parameters 
#' par_def3<-createParDef(selection=list(sd=2,mean=1:3,n=seq(10,50,10)))
#' 
#' generate(par_def3)

createParDef <-
function(selection=list(),banker=list()){
    x<-list(selection=list(),banker=list())
    class(x)<-"parameterDef"

    if (length(banker)>0)
        x<-do.call(function(...) setBanker.parameterDef(x,...),banker)
    if (length(selection)>0)
        x<-do.call(function(...) setSelection.parameterDef(x,...),selection)
    x
}
