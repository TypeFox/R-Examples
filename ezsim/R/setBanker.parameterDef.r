#' setBanker sets the scalar parameters of a parameterDef object. setBanker are "call by reference", so assignment is not needed to update the parameterDef object. In other words, they will overwrite the value of its argument(parameterDef object).
#' parameterDef is a short hand of "parameter definition". It defines parameters used by the \code{\link{dgp}} which is the most important part of a simulation. For each simulation,There is a particular set of parameter. parameterDef allow us to define several parameters for different simulation at once. There are two types of parameter in parameterDef, scalar parameters and other parameters.
#' Scalar parameters must be a scalar. Any vectors or matrix is regarded as a sequence of scalar parameters. For example, n=seq(10,50,10), first simulation  takes n=10, second simulation takes n=20 and so on.
#' Other parameters can be anything and it is banker over the scalar parameters.
#' For example, we would like to know how would the sample size affect the variance of the sample mean of normally distributed variable. We can set n=seq(10,50,10), mean=1 and sd=2.  (see example)
#' @name setBanker.parameterDef
#' @aliases setBanker.parameterDef
#' @title Set a parameterDef Object.
#' @usage \method{setBanker}{parameterDef}(x,...)
#' @param x A parameterDef object
#' @param \dots Variables to be added to a parameterDef object
#' @return A parameterDef object
#' @author TszKin Julian Chan \email{ctszkin@@gmail.com}
#' @export
#' @seealso \code{\link{setSelection.parameterDef}},\code{\link{createParDef}},\code{\link{evalFunctionOnParameterDef}},\code{\link{generate.parameterDef}}
#' @keywords parameterDef
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
#' 
#' # 
#' par_def4<-createParDef(selection=list(mean=1,sd=2,n=seq(10,50,10)))
#' setBanker(par_def4,some_matrix=matrix(1:4,nrow=2),some_vector=1:6)
#' par_def4
#' generate(par_def4)

setBanker.parameterDef <-
function(x,...){
    temp<-list(...)
    i<-NULL
    if (length(temp)==0)
        stop("Nothing to set")
    
    # remove parameter in x$selection if parameter appear in ...
    if ( length(x$selection)!=0   && length(intersect(names(temp),names(x$selection)))>0){
        warning("Same variables name as selection parameter. selection paramter will be removed")
        x$selection <- x$selection[[names(x$selection) %in% names(temp)]] 
    }

    # concat those parameter which is not repeated
    repeated<-names(temp) %in% names(x$banker)
    x$banker<-c(x$banker,temp[!repeated])
    
    # for those who are repeated, replace the original value
    if (sum(repeated)>0){
        warning(paste("\nOriginal value of ",names(temp)[repeated] ," is replaced"))
        foreach (i = names(temp)[repeated]) %do% {
            x$banker[[i]]<-temp[[i]]
        }        
    }

		x
}
