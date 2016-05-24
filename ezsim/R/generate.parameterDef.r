#' Generate Parameters from a parameterDef Object.
#' The selection parameters in parameterDef is expanded and concatenated with the banker parameters.
#' @name generate.parameterDef
#' @aliases generate.parameterDef
#' @title Generate Parameter
#' @usage \method{generate}{parameterDef}(x,...)
#' @param x A parameterDef Object
#' @param \dots unused
#' @return \item{other_parameters}{A list of other_parameters} \item{scalar_parameters}{A data.frame of scalar_parameters} \item{parameter_list}{A list of all parameters}
#' @author TszKin Julian Chan \email{ctszkin@@gmail.com}
#' @export
#' @seealso \code{\link{setBanker.parameterDef}}, \code{\link{setSelection.parameterDef}}, \code{\link{evalFunctionOnParameterDef}}, \code{\link{generate.parameterDef}}
#' @examples         
#' par_def1<-createParDef(selection=list(mean=1,sd=2,n=seq(10,50,10)))
#' generate(par_def1)
#' par_def2<-createParDef(selection=list(sd=2,mean=1:3,n=seq(10,50,10)))
#' generate(par_def2)
generate.parameterDef <-
function(x,...){
    if (length(x$selection)>0){
        selection_parameters<-expand.grid(x$selection)
        apply(selection_parameters,1,function(xx) { 
            temp<-c(x$banker,xx)
            class(temp)<-c("parameters","list")
            temp
        })
    } else {
        temp<-x$banker
        class(temp)<-c("parameters","list")
        list(temp)
    }
}
