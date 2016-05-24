##' method function
##'
##' A function to generate a roxygen template for a method of a generic S3 function. Normally, this would be called from
##' the function generic, see ?generic 
##'
##' @param meth character, the name of the method
##' @param gen character the associated generic method
##' @param oname name of object
##' @return a roxygen template for the method.
##' @export

method <- function(meth,gen,oname="obj"){
    
    cat("##' ",gen,".",meth," function\n",sep="")
    cat("##'\n")
    cat("##' \n")
    cat("##'\n")
    cat("##' @method ",gen," ",meth,"\n",sep="")
    cat("##' @param ",oname," an ",meth," object\n",sep="")    
    cat("##' @param ... additional arguments\n",sep="")  
    cat("##' @return ...\n",sep="")
    cat("##' @export\n",sep="")
    cat("\n",sep="")
    cat(gen,".",meth," <- function(",oname,",...){\n",sep="")
    cat("\n",sep="")
    cat("}\n")

}
