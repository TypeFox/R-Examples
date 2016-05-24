##' generic function
##'
##' A function to generate roxygen templates for generic funtions and associated methods.
##'
##' @param gen character string giving the name of an S3 generic.
##' @param methods character vector: a list of methods for which to provide templates
##' @param sp the amont of space to put in between functions
##' @param oname name of the generic object
##' @return roxygen text printed to the console.
##' @export

generic <- function(gen,methods=NULL,sp=3,oname="obj"){
    
    cat("##' ",gen," function\n",sep="")
    cat("##'\n")
    cat("##' \n")
    cat("##'\n")
    cat("##' @param ",oname," an object\n",sep="")    
    cat("##' @param ... additional arguments\n",sep="")  
    cat("##' @return method ",gen,"\n",sep="")
    cat("##' @export\n",sep="")
    cat("\n")
    cat(gen," <- function(",oname,",...){\n",sep="")
    cat("    UseMethod(\"",gen,"\")\n",sep="")
    cat("}\n")
    ct <- 0; while(ct<sp){cat("\n");ct<-ct+1}
    if(!is.null(methods)){
        for(m in methods){
            method(meth=m,gen=gen,oname=oname)
            ct <- 0; while(ct<sp){cat("\n");ct<-ct+1}
        }
    }
}
