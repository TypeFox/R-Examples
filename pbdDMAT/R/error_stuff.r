# error checking for the lazy man
must.be <- function(x, type)
{
    Rstuff <- c("character", "numeric", "integer", "double", "logical", "matrix", "data.frame", "vector")
    mystuff <- c("dmat", "ddmatrix", "ddvector")
    
    type <- pbdMPI::comm.match.arg(type, c(Rstuff, mystuff))
    
    nm <- deparse(substitute(x))
    
    fun <- eval(parse(text=paste("is.", type, sep="")))
    
    if (!fun(x))
        comm.stop(paste("argument '", nm, "' must be of type ", type, sep=""), call.=FALSE)
    
    invisible(TRUE)
}



