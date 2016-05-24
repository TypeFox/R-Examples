utils::globalVariables("biocLite")
#' Install or update Bioconductor and CRAN packages
#' 
#' This is function simply a wrapper for \code{biocLite}. It first sources the code 
#' from the Bioconductor site then calls \code{biocLite}.
#' 
#' @param ... arguments passed on to \code{biocLite}
#' 
#' @author Rafael A. Irizarry
#' 
#' @details Note that once you run this function in a session, you no
#' longer need to call since 
#' you can call \code{biocLite} directly.
#' 
#' 
#' @examples
#' install_bioc("affy")
install_bioc <- function(...){
    internet <- try(source("http://bioconductor.org/biocLite.R"), silent=TRUE) 
    if(!class(internet)=="try-error") 
        { 
            biocLite(...)
        } else 
            { 
                stop("connection to http://bioconductor.org not successfull\n") 
            } 
}

