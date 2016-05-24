full.FindIt <- function(object){
    ## Make the full factorial design matrix as the target population.
    comb.list <- apply(object$treat.orig,2,unique) 
    full <- as.data.frame(expand.grid(comb.list))
    invisible(full)
}
