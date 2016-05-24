print.aliases <- function(x, ...){
   ### function for printing output objects from function aliases in nicer format
   if (!"aliases" %in% class(x)) stop("applicable for objectso of class aliases only")
   if ("aliases" %in% names(x)){
     als <- which(sapply(x$aliases,"length")>1)
     if (length(als)>0){
        alprint <- list(legend=x$legend, aliases=x$aliases[als])
        alprint$aliases <- t(t(sapply(alprint$aliases, "paste", collapse=" = ")))
        colnames(alprint$aliases) <- ""
        rownames(alprint$aliases) <- rep("",nrow(alprint$aliases))
        if(!is.null(alprint$legend)) names(alprint$legend)<-rep("", length(alprint$legend))
        }
     else alprint <- list(legend=NULL, aliases="no aliasing in the model")
     if (is.null(alprint$legend)) print(alprint$aliases,quote=FALSE, ...) else
     print(alprint,quote=FALSE, ...)
   }
   else print.default(x, quote=FALSE, ...)
}