summary.rrs <- function(object, ...){
    if(!object$log) cat("Relative risk function.\n\n")
    else cat("Log-Relative risk function.\n\n")
    
    
    if(!object$log) cat("Estimated risk range ",min(as.vector(object$rsM),na.rm=T)," to ",max(as.vector(object$rsM),na.rm=T),".\n",sep="")
    else cat("Estimated log-risk range ",min(as.vector(object$rsM),na.rm=T)," to ",max(as.vector(object$rsM),na.rm=T),".\n",sep="")
    cat(sum(!is.na(as.vector(object$rsM))),"grid cells out of",prod(dim(object$rsM)),"fall inside study region.\n")
    
    #cat("Surface (Z) summary:\n")
#    print(summary(as.vector(object$rsM)))
#    
    cat("\n--Numerator (case) density--\n")
    summary.bivden(object$f)
    cat("\n--Denominator (control) density--\n")
    summary.bivden(object$g)
    
}
    
