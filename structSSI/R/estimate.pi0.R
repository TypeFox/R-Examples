estimate.pi0 <- function(pvalues, method, alpha = 0.05, lambda = 0.5){
    method <- tolower(method)
    matched.method <- match.arg(method, c("tst", "lsl", "storey"))
    if(matched.method == "tst"){
        return(pi0.tst(pvalues, alpha))
    } else if(matched.method == "lsl"){
        return(pi0.lsl(pvalues))
    } else if(matched.method == "storey"){
        return(pi0.tail.p(lambda, pvalues))
    }
}
