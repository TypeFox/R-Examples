summary.simexaft <-
function (object, ...) 
{
    p.names <- names(object$coefficients)
    est <- object$coefficients
    est.table <- list()
      
    se <- object$se
    pvalue <-object$pvalue
    est.table <- cbind(est, se, pvalue)
    dimnames(est.table) <- list(p.names, c("Estimate","Std. Error", "P value"))
       
    ans <- list()
    class(ans) <- "summary.simaxaft"
    ans$coefficients <- est.table
    ans$call <- object$call
    ans$scalereg <- object$scalereg
    ans$extrapolation <- object$extrapolation
    ans$SIMEXvariable <- object$SIMEXvariable
    ans

}
