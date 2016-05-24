# wrapper for stepAIC in the MASS package

# last modified 2014-08-04 by J. Fox

stepwise <- function(mod, 
    direction=c("backward/forward", "forward/backward", "backward", "forward"), 
    criterion=c("BIC", "AIC"), ...){
    criterion <- match.arg(criterion)
    direction <- match.arg(direction)
    cat("\nDirection: ", direction)
    cat("\nCriterion: ", criterion, "\n\n")
    k <- if (criterion == "BIC") log(nrow(model.matrix(mod))) else 2
    rhs <- paste(c("~", deparse(formula(mod)[[3]])), collapse="")
    rhs <- gsub(" ", "", rhs)
    if (direction == "forward" || direction == "forward/backward")
        mod <- update(mod, . ~ 1)
    if (direction == "backward/forward" || direction == "forward/backward") direction <- "both"
    lower <- ~ 1
    upper <- eval(parse(text=rhs))   
    stepAIC(mod, scope=list(lower=lower, upper=upper), direction=direction, k=k, ...)
}