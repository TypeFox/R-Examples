coxme <- function(formula,  data, 
        weights, subset, na.action, init, 
        control, ties= c("efron", "breslow"),
        varlist, vfixed, vinit, sparse=c(50,.02),
        x=FALSE, y=TRUE, 
        refine.n=0, random, fixed, variance,  ...) {

    time0 <- proc.time()    #debugging line
    ties <- match.arg(ties)
    Call <- match.call()






}   
