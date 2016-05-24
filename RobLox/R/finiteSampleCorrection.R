###############################################################################
## Function for finite-sample correction of the neighborhood radius
###############################################################################
finiteSampleCorrection <- function(r, n, model = "locsc"){
    if(model == "locsc" & r >= 1.74) return(r)
    if(model %in% c("loc", "sc") & r >= 3.0) return(r)
    if(n == 1) return(Inf)
    if(n == 2) return(Inf)

    eps <- r/sqrt(n)
    ns <- c(3:50, seq(55, 100, by = 5), seq(110, 200, by = 10), 
            seq(250, 500, by = 50))
    epss <- c(seq(0.001, 0.01, by = 0.001), seq(0.02, to = 0.5, by = 0.01))
    if(n %in% ns){
        ind <- ns == n
    }else{
        ind <- which.min(abs(ns-n))
    }
    if(model == "locsc")
        return(max(r, approx(x = epss, y = .finiteSampleRadius.locsc[,ind], xout = eps, rule = 2)$y))
    if(model == "loc")
        return(max(r, approx(x = epss, y = .finiteSampleRadius.loc[,ind], xout = eps, rule = 2)$y))
    if(model == "sc")
        return(max(r, approx(x = epss, y = .finiteSampleRadius.sc[,ind], xout = eps, rule = 2)$y))
    else
        stop("argument 'model' has to be 'locsc', 'loc' or 'sc'")
}
