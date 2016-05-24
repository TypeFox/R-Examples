"mctest" <-
function(x, distn, parm, H, sim, tol, STATISTIC, estfun) 
{
    qdist <- get(paste("q", substring(distn, 2), sep = ""), mode = "function")
    pdist <- get(distn, mode = "function")    
    cdist <- cdens(paste("p", substring(distn, 2), sep = ""), H)
    rdist <- function(m, parm){
        FH <- do.call("pdist", c(list(H), parm))
        do.call("qdist", c(list(runif(m) * (1-FH) + FH), parm))
    }
    TS <- function(zH, z) 
        STATISTIC(n, zH, z, j)
            
    ll0 <- function(p, v) {
        names(p) <- nm
        d  <- do.call("cdist", c(list(v), p))
        res <- -sum(log(d[d>0]))
        if (is.na(res)) res <- errvalue
        res
    }   
     
    ll <- function(p, v) tryCatch(ll0(p, v), 
            warning = function(w) errvalue, error = function(e) errvalue)
    
    S <- function(v) {                    
        if (is.estfun) {
            x <- v
            fit <- eval(parse(text = estfun))
        }
        else {   
            est <- optim(parm, ll, method = "BFGS", hessian = FALSE, v = v)$par
            fit <- as.list(est)}
        
        zH  <- do.call("pdist", c(H, fit))   
        zj  <- do.call("pdist", c(list(sort(v)), fit))            
        TS(zH, zj)
    }
    
    is.estfun <- ifelse(!is.na(estfun), TRUE, FALSE)
    if (min(x) < H) 
       stop("'min(x)' must be greater or equal to 'H'")
    nm <- names(parm) 
    n <- length(x)
    j <- c(1:n)
    errvalue <- -log(.Machine$double.xmin) * n
    TS0 <- TS(do.call("pdist", c(list(H), parm)), 
        do.call("pdist", c(list(sort(x)), as.vector(parm))))
    if (!is.finite(TS0))
       stop("test statisic value can't be calculated")
    sz <- sapply(1:sim, function(y) rdist(n, parm))
    pval <- .Call("pvalue", as.list(1:sim), function(k) S(sz[,k]), new.env(), TS0, tol)
    return(list(p.value = pval[2], TS = TS0, sim = pval[1]))  
}
