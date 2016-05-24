# Maximum likelihood based on the full log-likelihood (personal level)
ml.gb2 <- function (z, w=rep(1, length(z)), method = 1, hess = FALSE){
fn <- function(x, z, w) {
        a <- x[1]
        b <- x[2]
        p <- x[3]
        q <- x[4]
        return(-loglp.gb2(z, a, b, p, q, w))
    }
    gr <- function(x, z, w) {
        a <- x[1]
        b <- x[2]
        p <- x[3]
        q <- x[4]
        return(-scoresp.gb2(z, a, b, p, q, w))
    }
    
    x0 <- fisk(z, w)
    
    opt1 <- optim(x0, fn, gr, z, w, method = "BFGS", control = list(parscale = x0, 
        pgtol = 1e-08), hessian = hess)
    if (method != 2) 
        return(list(opt1 = opt1))
    if (method == 2) {
        opt2 <- optim(x0, fn, gr, z, w, method = "L-BFGS-B", 
            lower = 0, control = list(parscale = x0, pgtol = 1e-08), 
            hessian = hess)
        return(list(opt2 = opt2))
    }
}

# Maximum likelihood based on the full log-likelihood (household level)
mlh.gb2 <- function (z, w=rep(1, length(z)), hs=rep(1, length(z)), method = 1, hess = FALSE) 
{
    fn <- function(x, z, w, hs) {
        a <- x[1]
        b <- x[2]
        p <- x[3]
        q <- x[4]
        return(-loglh.gb2(z, a, b, p, q, w, hs))
    }
    gr <- function(x, z, w, hs) {
        a <- x[1]
        b <- x[2]
        p <- x[3]
        q <- x[4]
        return(-scoresh.gb2(z, a, b, p, q, w, hs))
    }
    
    x0 <- fiskh(z, w, hs)
    
   opt1 <- optim(x0, fn, gr, z, w, hs, method = "BFGS", control = list(parscale = x0, 
        pgtol = 1e-08), hessian = hess)
    if (method != 2) 
        return(list(opt1 = opt1))
    if (method == 2) {
        opt2 <- optim(x0, fn, gr, z, w, hs, method = "L-BFGS-B", 
            lower = 0, control = list(parscale = x0, pgtol = 1e-08), 
            hessian = hess)
        return(list(opt2 = opt2))
    }
}