"power.grouped" <-
function(n = NULL, m, X = NULL, theta, sigma, sign.level = 0.05, 
                            type.power = c("conditional", "marginal"), MC.iter = 10,
                            type.lik = c("original", "approximate"), gr.mech = c("rounding", "equispaced"), 
                            dist.t, dist.x, power. = NULL, limits = c(10, 3000)){
    if(sum(sapply(list(n, power.), is.null)) != 1)
        stop("exactly one of the ", sQuote("n"), " and ", sQuote("power."), " must be NULL")
    if(is.null(n) && type.power == "conditional")
        stop(sQuote("type.power"), " must be ", sQuote("marginal"), ": sample size calculations are based on the marginal power formulas")
    if(is.null(n) && !is.null(X)) n <- nrow(X)
    if(missing(dist.x) && length(theta) != 2 && type.power == "marginal")
        stop("distributional assumptions about the extra covariates must be made")
    if(is.null(X) && missing(dist.t))                           
        stop("exactly one of the ", sQuote("X"), " and ", sQuote("dist.t"), " must be given")        
    if(!is.null(sign.level) && !is.numeric(sign.level) || any(sign.level < 0 | sign.level > 1))
        stop(sQuote("sign.level"), " must be numeric in [0, 1]")                    
    if(!is.null(X) && nrow(X) != n) 
        stop(sQuote("n"), " does not equal the number of the rows of X")
    gr.mech <- match.arg(gr.mech)   
    type.power <- match.arg(type.power)        
    type.lik <- match.arg(type.lik)          
    if(!is.null(X) && type.power == "marginal")
        stop("the type of power calculation cannot be marginal when ", sQuote("X"), " is given")
    if(is.null(X) && type.power == "conditional")
        stop("the type of power calculation cannot be conditional")    
    tpower <- switch(type.power, conditional = 1, marginal = 2)
    tlik <- switch(type.lik, original = 1, approximate = 2)
    pow. <- if(tpower == 1 && tlik == 1)
                quote({cond.power.orig(X = X, theta = theta, sigma = sigma, m = m, 
                    a = sign.level, grouping.mech = gr.mech)}) else
            if(tpower == 1 && tlik == 2)
                quote({cond.power.approx(X = X, theta = theta, sigma = sigma, m = m,
                    a = sign.level, grouping.mech = gr.mech)}) else
            if(tpower == 2 && tlik == 1)
                quote({marg.power.orig(MC = MC.iter, n = n, m = m, theta = theta, sigma = sigma,
                    a = sign.level, dist.t = dist.t, dist.x = dist.x, grouping.mech = gr.mech)[[1]]}) else
            if(tpower == 2 && tlik == 2)
                quote({marg.power.approx(MC = MC.iter, n = n, m = m, theta = theta, sigma = sigma, 
                    a = sign.level, dist.t = dist.t, dist.x = dist.x, grouping.mech = gr.mech)[[1]]})
    nn <- if(tlik ==1 ) quote({uniroot(samplesize.orig, interval = limits, p. = power., MC = MC.iter, theta = theta, 
                                         sigma = sigma, m = m, a = sign.level, dist.t = dist.t, dist.x = dist.x, 
                                         grouping.mech = gr.mech, tol = 1e-01, maxiter = 20)$root}) else
          if(tlik == 2) quote({uniroot(samplesize.approx, interval = limits, p. = power., MC = MC.iter, theta = theta, 
                                         sigma = sigma, m = m, a = sign.level, dist.t = dist.t, dist.x = dist.x, 
                                         grouping.mech = gr.mech, tol = 1e-01, maxiter = 20)$root}) 
    if(is.null(power.)) power. <- eval(pow.) 
    if(is.null(n)) n <- eval(nn)
    structure(list(n = round(n), m = m, theta = theta, sigma = sigma, sign.level = sign.level, 
                    type.power = type.power, type.lik = type.lik, power. = power., grouping.mech = gr.mech), 
			 	class = "power.grouped")   
              
}
