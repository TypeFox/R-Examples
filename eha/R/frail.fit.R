frail.fit <- function(X,
                      Y, rs, strats, offset, init, max.survs,
                      frailty, control
                      ){

    ## NOTE: Without intercept!
    ## Parameters: 1,..., n.rs first, then 1,..., p.
    nn <- NROW(X)
    ncov <- NCOL(X)
    
    if (missing(strats) || is.null(strats)) 
        strats <- rep(1, nn)
    
    if (missing(rs) || is.null(rs)){
        rs <- risksets(Y, strata = strats, max.survs)
    }

    if (max(rs$riskset) > nn) stop("Riskset does not match data")
    
    weg <- (abs(rs$size - rs$n.events) > 0.01) ## They are integers ?!!
    rs$riskset <- rs$riskset[rep(weg, rs$size)]
    rs$eventset <- rs$eventset[rep(weg, rs$n.events)]
    rs$n.events <- rs$n.events[weg]
    rs$size <- rs$size[weg]

    n.rs <- length(rs$size)

    ev <- numeric(sum(rs$size))

    start <- 1

    for (i in 1:n.rs){
        ev[start:(start + rs$n.events[i] - 1)] <- 1
        start <- start + rs$size[i]
    }
    
    haz <- rep(1:length(rs$size), rs$size)

    frailty <- frailty[rs$riskset]

    frailty <- as.integer(frailty)

    #ord <- order(frailty + seq(0, 0.9, length = length(frailty))) #integers
    # This is faster (?!?):
    ord <- sort(frailty, method = "quick", index.return = TRUE)$ix

    fam.size <- table(frailty)

    haz <- haz[ord]
    ev <- ev[ord]  # Binary Response
    rs$riskset <- rs$riskset[ord]

    n <- length(ev)
    
    if (missing(offset) || is.null(offset)) 
        offset <- rep(0, nn)
    
    if (missing(init) || is.null(init)) 
        init <- rep(0, ncov)
    
    if (missing(control)){
        control <- list(eps=1.e-8,
                        maxiter = 10,
                        n.points = 12,
                        trace = FALSE)
    }else{
        if (!is.numeric(control$eps)){
            stop("Error in control = list(eps = ...) ")
        }else{
            if (control$eps <= 0) stop("control$eps must be strictly positive")
        }
        if (!is.numeric(control$maxiter)){
            stop("Error in control = list(maxiter = ...) ")
        }else{
            if (as.integer(control$maxiter) <= 0)
                stop("control$maxiter must be positive")
        }
        if (!is.numeric(control$n.points)){
            stop("Error in control = list(n.points = ...) ")
        }else{
            if (as.integer(control$n.points) <= 0)
                stop("control$maxiter must be positive")
        }
        if (!is.logical(control$trace)) stop("control$trace must be logical")
    }

    n.rs <- max(haz)
    p <- ncov
    start.coef <- numeric(p + n.rs)
    start.coef[1:n.rs] <- log(-log(1 - rs$n.events / rs$size))

    start.coef[(n.rs + 1):(n.rs + p)] <- init

    printlevel <- control$trace
      ## NOTE: silent == TRUE ===> printlevel = 0
    iter <- control$maxiter

    fam <- 1 ##logit(link = cloglog)
    method <- 1
    start.sigma <- 1
    n.fam <- length(fam.size)

    X <- scale(X, center = TRUE, scale = FALSE)
    if (control$trace)
        cat("Go into [frail_ml]\n")
    bdim = n.rs + p + 1
    fit <- .C("frail_ml",
              as.integer(fam), ## logit(cloglog)
              as.integer(method), ## "vmmin"
              as.integer(p), #leading dimension of t(X)
              as.integer(nn), #second dimension of t(X)
              as.integer(n.rs), # No. of risksets (extra parameters).
              as.integer(rs$riskset - 1),
              as.double(start.coef), 
              as.double(start.sigma),
              as.double(t(X)),       ### Note CAREFULLY (03-01-09)!!!
              as.integer(ev),    ## "Y"
              as.integer(haz - 1),  ## New, for hazard contributions
              as.double(offset),
              as.integer(fam.size),
              as.integer(n.fam),
              as.integer(control$n.points),  ## For Gauss-Hermite
              as.double(control$eps),
              as.integer(control$maxit),
              as.integer(control$trace),
              beta = double(p),  ## Return values from here.
              hazards = double(n.rs),
              sigma = double(1),
              sigma.sd = double(1),
              loglik = double(2), ## Changed from 1 to 2, May 8, 2003
              variance = double(bdim * bdim),
              frail = double(n.fam),
           ##   mu = double(n),
              convergence = integer(1),
              fail = integer(1),
              ## DUP = FALSE,
              PACKAGE = "eha"
              )  

    if (fit$fail){
        warning("Singular hessian (returned); no S.E. estimates.")
    }

    var <- matrix(fit$variance, ncol = bdim)[(n.rs + 1):(n.rs + p),
                                (n.rs + 1):(n.rs + p), drop = FALSE]

    list(coefficients = fit$beta,
         var = var,
         sigma = fit$sigma,
         sigma.sd = fit$sigma.sd,
         loglik = fit$loglik,
         ##score = fit$sctest,
         linear.predictors = X %*% fit$beta,
         residuals = resid,
         hazards = fit$hazards,
         frailty = fit$frail,
         means = apply(X, 2, mean),
         conver = fit$conver,
         fail = fit$fail,
         iter = fit$iter
       )
       
  }
