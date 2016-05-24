dose.distr <-
  function(f, pars, beta, cov, cells, dics, m.prior="gamma", d.prior="uniform", prior.param=c(0,"Inf"), stdf=6)
  { 
    pos   <- 1
    envir <- as.environment(pos)
    if (class(f)!="expression") stop("Dose response curve should be an expression")
    if (length(pars)!=length(beta)) stop("Wrong number of parameters")
    if (d.prior!="gamma" & d.prior!="uniform") stop("Wrong dose prior distribution")
    if (m.prior!="normal" & m.prior!="gamma") stop("Wrong mean prior distribution")
    if (d.prior=="uniform" & (length(prior.param)!=2 | prior.param[2] <= prior.param[1])) stop ("Wrong parameter specification for uniform prior")
        
    for (i in 1:length(beta)) 
    {
      assign(eval(parse(text = paste0("pars[", i, "]"))), beta[i], envir=envir)
    }
    
    assign("u.local2", u(f, x, beta), envir=envir)
    assign("v.local2", v(f, x, beta), envir=envir)

    u.local <- get("u.local2", envir=envir)
    v.local <- get("v.local2", envir=envir)

    # Doses sequence for the calibrative density
    pe <- uniroot(function(x) u.local(x) - dics/cells, c(.001, 100), extendInt="yes")$root
    nsim <- 1000
    if (d.prior=="uniform" & prior.param[2]!="Inf") rand <- runif(nsim, prior.param[1], prior.param[2])
    if (d.prior=="uniform" & prior.param[2]=="Inf") rand <- runif(nsim, 0, 50)
    if (d.prior=="gamma")   rand <- rgamma(nsim, prior.param[1]^2/prior.param[2]^2, prior.param[1]/prior.param[2]^2)
    s <- sum(sapply(1:length(rand), function(i) dnbinom(dics, u.local(rand[i])^2/diag((attr(v.local[[1]](rand[i])(rand[i]), "gradient") %*% cov %*% t(attr(v.local[[1]](rand[i])(rand[i]), "gradient")))), 1/(1 + cells*diag((attr(v.local[[1]](rand[i])(rand[i]), "gradient") %*% cov %*% t(attr(v.local[[1]](rand[i])(rand[i]), "gradient"))))/u.local(rand[i])))))
    expv <- sum(sapply(1:length(rand), function(i) rand[i]*dnbinom(dics, u.local(rand[i])^2/diag((attr(v.local[[1]](rand[i])(rand[i]), "gradient") %*% cov %*% t(attr(v.local[[1]](rand[i])(rand[i]), "gradient")))), 1/(1 + cells*diag((attr(v.local[[1]](rand[i])(rand[i]), "gradient") %*% cov %*% t(attr(v.local[[1]](rand[i])(rand[i]), "gradient"))))/u.local(rand[i])))))/s 
    exv2 <- sum(sapply(1:length(rand), function(i) rand[i]*rand[i]*dnbinom(dics, u.local(rand[i])^2/diag((attr(v.local[[1]](rand[i])(rand[i]), "gradient") %*% cov %*% t(attr(v.local[[1]](rand[i])(rand[i]), "gradient")))), 1/(1 + cells*diag((attr(v.local[[1]](rand[i])(rand[i]), "gradient") %*% cov %*% t(attr(v.local[[1]](rand[i])(rand[i]), "gradient"))))/u.local(rand[i])))))/s
    stde <- sqrt(exv2 - expv^2)
    xl <- pe - stdf*stde
    if (d.prior=="uniform" & xl < prior.param[1]) 
    {
      xl <- as.numeric(prior.param[1])
    }
    xu <- pe + stdf*stde
    if (prior.param[2] != Inf)
    {
      if (d.prior=="uniform" & xu > prior.param[2]) 
      {
        xu <- as.numeric(prior.param[2])
      }
    }
    if (xl < 0) xl <- 0
    if (xu < 0) xu <- 0
    x  <- seq(xl, xu, 0.001)
    l  <- length(x)
    cd <- numeric(l)
    
    if (m.prior=="gamma")
    {
      # Calibrative density calulation for the gamma mean prior and Uniform dose prior 
      pr <- attr(v.local[[1]](x)(x), "gradient")
      cd <- sapply(1:l, function(i) dnbinom(dics, a(x[i]), b(x[i])/(b(x[i]) +  cells)))
      cd <- approxfun(x, cd)
      cnorm <- integrate(cd, lower=xl, upper=xu)$value
      cd <- cd(x)/cnorm
      if (d.prior=="gamma")
      {
        cd <- cd*dgamma(x, prior.param[1]^2/prior.param[2]^2, prior.param[1]/prior.param[2]^2)
        cd <- approxfun(x, cd)
        cnorm <- integrate(cd, lower=xl, upper=xu)$value
        cd <- cd(x)/cnorm
      }
    }
    if (m.prior=="normal")
    {
      pr <- attr(v.local[[1]](x)(x), "gradient")
      cd <- sapply(1:l, function(i) dhermite(dics, cells*u.local(x[i])*(2-(1 + cells*(pr[i,] %*% cov %*% as.matrix(pr[i,]))/u.local(x[i]))), cells*u.local(x[i])*(cells*(pr[i,] %*% cov %*% as.matrix(pr[i,]))/u.local(x[i]))/2, 2))
      cd <- approxfun(x, cd)
      cnorm <- integrate(cd, lower=xl, upper=xu)$value
      cd <- cd(x)/cnorm
      if (d.prior=="gamma")
      {
        cd <- cd*dgamma(x, prior.param[1]^2/prior.param[2]^2, prior.param[1]/prior.param[2]^2)
        cd <- approxfun(x, cd)
        cnorm <- integrate(cd, lower=xl, upper=xu)$value
        cd <- cd(x)/cnorm
      }
    }
    res <- list()
    res[[1]] <- cd
    res[[2]] <- x
    attr(res,"class") <- "dose.radir"
    return(res)
  }