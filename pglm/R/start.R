starting.values <- function(family, link, vlink, rn, model, Kw, X, y, id, cl, start, other){
  # the relevant length for the start vector provided by the user is
  # either : 0 (nothing), K (the length of the beta vector or the
  # number of the additional parameter (1 most of the time, but 2 for
  # the negbin model)
  ls <- length(start)
  K <- ncol(X)
  
  if (family == "binomial"){
    if (model == "pooling"){
      # use Phi^-1(ybar) as starting value for the intercept and 0 for
      # the other coefficients
      if (!ls %in% c(0, K)) stop("irrelevant length for the start vector")
      if (ls == 0){
        interc <- ifelse(link == "probit", qnorm(mean(y)), - log(1 / mean(y) - 1))
        start <- c(interc, rep(0, ncol(X) - 1))
        names(start) <- colnames(X)
      }
    }
    if (model == "random"){
      if (!ls %in% c(0, 1, K + 1)) stop("irrelevant length for the start vector")
      if (ls <= 1){
        # the case where start is null or of length one (the sigma
        # parameter). In these cases, compute the pooling model
        startcl <- cl
        startcl[c('model', 'method', 'print.level')] <- c('pooling', 'nr', 0)
        glmest <- eval(startcl, parent.frame())
        if (ls == 0){
          # when start is null compute the starting value of sigma
          sigma <- lnl.binomial(coef(glmest), y = y, X = X, id = id,
                                link = link, model = "pooling", start.sigma = TRUE)
          start <- c(coef(glmest), sigma = sigma)
        }
        else start <- c(coef(glmest), sigma = start)
      }
    }
  }

  if (family == "ordinal"){
    J <- length(unique(y))
    if (model == "pooling"){
      if (!ls %in% c(0, K + J - 2)) stop("irrelevant length for the start vector")
      if (ls == 0){
        fy <- prop.table(table(y))
        fyc <- qnorm(cumsum(fy))
        beta0 <- - fyc[1]
        sup.coef <- fyc[2:(J-1)] + beta0
        start <- c(beta0, rep(0, K-1), sup.coef)
      }
      names(start) <- c(colnames(X), paste("mu", 1:(J - 2), sep = "_"))
    }
    if (model == "random"){
      sup.coef.names <- "sigma"
      if (!ls %in% c(0, 1, K + J - 1)) stop("irrelevant length for the start vector")
      if (ls <= 1){
        # the case where start is null or of length one (the sigma
        # parameter). In these cases, compute the pooling model
        startcl <- cl
        startcl[c('model', 'method', 'print.level')] <- c('pooling', 'nr', 0)
        glmest <- eval(startcl, parent.frame())
        if (ls == 0){
          # when start is null compute the starting value of sigma
          sigma <- lnl.ordinal(coef(glmest), y = y, X = X, id = id,
                               link = link, model = "pooling", start.sigma = TRUE)
          start <- c(coef(glmest), sigma = sigma)
        }
        else start <- c(coef(glmest), sigma = start)
      }
    }
  }
    
  if (family == "poisson"){
    if (model == "pooling"){
      # use Phi^-1(ybar) as starting value for the intercept
      if (!ls %in% c(0, K)) stop("irrelevant length for the start vector")
      if (ls == 0){
        start <- c(log(mean(y)), rep(0, K - 1))
        names(start) <- colnames(X)
      }
    }
    else{
      if (model == "within") relevant.values <- c(0, Kw)
      else relevant.values <- c(0, 1, K + 1)
      if (!ls %in% relevant.values) stop("irrelevant length for the start vector")
      if (ls <= 1){
        # the case where start is null or of length one (the sigma
        # parameter). In these cases, compute the pooling model
        startcl <- cl
        startcl$start <- NULL
        startcl[c('model', 'method', 'print.level')] <- c('pooling', 'nr', 0)
        glmest <- eval(startcl, parent.frame())
        if (ls == 0){
          if (model == "within") start <- coef(glmest)[Kw]
          else{
            if (FALSE){
              big <- 100
              init.val <- ifelse(other == "inv", big, 1/big)
              haty <- lnl.poisson(coef(glmest), y = y, X = X, id = id,
                                  link = link, model = "pooling")
              haty <- attr(haty, "fitted.values")
              print(head(haty))
              alphait <- y / haty
              alphai <- tapply(alphait, id, mean)
              print(head(alphai))
              print(var(alphai));stop()
              sigma <- lnl.poisson(c(coef(glmest), init.val), y = y, X = X, id = id,
                                   link = link, model = model)
              g <- apply(attr(sigma, "gradient"), 2, sum)
              g <- g[length(g)]
              h <- attr(sigma, "hessian")
              h <- h[nrow(h), nrow(h)]
              sigma <- init.val - g / h
              start <- c(coef(glmest), sigma = sigma)
            }
            start <- c(coef(glmest), sigma = 1.2)
          }
        }
        else start <- c(coef(glmest), sigma = start)
      }
    }
  }
    
  if (family == "negbin"){
    if (model == "pooling"){
      # Use the Poisson model for starting values
      startcl <- cl
      startcl[c('model', 'method', 'print.level', 'family')] <- c('pooling', 'nr', 0, poisson)
      glmest <- eval(startcl, parent.frame())
      haty <- lnl.poisson(coef(glmest), y = y, X = X, id = id,
                          link = link, model = "pooling")
      haty <- attr(haty, "fitted.values")
      res <- y - haty
      V <- var(res)
      E <- mean(y)
      if (vlink == 'nb1'){
        alpha <- V / E - 1
      }
      if (vlink == 'nb2'){
        alpha <- V / E^2 - 1 / E
##         print(alpha)
##         print(summary(lm(I(res^2/haty) ~ haty )))
##         print(summary(lm(I(res^2/haty) ~ haty -1)))
##         alpha <- coef(lm(I(res^2/haty)~ haty - 1)) - 1
##         print(alpha)
      }
      start <- c(coef(glmest), sigma = alpha)
    }
    else{
      if (model != "within") relevant.values <- c(0, K + 1L)
      else relevant.values <- c(0, 2, K + 3L)
      if (!ls %in% relevant.values) stop("irrelevant length for the start vector")
      if (ls <= 1){
        # the case where start is null or of length one (the sigma
        # parameter). In these cases, compute the pooling model
        startcl <- cl
        startcl[c('model', 'method', 'print.level', 'family')] <- c('pooling', 'nr', 0, poisson)
        glmest <- eval(startcl, parent.frame())
        if (ls == 0){
          if (model == "within") start <- coef(glmest)#[within.var]
          else start <- c(coef(glmest), a = 2, b = 1)
        }
        else start <- c(coef(glmest), a = 2, b = 1)
      }
    }
  }
    
  if (family == "gaussian"){
    if (is.null(other)) other <- "sd"
    if (is.null(start)){
      if (FALSE){
        startcl <- cl
        startcl[[1]] <- as.name("plm")
        startcl$family<- NULL
        startcl$model <- "random"
        glmest <- eval(startcl, parent.frame())
        theta <- glmest$ercomp$theta
        sigma2 <- glmest$ercomp$sigma2$idios
        sigmu <- glmest$ercomp$sigma2$id^0.5
        sigeps <- glmest$ercomp$sigma2$idiod^0.5
      }
      else{
        m <- match(c("formula", "data", "subset", "na.action"),names(cl),0)
        lmcl <- cl[c(1,m)]
        lmcl[[1]] <- as.name("lm")
        lmcl <- eval(lmcl, parent.frame())
        eb <- tapply(resid(lmcl), id, mean)[as.character(id)]
        sig2mu <- var(eb)
        ew <- resid(lmcl) - eb
        sig2eta <- var(ew)
        gamma <- sig2mu / sig2eta
      }
      if (other == "sd"){
        start <- c(coef(lmcl), sd.mu = sqrt(sig2mu), sd.eps = sqrt(sig2eta))
      }
      else{
        start <- c(coef(lmcl), gamma = gamma, sig2eta = sig2eta)
      }
    }
  }

  if (family == "tobit"){
    if (model == "pooling"){
      if (!ls %in% c(0, K + 1)) stop("irrelevant length for the start vector")
      if (ls == 0){
        m <- match(c("formula", "data", "subset", "na.action"),names(cl),0)
        lmcl <- cl[c(1,m)]
        lmcl[[1]] <- as.name("lm")
        lmcl <- eval(lmcl, parent.frame())
        sig2 <- deviance(lmcl) / df.residual(lmcl)
        if (is.null(other)) other <- "sd"
        if (other == "var") sigma <- sig2
        if (other == "sd") sigma <- sqrt(sig2)
        if (other == "lsd") sigma <- log(sqrt(sig2))
        start <- c(coef(lmcl), sd.eps = sigma)
      }
    }
    else{
      if (ls <= 1){
        startcl <- cl
        startcl$model <- "pooling"
        pglmest <- eval(startcl, parent.frame())
        thestart <- coef(pglmest)
        if (ls == 1){
          start <- c(thestart, start)
        }
        else{
          sigma <- lnl.tobit(coef(pglmest), y = y, X = X, id = id,
                               link = link, model = "pooling", start.sigma = TRUE)
          start <- c(thestart, sd.mu = sigma)
        }
      }
    }
  }
  start
}
