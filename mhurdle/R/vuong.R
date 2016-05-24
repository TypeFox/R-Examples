# simulated probabilities and quantiles for the weighted chi-squares
# distribution
pwchisq <- function(q, weights, lower.tail = TRUE, n = 1000){
    set.seed(100)
    K <- length(weights)
    e <- matrix(rnorm(n * K) ^ 2, n, K)
    wcs <- apply(e, 1, function(x) sum(x * weights))
    F <- ecdf(wcs)
    ifelse(lower.tail, F(q), 1 - F(q))
}
qwchisq <- function(p, weights, lower.tail = TRUE, n = 1000){
    K <- length(weights)
    e <- matrix(rnorm(n * K) ^ 2, n, K)
    wcs <- apply(e, 1, function(x) sum(x * weights))
    ifelse(lower.tail, quantile(wcs, p), quantile(wcs, 1 - p))
}  

vuongtest <- function(x, y,
                      type = c("non-nested", "nested", "overlapping"),
                      hyp = FALSE,
                      variance = c("centered", "uncentered"),
                      matrix = c("large", "reduced")
                      ){
    type <- match.arg(type)
    variance <- match.arg(variance)
    matrix <- match.arg(matrix)
    
    if (matrix == "reduced" && type != "nested")
        stop("the reduced matrix is only relevant for nested models")
    
    data.name <- c(
        paste(deparse(substitute(x))),
        paste(deparse(substitute(y)))
        )
    set.seed(100)
    
    ###  for convenience, call f the larger model and g the other one if
    ###  the models are nested, else call the first model f and g for
    ###  consistency with vuong paper. Check also that the models are
    ###  really nested, otherwise break
    f <- x
    g <- y
    if (type == "nested"){
        if (length(coef(x)) < length(coef(y))){
            f <- y
            g <- x
        }
        nestedOK <- prod(names(coef(g)) %in% names(coef(f)))
        if (! nestedOK || (length(coef(f)) == length(coef(g)) ))
            stop("the two models are not nested")
    }
    
    lf <- as.numeric(f$logLik) ; lg <- as.numeric(g$logLik)
    logLf <- sum(lf)           ; logLg <- sum(lg)
    Kf <- length(coef(f))      ;  Kg <- length(coef(g))
    n <- nrow(model.frame(f))
    
    if (nrow(model.frame(g)) != n) stop("the number of observations of the two models differ")
    
    gradf <- f$gradient
    gradg <- g$gradient
    Bf <- crossprod(gradf) / n
    Bg <- crossprod(gradg) / n
    Bfg <- t(gradf) %*% gradg / n
    Af1 <- - vcov(f) * n
    Ag1 <- - vcov(g) * n
  
    ### compute the likelihood ratio
    LR <- logLf - logLg
  
    ### compute the variance
    w2 <- ifelse(variance == "centered",
                 1 / n * sum( (lf - lg) ^ 2) - (1 / n * LR) ^ 2,
                 1 / n * sum( (lf - lg) ^ 2)
                 )
    

    #### construct the large or reduced matrix and its eigen values
    if (matrix == "large"){
        W <- rbind(cbind( -     Bf %*% Af1,  - Bfg %*% Ag1),
                   cbind(   t(Bfg) %*% Af1,   Bg  %*% Ag1)
                   )
        Z <- eigen(W)$values
    }
    else{
        common.coef <- names(coef(f)) %in% names(coef(g))
        D <- t(diag(1, Kf)[common.coef, ])
        W <- Bf %*% (D %*% Ag1 %*% t(D) - Af1)
        Z <- eigen(W)$values
    }

    ### non nested test ; only the version with wrong specification
    ### hypothesis is implemented
    if (type == "non-nested"){
        if (hyp) stop("this non-nested test is not implemented")
        statistic <- c(z = LR / sqrt(n * w2))
        method <- "Vuong Test (non-nested)"
        # if the stat is negative, the p-value is the lower tail,
        # otherwise it is the upper tail
        lower.tail <- statistic < 0
        pval <- pnorm(statistic, lower.tail = lower.tail)
        parameter <- NULL
    }

    ### nested test
    if (type == "nested"){
        method <- "Vuong Test (nested)"
        statistic <- 2 * LR
        if (! hyp){
            parameter <- c(sev = sum(Z))
            names(statistic) <- "wchisq"
            pval <- pwchisq(statistic, Z, lower.tail = FALSE)
        }
        else{
            parameter <- c(df = Kf - Kg)
            names(statistic) <- "chisq"
            pval <- pchisq(statistic, parameter, lower.tail = FALSE)
        }
    }
    
    #### overlapping test
    if (type == "overlapping"){
        method <- "Vuong Test (overlapping)"
        if (! hyp){
        ### test first the hypothesis that w^2=0
            statistic <- c(wchisq = n * w2)
            pval <- pwchisq(statistic, Z ^ 2 , lower.tail = FALSE)
            parameter <- c(sev = sum(Z ^ 2))
        }
        else{
        ### In this case, the LR statistic can be either positive or
        ### negative, depending on the order of the models. The test is
        ### two-sided so that the p-value is twice the lower tail if the
        ### statistic is negative and twice the upper tail otherwise
            statistic <- c(wchisq = 2 * LR)
            lower.tail <- statistic < 0
            pval <- 2 * pwchisq(statistic, Z, lower.tail = lower.tail)
            parameter <- c(sev = sum(Z))
        }
    }
    if (length(data.name) > 1) data.name <- paste(data.name, collapse = "-")
    result <- list(statistic = statistic,
                   method = method,
                   p.value = pval,
                   data.name = data.name,
                   parameter = parameter)
    class(result) <- "htest"
    result
}

