iccbin <- function(n, y, data, method = c("A", "B", "C"), nAGQ = 1, M = 1000) {
    CALL <- match.call()
    method <- match.arg(method)

    resp <- c(deparse(substitute(n)), deparse(substitute(y)))
    datan <- data[ , resp]
    if(any(datan[ , 2] > datan[ , 1]))
        stop("Some ", deparse(substitute(y)), " were > ", deparse(substitute(n)), ".")
    if(any(datan[ , 1] <= 0))
        warning("Data with ", deparse(substitute(n)), " <= 0 were discarded.")
    names(datan) <- c("n", "y")

    datan <- datan[datan$n > 0, ]
    databin <- splitbin(cbind(y, n - y) ~ 1, datan)$tab

    # Binomial-Gaussian glmm model
    if(method %in% c("A", "B")) {
##        require(lme4)
        fm <- lme4::glmer(formula = y ~ 1 + (1 | idbin), family = binomial, data = databin, nAGQ = nAGQ)
        b0 <- lme4::fixef(fm)
        s2_u <- as.vector(lme4::VarCorr(fm)[[1]])
        if(method == "A") {
            # method A - model linearization
            p <- exp(b0) / (1 + exp(b0))
            nu1 <- p * (1 - p)
            nu2 <- s2_u * p^2 * (1 + exp(b0))^(-2)
            M <- NA
            } else {
            # method B - Monte Carlo
                u <- rnorm(n = M, mean = 0, sd = sqrt(s2_u))
                p <- exp(b0 + u) / (1 + exp(b0 + u))
                nu1 <- mean(p * (1 - p))
                nu2 <- var(p)
                }
        rho <- nu2 / (nu1 + nu2)
        }

    # Gaussian lmm model
    # Moments (ANOVA)
    if(method == "C") {
        k <- nrow(datan) ; n <- datan$n ; N <- sum(n)
        nA <- (N - sum(n^2) / N) / (k - 1)
        fm <- lm(y ~ idbin, data = databin)
        a <- anova(fm)
        MSA <- a[1, 3]
        MSE <- a[2, 3]
        rho <- (MSA - MSE) / (MSA + (nA - 1) * MSE)
        #ftest <- c(F.value = a[1, 4], df.num = a[1, 1], df.denom = a[2, 1], P = a[1, 5])
        #rho <- (ftest[1] - 1) / (ftest[1] + nA - 1)
        rho <- as.numeric(rho)
        nAGQ <- M <- NA
        }

    ## output
    features <- c(method = method, nAGQ = nAGQ, M = M)
    new(Class = "iccbin", CALL = CALL, features = features, rho = rho)
    }
