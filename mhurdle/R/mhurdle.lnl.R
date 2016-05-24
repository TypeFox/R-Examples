# This function computes the log-likelihood for the hurdles models and
# returns also as attributes, if required, the gradient, the fitted
# values and the score vector.

mhurdle.lnl <- function(param, X1, X2, X3, X4, y, gradient = FALSE,
                        fitted = FALSE, dist = NULL, corr = NULL){
  #  Extract the elements of the model
    h1 <- !is.null(X1) ;  K1 <- ifelse(is.null(X1), 0, ncol(X1))
    h3 <- !is.null(X3) ;  K3 <- ifelse(is.null(X3), 0, ncol(X3))
    K2 <- ncol(X2)
    beta2 <- param[(K1 + 1):(K1 + K2)]
    bX2 <- as.numeric(crossprod(t(X2), beta2))
    if (h1){
        beta1 <- param[1:K1]
        bX1 <- as.numeric(crossprod(t(X1), beta1))
        Phi1 <- pnorm(bX1) ; phi1 <- dnorm(bX1)
    }
    else{
        bX1 <- beta1 <- NULL
        Phi1 <- 1 ; phi1 <- 0;
    }
    if (h3){
        beta3 <- param[(K1 + K2 + 1):(K1 + K2 + K3)]
        bX3 <- as.numeric(crossprod(t(X3), beta3))
        Phi3 <- pnorm(bX3) ; phi3 <- dnorm(bX3)
    }
    else{
        bX3 <- beta3 <- NULL
        Phi3 <- 1 ; phi3 <- 0
    }
    if (is.null(X4)) K4 <- 1 else K4 <- ncol(X4)
    beta4 <- param[(K1 + K2 + K3 + 1):(K1 + K2 + K3 + K4)]
    if (is.null(X4)) sigma <- beta4 else sigma <- as.numeric(exp(crossprod(t(X4), beta4)))
    rho1 <- rho3 <- 0
    if (!is.null(corr)){
        if (corr == 'h1'){
            rho1 <- param[K1 + K2 + K3 + K4 + 1]
            if (rho1 < -1) rho1 <- - 0.99
            if (rho1 >  1) rho1 <-   0.99
        }
        if (corr == 'h3'){
            rho3 <- param[K1 + K2 + K3 + K4 + 1]
            if (rho3 < -1) rho3 <- - 0.99
            if (rho3 >  1) rho3 <-   0.99
        }
    }
    if (dist == "bc"){
        # in case of a box-cox distribution, the truncatation point is
        # bX2 + 1 / lambda < 0
        lambda <- param[K1 + K2 + K3 + K4 + (!is.null(corr)) + 1]
        Phi2 <- pnorm((bX2 + 1 / lambda) / sigma)
        phi2 <- dnorm((bX2 + 1 / lambda) / sigma)
        Phi12 <- mypbivnorm(bX1, (bX2 + 1 / lambda) / sigma , rho1)
        Phi23 <- mypbivnorm( (bX2 + 1 / lambda) / sigma, bX3, rho3)
    }
    else{
        Phi2 <- pnorm(bX2 / sigma)
        phi2 <- dnorm(bX2 / sigma)
        Phi12 <- mypbivnorm(bX1, bX2 / sigma , rho1)
        Phi23 <- mypbivnorm(bX2 / sigma, bX3, rho3)
    }
    if (dist == "ihs") lambda <- param[K1 + K2 + K3 + K4 + (!is.null(corr)) + 1]
    Ty <- switch(dist,
                 "ln" = log2(y) + log(Phi3),
                 "bc" = ((y * Phi3) ^ lambda - 1) / lambda,
                 "ihs" = log(lambda * y * Phi3 + sqrt(1 + (lambda  * y * Phi3) ^ 2)) / lambda,
                 y * Phi3
                 )
    # logarithm of the jacobian
    lnJ <- switch(dist,
                 "ln" = - log2(y),
                 "bc" = (lambda - 1) * log2(y) + lambda * log(Phi3),
                 "ihs" = - 0.5 * log(1 + (lambda * Phi3 * y) ^ 2) + log(Phi3),
                 log(Phi3)
                 )

    lnJlb <- switch(dist,
                    "bc" = log2(y) + log(Phi3),
                    "ihs" = - lambda * y ^ 2 * Phi3 ^ 2 / (1 + (lambda * y * Phi3) ^ 2)
                    )
    resid <- Ty - bX2
    
    z <- function(x, resid, rho){
        if (is.null(x)) result <- list(f = 100, g = 0, h = 0)
        else{
            f <- (x + rho / sigma * resid) / sqrt(1 - rho ^ 2) 
            g <- x * rho * (1 - rho ^ 2) ^ - 1.5 + ((1 - rho ^ 2) ^ - 0.5 + rho ^ 2 * (1 - rho ^ 2) ^ - 1.5) * resid / sigma
            h <- x * ((1 - rho ^ 2) ^ -1.5 + 3 * rho ^ 2 * (1 - rho ^ 2) ^ -2.5) +
                resid / sigma * (3 * rho * (1 - rho ^ 2) ^ -1.5 + 3 * rho ^ 3 * (1 - rho ^ 2) ^ -2.5)
            result <- list(f = f, g = g, h = h)
        }
        result
    }
    z1 <- z(bX1, resid, rho1)
    z3 <- z(bX3, resid, rho3)
    
    lnL.null <- switch(dist,
                       "tn" = log(1 - Phi12$f * Phi23$f / Phi2 ^ 2),
                       "ln" = log(1 - Phi1 * Phi3),
                       log(1 - Phi12$f * Phi23$f / Phi2)
                       )
    lnL.pos <-
        - log(sigma) +
            dnorm(resid / sigma, log = TRUE) +
                pnorm(z1$f, log.p = TRUE) +
                    pnorm(z3$f, log.p = TRUE) +
                        lnJ - (dist == "tn") * log(Phi2)
    
    lnL <- lnL.null * (y == 0) + lnL.pos * (y != 0)
    
    if (gradient){
        gradi <- c()
        if (h1){
            lnL.beta1 <- (y == 0)*(
                              switch(dist,
                                     "tn"  = - (Phi12$a * Phi23$f)/(Phi2 ^ 2 - Phi12$f * Phi23$f),
                                     "ln"  = - phi1 * Phi3 / (1 - Phi1 * Phi3),
                                     - (Phi12$a * Phi23$f)/(Phi2 - Phi12$f * Phi23$f))
                              ) +
                                  (y != 0) * (mills(z1$f) / sqrt(1 - rho1 ^ 2))
            gradi <- cbind(gradi, lnL.beta1 * X1)
        }
        lnL.beta2 <- (y == 0)*(
                          switch(dist,
                                 "tn" =  (2 * Phi2 * phi2 - Phi12$b * Phi23$f - Phi12$f * Phi23$a) /
                                 (Phi2 ^ 2 - Phi12$f * Phi23$f) / sigma - 2 * mills(bX2 / sigma) / sigma,
                                 "ln" = 0,
                                 "bc" = (phi2 - Phi12$b * Phi23$f - Phi12$f * Phi23$a) /
                                 (Phi2 - Phi12$f * Phi23$f) / sigma - mills((bX2 + 1 / lambda)/ sigma) / sigma,
                                 (phi2 - Phi12$b * Phi23$f - Phi12$f * Phi23$a) /
                                 (Phi2 - Phi12$f * Phi23$f) / sigma - mills(bX2 / sigma) / sigma)
                          ) +
                              (y != 0) * (
                                   resid / sigma ^ 2 -  mills(z1$f) * rho1 / sigma / sqrt(1 - rho1 ^ 2) -
                                   mills(z3$f) * rho3 / sigma / sqrt(1 - rho3 ^ 2) -
                                   (dist == "tn") * mills(bX2 / sigma) / sigma
                                   )
        gradi <- cbind(gradi, lnL.beta2 * X2)
        
        if (h3){
            Ty3 <- switch(dist,
                          "ln" = mills(bX3),
                          "bc" = (y * Phi3) ^ lambda  * mills(bX3),
                          "ihs" = y * phi3 / sqrt( 1 + (lambda * y * Phi3) ^ 2),
                          y * phi3
                          )    
            lnJ3 <- switch(dist,
                           "ln" = 0,
                           "bc" = lambda * mills(bX3),
                           "ihs" = - phi3 * Phi3 * lambda ^ 2 * y ^ 2 / (1 + (lambda * y * Phi3) ^ 2) + mills(bX3),
                           mills(bX3)
                          )

            lnL.beta3 <- (y == 0) * (
                              switch(dist,
                                     "tn" = - (Phi12$f * Phi23$b)/(Phi2 ^ 2 - Phi12$f * Phi23$f),
                                     "ln" = - Phi1 * phi3 / (1 - Phi1 * Phi3),
                                     - (Phi12$f * Phi23$b)/(Phi2 - Phi12$f * Phi23$f))
                              ) + 
                                  (y != 0) * (
                                       - resid / sigma ^ 2 * Ty3
                                       + mills(z1$f) * rho1 / sqrt(1 - rho1 ^ 2) / sigma * Ty3
                                       + mills(z3$f) * (1 + rho3 / sigma * Ty3) / sqrt(1 - rho3 ^ 2)
                                       + lnJ3
                                       )
            gradi <- cbind(gradi, lnL.beta3 * X3)
        }
        
        lnL.sigma <- (y == 0)*(
                          switch(dist,
                                 "tn" =  - (2 * Phi2 * phi2 - Phi12$b * Phi23$f - Phi12$f * Phi23$a) /
                                 (Phi2^2 - Phi12$f * Phi23$f) * bX2 / sigma^2 + 2 * mills(bX2 / sigma) * bX2 / sigma^2,
                                 "ln" = 0,
                                 "bc" = - (phi2 - Phi12$b * Phi23$f - Phi12$f * Phi23$a) /
                                 (Phi2 - Phi12$f * Phi23$f) * (bX2 + 1 / lambda) / sigma^2 +
                                 mills((bX2 + 1 / lambda) / sigma) * (bX2 + 1 / lambda)/ sigma^2,
                                 - (phi2 - Phi12$b * Phi23$f - Phi12$f * Phi23$a) /
                                 (Phi2 - Phi12$f * Phi23$f) * bX2 / sigma^2 + mills(bX2 / sigma) * bX2 / sigma^2
                                 )
                          ) +
                              (y != 0) * (
                                   - 1 / sigma + resid ^ 2 / sigma ^ 3 - mills(z1$f) * rho1 / sqrt(1 - rho1 ^ 2) * resid / sigma ^ 2 -
                                   mills(z3$f) * rho3 / sqrt(1 - rho3 ^ 2) * resid / sigma ^ 2 +
                                   (dist == "tn") * mills(bX2 / sigma) * bX2 / sigma ^ 2
                                   )
        if (!is.null(X4)) lnL.sigma <- lnL.sigma * sigma * X4
        gradi <- cbind(gradi, lnL.sigma)
        
        if (h1 && !is.null(corr) && corr == 'h1'){
            lnL.rho1 <- (y == 0) * (
                             switch(dist,
                                    "tn" = - Phi12$rho * Phi23$f / (Phi2^2 - Phi12$f * Phi23$f),
                                    "ln" = 0,
                                    - Phi12$rho * Phi23$f / (Phi2   - Phi12$f * Phi23$f)
                                    )
                             ) +
                                 (y != 0) * (
                                      mills(z1$f) * z1$g
                                      )
            gradi <- cbind(gradi, lnL.rho1)
        }
        
        if (h3 && !is.null(corr) && corr == 'h3'){
            lnL.rho3 <- (y == 0) * (
                             switch(dist,
                                    "tn" = - Phi12$f * Phi23$rho / (Phi2^2 - Phi12$f * Phi23$f),
                                    "ln" = 0,
                                    - Phi12$f * Phi23$rho / (Phi2 - Phi12$f * Phi23$f)
                                    )
                             ) +
                                 (y != 0) * (
                                      mills(z3$f) * z3$g
                                      )
            gradi <- cbind(gradi, lnL.rho3)
        }

        if (dist == "bc"){
            Tylb <- (log(Phi3 * y) * (Phi3 * y) ^ lambda * lambda - ( (Phi3 * y) ^ lambda - 1)) / lambda ^ 2
            lnL.lambda <- vector(mode="numeric", length = length(y))
            lnL.lambda[y == 0] <- (( (phi2 - Phi12$b * Phi23$f - Phi12$f * Phi23$a) /
                           (Phi2 - Phi12$f * Phi23$f) / sigma - mills((bX2 + 1 / lambda) / sigma) / sigma
                           ) * (- 1 / lambda ^ 2))[y == 0]
            lnL.lambda[y != 0] <- (( -resid / sigma ^ 2 +  mills(z1$f) * rho1 / sigma / sqrt(1 - rho1 ^ 2) +
                           mills(z3$f) * rho3 / sigma / sqrt(1 - rho3 ^ 2) ) * Tylb + lnJlb )[y != 0]
            gradi <- cbind(gradi, lnL.lambda)
        }

        if (dist == "ihs"){
            Tylb <- (y * Phi3) / lambda / sqrt(1 + (lambda * y * Phi3) ^ 2) - Ty / lambda
            lnL.lambda <- vector(mode = "numeric", length = length(y))
            lnL.lambda[y != 0] <- (( -resid / sigma ^ 2 +  mills(z1$f) * rho1 / sigma / sqrt(1 - rho1 ^ 2) +
                           mills(z3$f) * rho3 / sigma / sqrt(1 - rho3 ^ 2) ) * Tylb + lnJlb)[y != 0]
            gradi <- cbind(gradi, lnL.lambda)
        }
        attr(lnL, "gradient") <- gradi
    }
    if (fitted){
        P0 <- exp(lnL.null)
        if (dist != "ln"){
            if (h3) Psi23 <- rho3 * phi3 * pnorm( (bX2 / sigma - rho3 * bX3) / sqrt(1 - rho3 ^ 2)) +
                phi2 * pnorm( (bX3 - rho3 * bX2 / sigma) / sqrt(1 - rho3 ^ 2))
            else Psi23 <- phi2
            if (h1) Psi21 <- rho1 * phi1 * pnorm( (bX2 / sigma - rho1 * bX1) / sqrt(1 - rho1 ^ 2)) +
                phi2 * pnorm( (bX1 - rho1 * bX2 / sigma) / sqrt(1 - rho1 ^ 2))
            else Psi21 <- phi2
            Econd <- bX2 / Phi3 + sigma * (Psi23 * Psi21 * Phi2) / (Phi12$f * Phi23$f * Phi3 * phi2)
        }
        else{
            if (h3) Psi23 <- pnorm(bX3 + sigma * rho3) else Psi23 <- 1
            if (h1) Psi21 <- pnorm(bX1 + sigma * rho1) else Psi21 <- 1
            Econd <- exp(bX2 + sigma ^ 2 / 2) * Psi21 * Psi23 / (Phi1 * Phi3 ^ 2)
        }
        attr(lnL, "fitted") <- cbind("P(y=0)" = P0, "E(y|y>0)" = Econd)
    }
    lnL
}

# Compute the estimation of hurdle models in the cases where it can be
# done using two independent estimations (a binomial logit model and a
# normal/log-normal/truncated linear model). This is relevant for
# uncorrelated models with selection

fit.simple.mhurdle <- function(X1, X2, y, dist = NULL){
  probit <- glm(y != 0 ~ X1 - 1, family = binomial(link = "probit"))
  lin <- switch(dist,
                "ln" = lm(log(y) ~ X2 - 1, subset = y != 0),
                "n" = lm(y ~ X2 - 1, subset = y != 0),
                "tn" = truncreg(y ~ X2 - 1, subset = y != 0)
                )
  df <- df.residual(lin)
  np <- sum(y != 0)
  K1 <- ncol(X1)
  beta1 <- coef(probit)
  bX1 <- as.numeric(crossprod(beta1, t(X1)))
  K2 <- ncol(X2)
  if (dist == "tn"){
    sigma <- coef(lin)[ncol(X2)+1]
    beta2 <- coef(lin)[-(ncol(X2)+1)]
  }
  else beta2 <- coef(lin)
  bX2 <- as.numeric(crossprod(beta2, t(X2)))
  L.null <- (y == 0) * log(1 - pnorm(bX1))
  if (dist == "ln"){
    logy <- rep(0, length(y))
    logy[y != 0] <- log(y[y != 0])
    resid <- (logy - bX2)
  }
  else resid <- y - bX2
  scr <- sum(resid[y != 0] ^ 2)

  if (dist != "tn") sigma <- sqrt(scr / np)
  
  mills1 <- mills(bX1)
  mills2 <- mills(bX2 / sigma) / pnorm(bX2 / sigma)
  mills1m <- mills(- bX1)
  
  L.pos <- switch(dist,
                  "ln" = (y != 0) * (- logy + pnorm(bX1, log.p = TRUE) + dnorm(resid / sigma, log = TRUE) - log(sigma)),
                  "n" = (y != 0) * (pnorm(bX1, log.p = TRUE) + dnorm(resid / sigma, log = TRUE) - log(sigma)),
                  "tn" = (y != 0) * (pnorm(bX1, log.p = TRUE)
                           + dnorm(resid / sigma, log = TRUE) - log(sigma)
                           - pnorm(bX2 / sigma, log.p = TRUE))
                  )
  gbX1 <- switch(dist,
                "ln" = (y == 0) * (- mills1m) + (y != 0) * mills1,
                "n" = (y == 0) * (- mills1m) + (y != 0) * mills1,
                "tn" = (y == 0) * (- mills1m) + (y != 0) * mills1
                )

  gbX2 <- switch(dist,
                "ln" = (y != 0) * (resid / sigma^2),
                "n" = (y != 0) * (resid / sigma^2),
                "tn" = (y != 0) * (resid / sigma^2 - 1 / sigma * mills2)
                ) 

  gsigma <- switch(dist,
                "ln" = (y != 0) * (resid^2 / sigma^3 - 1/ sigma),
                "n" = (y != 0) * (resid^2 / sigma^3 - 1/ sigma),
                "tn" = (y != 0) * (resid^2 / sigma^3 - 1/ sigma + bX2 / sigma^2 * mills2)
                )

  gradi <- cbind(gbX1 * X1, gbX2 * X2, as.numeric(gsigma))
  dss <- - 3 * scr / sigma ^ 4 + sum(y != 0) / sigma ^ 2

  if (dist == "tn"){
    vcov <- bdiag(vcov(probit), vcov(lin))
    coef <- c(coef(probit), coef(lin))
  }
  else{    
    vcov <- bdiag(vcov(probit), vcov(lin) / np * df, - 1 / dss)
    coef <- c(coef(probit), coef(lin), sigma)
  }
  fit <- cbind(zero = bX1, pos = bX2)

  other.coef <- c("sd")

  coef.names <- list(h1    = colnames(X1),
                     h2    = colnames(X2),
                     sd    = other.coef)

  fitted <- attr(mhurdle.lnl(coef, X1 = X1, X2 = X2, X3 = NULL, X4 = NULL, y = y,
                             gradient = FALSE, fitted = TRUE,
                             dist = dist, corr = NULL), "fitted")
  logLik <- structure(sum(L.null + L.pos), df = length(coef), nobs = length(y), class = "logLik")
  result <- list(coefficients = coef, 
                 vcov = vcov,
                 fitted.values = fitted,
                 logLik = logLik,
                 gradient = gradi,
                 model = NULL,
                 formula = NULL,
                 coef.names = coef.names,
                 call = NULL
                 )
  
  class(result) <- c("mhurdle","maxLik")
  result
}

# Compute the "naive" model, i.e. a model with no explanatory
# variables.  Full version with correlation ; not used because
# identification problems

if (FALSE){
lnl.naive <- function(param, dist = c("ln", "tn", "n"), moments,
                     h1 = TRUE, h3 = FALSE,
                     which = c("all", "zero", "positive")){
  dist <- match.arg(dist)
  which <- match.arg(which)
  n <- moments[1]
  ym <- moments[2]
  s2 <- moments[3]
  if (h1){
    alpha1 <- param[1]
    alpha2 <- param[2]
    param <- param[-c(1,2)]
  }
  else{
    alpha2 <- param[1]
    param <- param[-1]
  }
  if (h3){
    alpha3 <- param[1]
    param <- param[-1]
  }
  sigma <- param[1]
  
  if (length(param) == 2) rho <- param[2] else rho <- 0
  if (rho < - 1) rho <- - 0.99
  if (rho > 1) rho <- 0.99
  rho1 <- rho
  if (h1){
    Phi1 <- pnorm(alpha1)
    phi1 <- dnorm(alpha1)
  }
  else Phi1 <- 1
  if (h3){
    Phi3 <- pnorm(alpha3)
    phi3 <- dnorm(alpha3)
  }
  else Phi3 <- 1
  Phi2 <- pnorm(alpha2/sigma)
  phi2 <- dnorm(alpha2/sigma)
  scr <- ifelse(dist == "ln",
                s2 + (ym + log(Phi3) - alpha2)^2,
                Phi3^2*(s2 + (ym - alpha2/Phi3)^2)
                )
  if (!rho){
    Pbiv <- Phi1 * Phi2
    Phi1bis <- Phi1
    s2term <- 0
  }
  else{
    Pbiv <- mypbivnorm(alpha1, alpha2/sigma, rho)$f
    zo <- switch(dist,
                 "ln" = ym + log(Phi3) - alpha2 ,
                 Phi3 * (ym - alpha2/Phi3)
                 )
    millso <- dnorm(zo)/pnorm(zo)
    Phi1bis <-(alpha1 + rho/sigma * zo)/sqrt(1 - rho^2)
    s2term <- 0.5 * s2 * (rho / (sigma * sqrt(1 - rho^2) * Phi3^(dist != "ln")))^2 *
      millso * (zo + millso)
  }
  P0 <- switch(dist,
               "ln" = 1 - Phi1 * Phi3,
               "tn" = 1 - Pbiv/Phi2 * Phi3,
               "n" = 1 - Pbiv * Phi3
               )
  lnPos <-
    -log(sigma) - 0.5 * log(2*pi) -
      scr/(2*sigma^2) +
        log(Phi3) +
          (log(Phi1bis)+s2term) * h1 -
            ym * (dist == "ln") +
              log(Phi3) * (dist != "ln") -
                log(Phi2) * (dist == "tn")
  switch(which,
         "all" = n * log(P0) + (1 - n) * lnPos,
         "zero" = n * log(P0),
         "positive" = (1 - n) * lnPos)
}
}
# Version without correlation

if (TRUE){
    lnl.naive <- function(param, dist = c("ln", "tn", "n"), moments,
                          h1 = TRUE, h3 = FALSE){
        dist <- match.arg(dist)
        n <- moments[1]
        ym <- moments[2]
        s2 <- moments[3]
        if (h1){
            alpha1 <- param[1]
            alpha2 <- param[2]
            param <- param[-c(1,2)]
        }
        else{
            alpha2 <- param[1]
            param <- param[-1]
        }
        if (h3){
            alpha3 <- param[1]
            param <- param[-1]
        }
        sigma <- param[1]
        
        if (h1){
            Phi1 <- pnorm(alpha1)
            phi1 <- dnorm(alpha1)
        }
        else Phi1 <- 1
        if (h3){
            Phi3 <- pnorm(alpha3)
            phi3 <- dnorm(alpha3)
        }
        else Phi3 <- 1
        Phi2 <- pnorm(alpha2/sigma)
        phi2 <- dnorm(alpha2/sigma)
        scr <- ifelse(dist == "ln",
                      s2 + (ym + log(Phi3) - alpha2)^2,
                      Phi3^2*(s2 + (ym - alpha2/Phi3)^2)
                      )
        Pbiv <- Phi1 * Phi2
        Phi1bis <- Phi1
        s2term <- 0
        P0 <- switch(dist,
                     "ln" = 1 - Phi1 * Phi3,
                     "tn" = 1 - Pbiv/Phi2 * Phi3,
                     "n" = 1 - Pbiv * Phi3
                     )
        lnPos <-
            -log(sigma) - 0.5 * log(2*pi) -
                scr/(2*sigma^2) +
                    log(Phi3) +
                        (log(Phi1bis)+s2term) * h1 -
                            ym * (dist == "ln") +
                                log(Phi3) * (dist != "ln") -
                                    log(Phi2) * (dist == "tn")
        
        n * log(P0) + (1 - n) * lnPos
    }
}
