WYht <-
function(X, f, theta, cmat, conf.level, alternative, R, args)
  {
    bargs <- args
    dat <- estThetaRow(X = X, f = f, theta = theta)
    linmod <- lm(dat$inds ~ dat$fac -1)
    res <- residuals(linmod)
    resdat <- data.frame(x = res, y = dat$fac)
    names(resdat) <- c("res", "fac")
    coeff <- coefficients(linmod)

    estres <- function(X, f, coeff, cmat)
      {
        eres <- coeff
        ssq <- tapply(X, f, FUN = function(x){sum((x - mean(x))^2)})
        dfd <- tapply(X, f, FUN = function(x){length(x) - 1})
        evar <- sum(ssq) / sum(dfd)
        ni <- tapply(X, f, FUN = length)
        estC <- cmat %*% eres
        varC <- (cmat^2) %*% (evar / ni)
        return(list(
                    estC = estC,
                    varC = varC
                    )
               )

      }

    EST <- estres(X = dat$inds, f = dat$fac, coeff = coeff, cmat = cmat)

    teststat.org <- EST$estC / sqrt(EST$varC)

    ResTeststat <- function(X, i, f, cmat)
      {
        RNEW <- X[i]
        eresB <- tapply(X = RNEW, INDEX = f, FUN = mean)
        ssqB <- tapply(X = RNEW, INDEX = f, FUN = function(x){sum((x - mean(x))^2)})
        dfdB <- tapply(X = RNEW, INDEX = f, FUN = function(x){length(x) - 1})
        evarB <- sum(ssqB) / sum(dfdB)
        niB <- tapply(X = RNEW, INDEX = f, FUN = length)
        estCB <- cmat %*% eresB
        varCB <- (cmat^2) %*% (evarB / niB)
        teststatB <- estCB / sqrt(varCB)
        return(teststatB = teststatB)
      }

    bargs$data <- resdat$res
    bargs$R <- R
    bargs$statistic <- ResTeststat
    bargs$f <- resdat$fac
    bargs$cmat <- cmat
    if(is.null(bargs$sim))
      {
        bargs$sim <- "ordinary"
      }
    if(is.null(bargs$stype))
      {
        bargs$stype <- "i"
      }

    bootout <- do.call("boot", bargs)

    matraw <- matrix( c( teststat.org, bootout$t ), byrow = TRUE, ncol = ncol( bootout$t ), dimnames = NULL)

    alpha <- 1 - conf.level
    switch(alternative,
           "two.sided" =
           {
             maxabsT <- apply(X = bootout$t, MARGIN = 1, FUN = function(x){max(abs(x))})
             count <- sapply( lapply( X = teststat.org, FUN = function( x ){
               maxabsT >= abs( x )
             }), FUN = sum )

             countraw <- apply( apply( X = matraw, MARGIN = 2, FUN = function( x ){
               abs( x[2:length( x )] ) >= abs( x[1] )
             }), MARGIN = 2, FUN = sum)

             pval <- count / R
             pvalraw <- countraw / R

             quant <- quantile(maxabsT, probs = 1-alpha, na.rm = TRUE)
             LOWER <- EST$estC - quant * sqrt(EST$varC)
             UPPER <- EST$estC + quant * sqrt(EST$varC)
           },
           "less" =
           {
             maxT <- apply(X = bootout$t, MARGIN = 1, FUN = max)

             count <- sapply( lapply( X = teststat.org, FUN = function( x ){
               maxT >= x
             }), FUN = sum )

             countraw <- apply( apply( X = matraw, MARGIN = 2, FUN = function( x ){
               x[2:length( x )] >=  x[1]
             }), MARGIN = 2, FUN = sum)

             pval <- count / R
             pvalraw <- countraw / R

             quant <- quantile(maxT, probs = 1-alpha, na.rm = TRUE)
             LOWER <- NA
             UPPER <- EST$estC + quant * sqrt(EST$varC)
           },
           "greater" =
           {
             minT <- apply(X = bootout$t, MARGIN = 1, FUN = min)

             count <- sapply( lapply( X = teststat.org, FUN = function( x ){
               minT <= x
             }), FUN = sum )

             countraw <- apply( apply( X = matraw, MARGIN = 2, FUN = function( x ){
               x[2:length( x )] <=  x[1]
             }), MARGIN = 2, FUN = sum)

             pval <- count / R
             pvalraw <- countraw / R

             quant <- quantile(minT, probs = alpha, na.rm = TRUE)
             LOWER <- EST$estC + quant * sqrt(EST$varC)
             UPPER <- NA
           })

    conf.int <- cbind(EST$estC, LOWER, UPPER)
    colnames(conf.int) <- cbind("estimate", "lower", "upper")

    p.value <- matrix(c( pval, pvalraw ), ncol = 2, dimnames = list(dimnames(cmat)[[1]], c("adj. p", "raw p")))

    return(list(conf.int = round(conf.int, 3), p.value = round(p.value, 3), conf.level = conf.level, alternative = alternative))
  }
