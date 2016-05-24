#' Estimate parameters from a specified model using pseudo-score
#'
#' @param psdesign An object of class psdesign
#' @param start Vector of starting values, if NULL, will come up with starting values
#' @param epsilon Convergence criteria
#' @param maxit Maximum number of iterations
#'
#' @export
#'
pseudo_score <- function(psdesign, start = NULL, epsilon = 1e-5, maxit = 50){

    if(psdesign$risk.model$model != "binary") stop("Only binary models supported for pseudo-score")

    if(is.null(start)){

      start <- rep(0, psdesign$nparam)

    }

    ## estimate P(Delta = 1 | Y, Z)
    aug <- psdesign$augdata
    del.lookup <- NULL
    for(i in c(0, 1)){ for(j in c(0, 1)){

      del.lookup <- rbind(del.lookup, data.frame(Y = i, Z = j, Pr = mean(!is.na(aug[aug$Y == i & aug$Z == j, "S.1"]))))

    }}

    del0 <- aug[is.na(aug$S.1), ]
    del1 <- aug[!is.na(aug$S.1), ]
    del1$glmweights <- 1
    if(is.null(psdesign$risk.model$args$model)){
      form <- Y ~ S.1 * Z
    } else {
      form <- eval(psdesign$risk.model$args$model)
    }

    ## check that all levels of BIP in del0 have a value in del1, otherwise imputations won't work

    bip0 <- unique(del0$BIP)
    bipmatch <- sapply(bip0, function(x) x %in% del1$BIP)
    if(!all(bipmatch)) stop("BIP has too many levels, is it categorical?")

    impute <- lapply(1:nrow(del0), function(i){

      df <- del0[i, ]
      ws <- del1[del1$BIP == df$BIP, "S.1"]
      df1 <- df[rep(1, length(ws)), ]
      df1$S.1 <- ws
      df1

    })

    pz1 <- mean(aug$Z == 1)
    pz0 <- mean(aug$Z == 0)

    if(is.null(psdesign$risk.model$args$risk)){
      link <- "logit"
    } else {
      link <- strsplit(as.character(psdesign$risk.model$args$risk), ".", fixed = TRUE)[[1]][2]
    }
    stopifnot(link %in% c("probit", "logit"))

    beta0 <- start
    niter <- 0

    repeat{

      weighteddf <- lapply(impute, function(df){

        risk <- psdesign$risk.function(data = df, beta = beta0)
        df1 <- df0 <- df
        df1$Z <- 1
        df0$Z <- 0

        prDelt <- (del.lookup[1, "Pr"] * (1 - psdesign$risk.function(data = df0, beta = beta0)) * pz0 +
                     del.lookup[2, "Pr"] * (1 - psdesign$risk.function(data = df1, beta = beta0)) * pz1 +
                     del.lookup[3, "Pr"] * (psdesign$risk.function(data = df0, beta = beta0)) * pz0 +
                     del.lookup[4, "Pr"] * (psdesign$risk.function(data = df1, beta = beta0)) * pz1)

        df$glmweights <- (risk / prDelt) / sum(risk / prDelt)
        df

      })

      glmdf <- do.call(rbind, weighteddf)
      glmdf <- rbind(glmdf, del1)

      fit <- glm(form, data = glmdf, weights = glmdf$glmweights, family = quasibinomial(link = link))

      beta1 <- fit$coefficients
      diffbeta <- sum(abs(beta0 - beta1))
      niter <- niter + 1

      beta0 <- beta1
      if(diffbeta < epsilon | niter > maxit) break


    }

    list(par = beta1, value = diffbeta, counts = niter, convergence = ifelse(niter > maxit, 11, 0))

}