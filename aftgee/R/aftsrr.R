Ln <- function(beta, other) {
  ##  other <- list(y, x, delta, clsize, sigma)
  Y <- other[[1]]
  X <- other[[2]]
  delta <- other[[3]]
  clsize <- other[[4]]
  sigma <- other[[5]]
  weights <- other[[6]]
  Z <- other[[7]]
  gw <- other[[8]]
  p  <- ncol(X)
  N <- nrow(X)
  n <- length(clsize)
  ln <- double(1)
  .C("lfun", as.double(beta), as.double(Y), as.double(X), as.double(delta),
     as.integer(clsize), as.double(sigma), as.integer(n), as.integer(p),
     as.integer(N), as.double(weights), as.double(gw),
     as.double(Z), out = as.double(ln), PACKAGE = "aftgee")$out
}


abargehanfun <- function(beta, Y, X, delta, clsize, sigma, weights, gw = rep(1, nrow(X))) {
    p <- ncol(X)
    N <- nrow(X)
    n <- length(clsize)
    a <- vector("double", p * p)
    matrix(.C("abargehanfun", as.double(beta), as.double(Y), as.double(X), as.double(delta),
              as.integer(clsize), as.double(sigma), as.integer(n), as.integer(p),
              as.integer(N), as.double(weights), as.double(gw),
              out = as.double(a), PACKAGE = "aftgee")$out, nrow = p)
}


abarlogfun <- function(beta, Y, X, delta, clsize, sigma, weights, pw = rep(1, nrow(X))) {
    p <- ncol(X)
    N <- nrow(X)
    n <- length(clsize)
    a <- vector("double", p * p)
    matrix(.C("abarlogfun", as.double(beta), as.double(Y), as.double(X), as.double(delta),
              as.integer(clsize), as.double(pw), as.double(sigma), as.integer(n),
              as.integer(p), as.integer(N), as.double(weights),
              out = as.double(a), PACKAGE = "aftgee")$out, nrow = p)
}

abarpwfun <- function(beta, Y, X, delta, clsize, sigma, weights, pw) {
    p <- ncol(X)
    N <- nrow(X)
    n <- length(clsize)
    a <- vector("double", p * p)
    pt1 <- matrix(.C("abarpwfun", as.double(beta), as.double(Y), as.double(X),
                     as.double(delta), as.integer(clsize), as.double(pw$fhat),
                     as.double(sigma), as.integer(n), as.integer(p), as.integer(N),
                     as.double(weights), out = as.double(a), PACKAGE = "aftgee")$out, p)
    ## pt1 <- uilogFun(beta, Y, X, delta, clsize, sigma, n, Z = rep(1, nrow(X)), weights, smooth = TRUE, constant = 0, s = 0, pw = pw$fhat)
    pt2 <- abarlogfun(beta, Y, X, delta, clsize, sigma, weights, pw$Shat)
    ## rep(1, p) %o% pt1 + pt2
    ## diag(pt1) + pt2
    pt1 + pt2
}

omegaFun <- function(beta, Y, X, delta, clsize, weights) {
  p <- ncol(X)
  N <- nrow(X)
  n <- length(clsize)
  omega <- vector("double", p * p)
  matrix(.C("omegafun", as.double(beta), as.double(Y), as.double(X),
            as.double(delta), as.integer(clsize), as.integer(n), as.integer(p),
            as.integer(N), as.double(weights),
            out = as.double(omega), PACKAGE = "aftgee")$out, nrow = p)
}

getRankName <- function(rankWeights, method) {
    rktemp <- rankWeights
    if (method == "nonsm") {
        if (rktemp == "logrank") {rankWeights <- "nslogrank"}
        if (rktemp == "PW") {rankWeights <- "nsPW"}
        if (rktemp == "GP") {rankWeights <- "nsGP"}
    }
    if (method == "monosm") {
        if (rktemp == "PW") {
            rankWeights <- "mPW"
            ## variance[variance == "ZLMB"] <- "sZLMB"
        }
        if (rktemp == "GP") {
            rankWeights <- "mGP"
            ## variance[variance == "ZLMB"] <- "sZLMB"
        }
        if (rktemp == "logrank") {rankWeights <- "mlogrank"}
    }
    if (method == "sm") {
        if (rktemp == "gehan") {rankWeights <- "gehan"}
        if (rktemp == "logrank") {rankWeights <- "logrank"}
        ## variance[variance == "ZLMB"] <- "sZLMB"
    }
    out <- list(rankWeights = rankWeights)
}

getSuv <- function(Y, X, beta, N, delta, weights) {
    en <- Y - X %*% beta
    Shat <- fhat <- NULL
    dummy <- 1:N
    ord <- order(en)
    ei <- en[ord]
    weightsi <- weights[ord]
    deltai <- delta[ord]
    repeats <- table(ei)
    ## Shati <- survfit(Surv(ei, deltai) ~ 1, weights = weightsi)$surv
    Shati <- exp(-1 * basehaz(coxph(Surv(ei, deltai)~1, weights = weightsi))$hazard)
    ## fhati <- diff(c(Shati[1], Shati, Shati[length(Shati)]), 2) / diff(c(ei[1], ei, ei[length(ei)]), 2)
    Shati <- rep(Shati, repeats)
    Shatlast <- rev(Shati)[1]
    ## Shatlast <- Shati[1]
    ## fhati <- rep(fhati, repeats)
    Shat[dummy[ord]] <- Shati
    ## fhat[dummy[ord]] <- fhati
    ## fhat <- ifelse(fhat < -1, -1, fhat)
    ## mu <- mean(subset(fhat, fhat > -1))
    ## fhat <- ifelse(fhat < -1, mu, fhat)
    ## list(Shat = Shat, fhat = fhat, Shatlast = Shatlast)
    list(Shat = Shat, Shatlast = Shatlast)
}

getWi <- function(Y, X, beta, N, delta, weights) {
    en <- Y - X %*% beta
    shat <- fhat <- NULL
    dummy <- 1:N
    ord <- order(en)
    ei <- en[ord]
    weightsi <- weights[ord]
    deltai <- delta[ord]
    repeats <- table(ei)
    Shati <- survfit(Surv(ei, deltai) ~ 1, weights = weightsi)$surv
    Shati <- rep(Shati, repeats)
    ## shat <- getSuv(Y, X, beta, N, delta, weights)$Shat
    what <- NULL
    whati <- cumsum(Shati)
    shat[dummy[ord]] <- Shati
    what[dummy[ord]] <- whati
    list(what = what, shat = shat)
}
   
getGehan <- function(Y, X, beta, N, delta, clsize, sigma, weights, smooth = FALSE) {
    p <- ncol(X)
    N <- nrow(X)
    n <- length(clsize)
    a <- vector("double", N)
    if (smooth == TRUE) {
        out <- matrix(.C("getgehan", as.double(beta), as.double(Y), as.double(X), as.integer(clsize),
                         as.double(sigma), as.integer(n), as.integer(p), as.integer(N), as.double(weights),
                         out = as.double(a), PACKAGE = "aftgee")$out, ncol = 1)
    }
    if (smooth == FALSE) {
        out <- matrix(.C("getnsgehan", as.double(beta), as.double(Y), as.double(X), as.integer(clsize),
                         as.double(sigma), as.integer(n), as.integer(p), as.integer(N), as.double(weights),
                         out = as.double(a), PACKAGE = "aftgee")$out, ncol = 1) ##+ 1
    }

    out
}

getPw <- function(Y, X, beta, N, delta, weights, rankWeights) {
    if (rankWeights == "logrank") {
        pw <- rep(1, nrow(X))
    }
    if (rankWeights == "PW") {
        pw <- getSuv(Y, X, beta, N, delta, weights)$Shat 
    }
    if (rankWeights == "GP") {
        pw <- getSuv(Y, X, beta, N, delta, weights)$Shat ^ (1 / ncol(X))
    }
    if (rankWeights == "eGP") {
        pw <- getSuv(Y, X, beta, N, delta, weights)
        pw <- ((pw$Shat - pw$Shatlast) / (1 - pw$Shatlast))^ (1 / ncol(X))
    }
    pw
}

getGw <- function(Y, X, beta, N, delta, clsize, sigma, weights, rankWeights, pw = NULL) {
    de <- getGehan(Y, X, beta, N, delta, clsize, sigma, weights)
    if (is.numeric(pw)) {
        gw <- pw / de
    }
    if (!is.numeric(pw)) {
        if (rankWeights == "gehan") {
            gw <- rep(1, nrow(X))
        }
        if (rankWeights == "logrank") {
            gw <- 1 / de
        }
        if (rankWeights == "PW") {
            ne <- getSuv(Y, X, beta, N, delta, weights)$Shat
            gw <- ne / de
        }
        if (rankWeights == "GP") {
            ne <- getSuv(Y, X, beta, N, delta, weights)$Shat ^ (1 / ncol(X))
            gw <- ne / de
        }
    }
    gw
}

getSmoothSuv <- function(Y, X, beta, N, delta, weights) {
    en <- Y - X %*% beta
    ik <- rep(1:N, each=N)
    jl <- rep(1:N, N)
    Shat <- fhat <- NULL
    dummy <- 1:N
    ord <- order(en)
    rij <- sqrt(diag((X[ord,] - X[dummy,]) %*% t((X[ord,] - X[dummy,]))))
    ei <- en[ord]
    weightsi <- weights[ord]
    deltai <- delta[ord]
    repeats <- table(ei)
    ## Shati <- survfit(Surv(ei, deltai) ~ 1, weights = weightsi)$surv
    ## assume no ties
    di <- rev(cumsum(rev(rep(1, N))))
    Z <- pnorm((en - ei) / rij)
    Z <- ifelse(is.na(Z) == T, 0, Z)
    hazi <- cumsum(deltai * Z / di)
    Shati <- exp(-1 * hazi)
    fhati <- diff(c(Shati[1], Shati, Shati[length(Shati)]), 2) / diff(c(ei[1], ei, ei[length(ei)]), 2)
    Shati <- rep(Shati, repeats)
    fhati <- rep(fhati, repeats)
    Shat[dummy[ord]] <- Shati
    fhat[dummy[ord]] <- fhati
    ## fhat <- ifelse(fhat < -1, -1, fhat)
    mu <- mean(subset(fhat, fhat > -1))
    fhat <- ifelse(fhat < -1, mu, fhat)
    list(Shat = Shat, fhat = fhat)
}

uiFun <- function(beta, Y, X, delta, clsize, sigma, n, Z, weights, smooth = TRUE, gw, constant = 0) {
  N <- nrow(X)
  p <- ncol(X)
  ans <- numeric(p)
  sn <- vector("double", p)
  ans <- .C("ufun", as.double(beta), as.double(Y), as.double(X), as.double(delta),
            as.integer(clsize), as.double(sigma), as.integer(n), as.integer(p),
            as.integer(N), as.double(Z), as.double(weights), as.double(gw),
            out = as.double(sn), PACKAGE = "aftgee")$out
  ans <- ans - constant
}

uilogFun <- function(beta, Y, X, delta, clsize, sigma, n, Z,
                     weights, constant = 0, pw = rep(1, nrow(X)), rankWeights, rkmethod) {
    N <- sum(clsize)
    p <- ncol(X)
    sn <- vector("double", p)
    ans <- numeric(p)
    n <- length(clsize)
    if (rkmethod != "nonsm") {
        ans <- .C("ulogfun", as.double(beta), as.double(Y), as.double(X),
                  as.double(delta), as.integer(clsize),
                  as.double(sigma), as.integer(n), as.integer(p),
                  as.integer(N), as.double(Z), as.double(weights),
                  as.double(pw), out = as.double(sn), PACKAGE = "aftgee")$out
    }
    if (rkmethod == "nonsm") {
        pw <- getPw(Y = Y, X = X, beta = beta, N = N, delta = delta,
                    weights = weights, rankWeights)
        ans <- .C("ulognsfun", as.double(beta), as.double(Y), as.double(X),
                  as.double(delta), as.integer(clsize),
                  as.double(sigma), as.integer(n), as.integer(p),
                  as.integer(N), as.double(Z), as.double(weights),
                  as.double(pw), out = as.double(sn), PACKAGE = "aftgee")$out
    }
    ans - constant
}

solveGehan <- function(beta, Y, X, delta, clsize, sigma, weights,
                       Z, gw = rep(1, nrow(X))) {
    tbeta <- system.time(temp <- nlm(Ln, p = beta,
                                     other = list(Y, X, delta, clsize, sigma, weights, Z, gw),
                                     fscale = 0.01))
    conv <- ifelse(temp$code == 1, 0, 1)
    out <- list(beta = temp$estimate, tbeta = tbeta, conv = conv)
}

solvePW <- function(beta, Y, X, delta, clsize, sigma, weights, Z,
                    rankWeights, method, constant, control = aftgee.control()){
    p <- ncol(X)
    N <- sum(clsize)
    if (method == "sm") {
        pw <- getPw(Y = Y, X = X, beta = beta, N = N, delta = delta,
                     weights = weights, rankWeights)
        btime <- system.time(temp <- BBsolve(beta, uilogFun, Y = Y, X = X, delta = delta,
                                             clsize = clsize, sigma = sigma, n = N, Z = Z,
                                             weights = weights, constant = 0, pw = pw,
                                             rankWeights = rankWeights, rkmethod = method,
                                             quiet = TRUE))
        out <- list(beta = temp$par, btime = btime[3])
    }
    if (method == "monosm") {
        gw <- getGw(Y = Y, X = X, beta = beta, N = N, delta = delta,
                     clsize = clsize, sigma = sigma, weights = weights,
                     rankWeights = rankWeights)
        btime <- system.time(temp <- solveGehan(beta, Y, X, delta, clsize, sigma, weights, Z, gw))
        out <- list(beta = temp$beta, btime = btime[3])
    }
    out$gw <- ifelse(method == "monosm", gw, rep(1, nrow(X)))
    out   
}

viEmp <- function(beta, Y, delta, X, id, weights = rep(1, nrow(X)), B = 500,
                  mb = TRUE, zbeta = FALSE, smooth = TRUE,
                  rankWeights = "gehan", gw = gw,
                  sigma = diag(ncol(X)), gpweight = 1){
  p <- ncol(X)
  clsize <- unlist(lapply(split(id, id), length))
  n <- sum(unlist(lapply(split(weights, id), unique)))
  N <- sum(weights)
  UnV <- zmat <- matrix(0, ncol = B, nrow = p)
  for (i in 1:B) {
    if ( mb == TRUE) {
      Z <- rep(rexp(length(clsize)), clsize)
    }
    if ( mb != TRUE) {
      Z <- rep(1, N)
    }
    if (zbeta == TRUE) {
      zb <- rnorm(p)
      newbeta <- beta + n ^ (-0.5) * zb
      zmat[,i] <- zb
    }
    if (zbeta != TRUE) {
      newbeta <- beta
    }
    sn <- vector("double", p)
    ## if (rankWeights %in% c("nsPW", "nsGP")) {
    ##     smooth <- FALSE
    ## }
    gpweight <- ifelse(rankWeights %in% c("nsGP", "GP"), gpweight, 1)
    if (smooth == TRUE) {
        if (rankWeights == "gehan") {
            UnV[,i] <- as.vector(.C("ufun", as.double(newbeta), as.double(Y),
                                    as.double(X), as.double(delta),
                                    as.integer(clsize), as.double(sigma),
                                    as.integer(length(clsize)), as.integer(p),
                                    as.integer(sum(clsize)), as.double(Z),
                                    as.double(weights), as.double(gw),
                                    out = as.double(sn), PACKAGE = "aftgee")$out) # / n
        }
        if (rankWeights %in% c("nslogrank", "logrank")) {
            UnV[,i] <- uilogFun(beta = newbeta, Y = Y, X = X, delta = delta,
                                clsize = clsize, sigma = sigma, n = length(clsize),
                                Z = Z, weights = weights, constant = 0,
                                pw = rep(1, sum(clsize)),
                                rankWeights = "logrank", rkmethod = "sm") # / n
        }
        if (rankWeights %in% c("mPW", "mGP", "mlogrank", "userdefined")) {
            if (rankWeights == "mGP") {
                gw <- getGw(Y = Y, X = X, beta = beta, N = n, delta = delta,
                             clsize = clsize, sigma = sigma, weights = weights,
                             rankWeights = "GP")
            }
            if (rankWeights == "mlogrank") {
                gw <- getGw(Y = Y, X = X, beta = beta, N = n, delta = delta,
                             clsize = clsize, sigma = sigma, weights = weights,
                             rankWeights = "logrank")
            }
            if (rankWeights == "mPW") {
                gw <- getGw(Y = Y, X = X, beta = beta, N = n, delta = delta,
                             clsize = clsize, sigma = sigma, weights = weights,
                             rankWeights = "PW")
            }
            UnV[,i] <- as.vector(.C("ufun", as.double(newbeta), as.double(Y),
                                    as.double(X), as.double(delta),
                                    as.integer(clsize), as.double(sigma),
                                    as.integer(length(clsize)), as.integer(p),
                                    as.integer(sum(clsize)), as.double(Z),
                                    as.double(weights), as.double(gw),
                                    out = as.double(sn), PACKAGE = "aftgee")$out) # / n
        }
        if (rankWeights %in% c("PW", "GP", "eGP")) {
            pw <- getPw(Y = Y, X = X, beta = newbeta, N = nrow(X),
                         delta = delta, weights = weights, rankWeights)
            UnV[,i] <- uilogFun(beta = newbeta, Y = Y, X = X, delta = delta,
                                clsize = clsize, sigma = sigma,
                                n = nrow(X), Z = Z, pw =pw,
                                weights = weights, constant = 0, rkmethod = "sm",
                                rankWeights = rankWeights) # / n
        }
    }
    if (smooth == FALSE) {
        if (rankWeights == "gehan") {
            UnV[,i] <- as.vector(.C("unsfun", as.double(newbeta), as.double(Y),
                                    as.double(X), as.double(delta),
                                    as.integer(clsize), as.double(sigma),
                                    as.integer(length(clsize)), as.integer(p),
                                    as.integer(sum(clsize)), as.double(Z),
                                    as.double(weights), as.double(gw),
                                    out = as.double(sn), PACKAGE = "aftgee")$out) # / n ## n and N?
        }
        if (rankWeights %in% c("logrank", "nslogrank")) {
            UnV[,i] <- uilogFun(beta = newbeta, Y = Y, X = X, delta = delta,
                                clsize = clsize, sigma = sigma,
                                n = length(clsize), Z = Z, weights = weights, constant = 0,
                                pw = rep(1, sum(clsize)),
                                rankWeights = "logrank", rkmethod = "nonsm") # / n
        }
        if (rankWeights %in% c("nsPW", "nsGP", "mPW", "mGP", "PW", "GP", "eGP")) {
            UnV[,i] <- uilogFun(beta = newbeta, Y = Y, X = X, delta = delta,
                                clsize = clsize, sigma = sigma,
                                n = length(clsize), Z = Z, weights = weights,
                                constant = 0, rkmethod = "nonsm",
                                pw = rep(1, sum(clsize)),
                                rankWeights = rankWeights) # / n
        }
    }
}
  vi <- var(t(UnV))
  list(vi = vi, zmat = zmat, UnV = UnV)
}


getSi <- function(beta, Y, delta, X, id, weights = rep(1, nrow(X)),
                   rankWeights = "gehan") {
    p <- ncol(X)
    clsize <- unlist(lapply(split(id, id), length))
    N <- sum(clsize)
    en <- Y - X %*% beta
    ik <- rep(1:N, each=N)
    jl <- rep(1:N, N)
    Shat <- NULL
    dummy <- 1:N
    ord <- order(en)
    ei <- en[ord]
    weightsi <- weights[ord]
    deltai <- delta[ord]
    repeats <- table(ei)
    Shati <- survfit(Surv(ei, deltai) ~ 1, weights = weightsi)$surv
    Shati <- rep(Shati, repeats)
    Shat[dummy[ord]] <- Shati
    xdif <- X[ik,] - X[jl,]
    edif <- en[ik] - en[jl]
    ind <- ifelse(edif <= 0, 1, 0)
    minEn <- ifelse(edif <= 0, ik, jl)
    si <- s <- NULL
    if (rankWeights == "gehan") {
        s <- weights[jl] * delta[ik] * ind * xdif + weights[jl] * log(Shat[minEn]) * xdif
        if (length(which(s == Inf)) > 0) {
            s <- ifelse(s == Inf, 0, s)
        }
        if (length(which(s == -Inf)) > 0) {
            s <- ifelse(s == -Inf, 0, s)
        }
        if (sum(is.na(s) > 0)) {
            s <- ifelse(is.na(s) == TRUE, 0, s)
        }
        s <- rowsum(s, ik)
        }
    if (rankWeights == "logrank") {
        haz <- -1 * log(Shat)
        haz <- ifelse(haz == Inf, 0, haz)
        haz <- ifelse(haz == -Inf, 0, haz)
        haz <- ifelse(is.na(haz) == TRUE, 0, haz)
        gamma1 <- rowsum(ind * X[jl,] * weights[jl], ik)
        gamma0 <- as.numeric(rowsum(ind * weights[jl], ik))
        si1i <-  si1 <- s2 <- matrix(0, nrow = N, ncol = p)
        si1 <- (X[1:N,] - gamma1 / gamma0)
        si1i <- si1[ord,]
        si1dif <- rbind(rep(0,p), diff(si1i))
        s2s <- apply(si1dif * haz[ord], 2, cumsum)
        s2[dummy[ord],] <- s2s
        s <- delta[1:N] * si1 - s2
    }
    s
}

getVi <- function(s, id, delta, weights, n) {
    clweights <- as.numeric(unlist(lapply(split(weights, id), unique)))
    s1 <- rowsum(s, group = id)
    si1 <- lapply(split(s1, 1:nrow(s1)), function(x) x %o% x)
    s11 <- mapply("*", si1, clweights, SIMPLIFY = FALSE)
    v1 <- Reduce("+", s11) / n
    v2 <- apply(s1 * clweights , 2, sum) / n
    v2 <- v2 %o% v2
    if (length(unique(weights)) == 1) {
        p <-  unique(weights) - 1
        vi <- p * (v1 - v2)
    }
    if (length(unique(weights)) > 1) {
        cweights <- unique(weights)
        cweights <- cweights[cweights != 1 & cweights != 0]
        ## vi <- (cweights - 1)* (v1 - v2)
        vi <- v1 - v2
    }
    list(vi = vi, v1 = v1, v2 = v2)
}


viClo <- function(beta, Y, delta, X, id, weights = rep(1, nrow(X)), B = 500,
                  gw = gw, rankWeights = "gehan", stratify = TRUE) {
    s <- v1 <- vi <- v2i <- NULL
    n <- sum(unlist(lapply(split(weights, id), unique)))
    p <- ncol(X)
    clsize <- unlist(lapply(split(id, id), length))
    clid <- unlist(lapply(clsize, function(x) 1:x))
    stra <- match(weights, unique(weights))
    dim <- unique(clsize)
    s <- getSi(beta = beta, Y = Y, delta = delta, X = X, id = id,
                weights = weights, rankWeights = rankWeights)
    if (stratify) {
        v1 <- getVi(s, id, delta, weights, n)$v1
        vi <- v1
        v2 <- matrix(0, ncol = p, nrow = p)
        if (length(unique(stra)) > 1) {
            for (i in 1:length(unique(stra))) {
                ns <- sum(unlist(lapply(split(weights[stra == i], id[stra == i]), unique)))
                ## weights at cluster level
                v2i <- getVi(s = s[stra == i, ], id = id[stra == i], delta = delta[stra == i],
                              weights = weights[stra == i],
                              n = ns
                              )$vi
                strPr <- ns / n
                v2 <- v2 + v2i * strPr
            }
            vi <- v1 + v2
        }
    }
    if (!(stratify)) {
        v1 <- getVi(s, id, delta, weights, n)$v1 # / sum(clsize)
        vi <- v1
        if (length(unique(stra)) > 1) {
            v2 <- getVi(s, id, delta, weights * (1 - delta), n)$vi # / sum(clsize)
            vi <- v1 + v2
        }
    }
    list(vi = as.matrix(vi), s = s)
}


huangFun <- function(beta, Y, delta, X, id, weights = rep(1, nrow(X)), B = 500,
                     vClose = TRUE, rankWeights = "gehan", gw = gw,
                     sigma = diag(ncol(X)), stratify = TRUE) {
  p <- ncol(X)
  clsize <- unlist(lapply(split(id, id), length))
  n <- sum(unlist(lapply(split(weights, id), unique)))
  N <- sum(weights)
  ## n <- length(clsize)
  ## N <- sum(clsize)
  betaTemp <- NULL
  UnMatV <- NULL
  if (vClose == TRUE) {
      vi <- viClo(beta, Y, delta, X, id, weights, B, rankWeights, stratify = stratify)$vi * n
  }
  if (vClose != TRUE) {
      vi <- viEmp(beta, Y, delta, X, id, weights, B, mb = TRUE, zbeta = FALSE,
                  smooth = FALSE, rankWeights = rankWeights, gw = gw)$vi
  }
  qq <- chol(vi)
  qq <- t(qq)
  newBeta <- NULL
  Z <- rep(1, sum(clsize)) #
  newBeta <- matrix(0, ncol = p, nrow = p)
  for ( i in 1:p) {
      if (rankWeights == "gehan") {
          bb <- BBsolve(beta, uiFun, Y = Y, X = X, delta = delta, clsize = clsize,
                        sigma = sigma, n = length(clsize), Z = Z, gw = gw,
                        weights = weights, smooth = TRUE, constant = qq[,i] * n ^ -0.5,
                        quiet = TRUE)
          newBeta[, i] <- bb$par
      }
      if (rankWeights == "logrank") {
          bb <- BBsolve(beta, uilogFun, Y = Y, X = X, delta = delta, clsize = clsize,
                        sigma = sigma, n = length(clsize), Z = Z,
                        weights = weights, smooth = TRUE, constant = qq[,i] * n ^ -0.5,
                        gpweight = 0, pw = rep(1, sum(clsize)),
                        rankWeights = rankWeights, quiet = TRUE)
          newBeta[, i] <- bb$par
      }
      if (rankWeights %in% c("PW", "Prentice-Wilcoxon", "GP")) {
          gpweight <- ifelse(rankWeights %in% c("GP", "nsGP"), gpweight, 1)
          bb <- BBsolve(beta, uilogFun, Y = Y, X = X, delta = delta, clsize = clsize,
                        sigma = sigma, n = length(clsize), Z = Z,
                        weights = weights, smooth = TRUE, constant = qq[,i] * n ^ -0.5,
                        gpweight = gpweight, pw = rep(1, sum(clsize)),
                        rankWeights = rankWeights, quiet = TRUE)
          newBeta[, i] <- bb$par
      }
  }
  dd <- newBeta - beta
  covmat <- n * t(dd) %*% dd
  list(covmat = covmat)
}

zlFun <- function(beta, Y, X, delta, id, weights = rep(1, nrow(X)),
                  B = 500, vClose = FALSE, rankWeights = "gehan", method = "sm",
                  gw = gw, stratify = TRUE,
                  sigma = diag(ncol(X))) {
    gpweight <- ifelse(rankWeights == "GP", 1/ncol(X), 1)
    rankWeights <- getRankName(rankWeights, method)$rankWeights
    ## if (sum(variance %in% "ZLMB") > 0) {
    ##     variance <- ifelse(rankWeights %in% c("mGP", "mPW"), "sZLMB", variance)
    ## }
    ## if (sum(variance %in% "ZLCF") > 0) {
    ##     variance <- ifelse(rankWeights %in% c("mGP", "mPW"), "sZLCF", variance)
    ## }    
    p <- ncol(X)
    clsize <- unlist(lapply(split(id, id), length))
    n <- sum(unlist(lapply(split(weights, id), unique)))
    smooth <- ifelse(method != "nonsm", TRUE, FALSE)
    UnMat <- zmat <- ahat <- unTime <- NULL
    An.inv <- 1
    
    UnV <- viEmp(beta, Y, delta, X, id, weights, B, mb = FALSE, zbeta = TRUE, smooth = smooth,
                 gw = gw, rankWeights = rankWeights, sigma = sigma,
                 gpweight = gpweight)
    
    zmat <- UnV$zmat
    UnV <- UnV$UnV
    if (vClose == TRUE) {
        vi <- viClo(beta, Y, delta, X, id, weights, B, rankWeights, stratify = stratify)$vi * n
    }
    if (vClose != TRUE) {
        vi <- viEmp(beta, Y, delta, X, id, weights, B, mb = TRUE, zbeta = FALSE, smooth = smooth,
                    gw = gw, rankWeights = rankWeights, gpweight = gpweight)$vi ## smooth
    }
    An <- matrix(0, nrow = p, ncol = p)
    for (i in 1:p) {
        An[i,] <- lm(UnV[i,] ~ matrix(t(zmat), ncol = p) - 1)$coef
    }
    if (qr(An)$rank != p) {
        covmat <- ginv(An) %*% vi %*% t(ginv(An))
        An.msg <- "An is singular"
        An.inv <- 0
    }
    if (qr(An)$rank == p) {
        covmat <- solve(An) %*% vi %*% solve(An)
        An.msg <- "An is nonsingular"
        An.inv <- 1
  }
    covmat <- covmat  / n
    covmat <- matrix(as.numeric(covmat), p)
    list(covmat = covmat, vi = vi, An.msg = An.msg, An.inv = An.inv)
}


isFun <- function(beta, Y, delta, X, id, weights = rep(1, nrow(X)), sigma, B = 500,
                  vClose = FALSE, gw = rep(1, nrow(X)), rankWeights = "gehan",
                  omega = FALSE, stratify = TRUE) {
    p <- ncol(X)
    clsize <- unlist(lapply(split(id, id), length))
    n <- sum(unlist(lapply(split(weights, id), unique)))
    ## n <- length(clsize)
    UnMat <- zmat <- ahat <- NULL
    An.inv <- 1
    if (omega == TRUE && rankWeights == "gehan") {
        vi <- omegaFun(beta, Y, X, delta, clsize, weights) * n## / n ^ 2
    }
    if (omega != TRUE && vClose == TRUE && rankWeights == "gehan") {
        vi <- viClo(beta, Y, delta, X, id, weights, B, gw = gw, stratify = stratify)$vi * n
    }
    if (omega != TRUE && vClose == TRUE && rankWeights == "logrank") {
        vi <- viClo(beta, Y, delta, X, id, weights, B, "logrank", stratify = stratify)$vi * n
    }
    if (omega != TRUE && vClose != TRUE) {
        vi <- viEmp(beta, Y, delta, X, id, weights, B, mb = TRUE, zbeta = FALSE,
                    gw = gw, rankWeights = rankWeights)$vi
    }
    if (rankWeights %in% c("gehan", "mPW", "mGP", "mlogrank")) {
        An <- abargehanfun(beta, Y, X, delta, clsize, sigma, weights, gw)
    }
    if (rankWeights == "logrank") {
        An <- abarlogfun(beta, Y, X, delta, clsize, sigma, weights)
    }
    if (rankWeights %in% c("PW", "GP", "Prentice-Wilcoxon")) {
        pw <- getSmoothSuv(Y, X, beta, n, delta, weights)
        An <- abarpwfun(beta, Y, X, delta, clsize, sigma, weights, pw)
    }

    if (qr(An)$rank != p) {
        covmat <- ginv(An) %*% vi %*% ginv(An)
        An.msg <- "An is singular"
        An.inv <- 0
    }
    if (qr(An)$rank == p) {
        covmat <- solve(An) %*% vi %*% solve(An)
        An.msg <- "An is nonsingular"
        An.inv <- 1
    }
    covmat <- matrix(as.numeric(covmat), p)
    list(covmat = covmat, vi = vi, An.msg = An.msg, An.inv = An.inv)
}


aftsrr <- function(formula, data, subset, id = NULL, contrasts = NULL, 
                   strata = NULL, weights = NULL,
                   rankWeights = "gehan", method = "sm",
                   variance = "ISMB",
                   B = 100, SigmaInit = NULL,
                   control = aftgee.control()) {
    rankNames <- c("gehan", "logrank", "PW", "GP", "userdefined")
    rankMethod <- c("sm", "nonsm", "monosm")
    varianceNames <- c("MB", "ZLCF", "sZLMB", "ZLMB", "sHCF", "sHMB", "ISCF", "ISMB")
    if (sum(variance %in% varianceNames) != length(variance %in% varianceNames))
         stop("Invalid variance estimates", call. = FALSE) 
    if (!(method %in% rankMethod)) 
        stop("Invalid method type", call. = FALSE)
    scall <- match.call()
    mnames <- c("", "formula", "data", "weights", "subset", "na.nation", "id", "strata")
    cnames <- names(scall)
    cnames <- cnames[match(mnames, cnames, 0)]
    mcall <- scall[cnames]
    mcall[[1]] <- as.name("model.frame")
    m <- eval(mcall, parent.frame())
    id <- model.extract(m, id)
    y <- model.extract(m, "response")
    if (ncol(y) > 2) 
        stop("aftsrr only supports Surv object with right censoring", call. = FALSE)
    if (is.null(id)) 
        id <- 1:nrow(y)
    N <- NROW(y)
    mterms <- attr(m, "terms")
    x <- model.matrix(mterms, m, contrasts) 
    weights <- model.extract(m, weights)
    strata <- model.extract(m, strata)
    if (is.null(weights))
        weights <- rep(1, N)
    stratify <- TRUE
    ##    if (is.null(strata))
    ##        stratify <- FALSE
    xnames <- colnames(x)
    if ("(Intercept)" %in% colnames(x)) {
        x <- as.matrix(x[,-1])
        xnames <- xnames[-1]
    }
    if(is.null(SigmaInit)) {
        if (ncol(x) == 1) 
            SigmaInit <- 1
        else
        SigmaInit <- diag(ncol(x))
    }
    ## check rankWeights 
    if (!is.numeric(rankWeights)) {
        pw <- NULL
        if (!(rankWeights %in% rankNames)) 
            stop("Invalid rankWeights weights", call. = FALSE)
    }
    if (is.numeric(rankWeights)) {
        if (length(rankWeights) != nrow(y)) 
            stop("rankWeights value length does not match with dimension", call. = FALSE)
        pw <- rankWeights
        rankWeights = "userdefined"
    }
    out <- aftsrr.fit(Y = log(y[,1]), delta = y[,2], X = x, id = id,
                        weights = weights, SigmaInit = SigmaInit,
                        variance = variance, B = B, rankWeights = rankWeights,
                        method = method, pw = pw, stratify = stratify, control = control)
    out$call <- scall
    out$vari.name <- xnames
    names(out$beta) <- xnames
    out$y <- y
    out$x <- as.matrix(x)
    out$B <- B
    class(out) <- "aftsrr"
    return(out)
}


aftsrr.fit <- function(Y, delta, X, id, weights = rep(1, nrow(X)),
                       variance = "ISMB", B = 100, rankWeights = "gehan", method = "sm",
                       SigmaInit = diag(ncol(X)), pw = NULL, stratify = TRUE,
                       control = aftgee.control()) {
    X <- as.matrix(X)
    p <- ncol(X)
    binit = "lm"
    if (is.numeric(binit)) {
        b0 <- binit
        if (length(b0) != p) {
            stop ("Initial value length does not match with the numbers of covariates", call. = FALSE)
        }
    }
    if (!(is.numeric(binit))) {
        b0 <- rep(1,p)
        if (binit == "lm") {
            b0 <- lm(Y ~ X - 1)$coef
        }
    }
    btemp1 <- b0
    btemp2 <- NULL
    clsize <- unlist(lapply(split(id, id), length))
    order <- unlist(lapply(clsize, function(x) 1:x))
    N <- sum(clsize)
    n <- length(clsize)
    Z <- rep(1, N)
    vMB <- vZLCF <- vZLMB <- vsHCF <- vsHMB <- vISCF <- vISMB <- bstep <- NaN
    
    ## Point Estimation
    pe <- aftsrr.pe(beta = btemp1, Y = Y, X = X, delta = delta,
                    clsize = clsize, sigma = SigmaInit, weights = weights,
                    Z = Z, pw = pw, rankWeights = rankWeights, method = method)
    tbeta <- pe$tbeta
    btemp2 <- pe$btemp2
    bconv <- pe$bconv
    bstep <- pe$bstep
    ## Variance Estimation
    if (B == 0) {
        vMB <- NULL
        out <- list(beta = btemp2, covmat = NULL,
                    convergence = bconv, bstep = bstep,
                    var.meth = variance)
    }
    ZLMB.An.inv <- ZLCF.An.inv <- ISMB.An.inv <- ISCF.An.inv <- js.An.inv <- 1
    bhist <- NULL
    if (sum(variance %in% "MB") > 0) {
        if (B > 0) {
            vtemp <- matrix(0, nrow = B, ncol = p)
            for (i in 1:B) {
                Z <- rep(rexp(length(clsize)), clsize)
                vtemp[i,] <- aftsrr.pe(beta = btemp1, Y = Y, X = X, delta = delta,
                                       clsize = clsize, sigma = SigmaInit, weights = weights,
                                       Z = Z, pw = pw, rankWeights = rankWeights,
                                       method = method)$btemp2
            }
            bhist <- vtemp
            vMB <- var(vtemp)
        }
    }
    if (sum(variance %in% c("ISMB", "ISCF", "ZLMB", "ZLCF", "sZLMB", "sZLCF", "sHMB", "sHCF")) > 0) {
        gw <- getGw(Y = Y, X = X, beta = btemp2, N = nrow(X), delta = delta,
                     clsize = clsize, sigma = SigmaInit, weights = weights,
                     rankWeights = rankWeights)
    }
    
    if (sum(variance %in% "ZLCF") > 0) {
        if (B > 0) {
            ## ZLCF only for gehan
            vZLCF <- zlFun(beta = btemp2, Y = Y, delta = delta, X = X,
                           id = id, weights = weights, B = B, vClose = TRUE,
                           rankWeights = rankWeights, gw = rep(1, nrow(X)),
                           method = method, stratify = stratify, sigma = SigmaInit)
            ZLCF.An.inv <- vZLCF$An.inv
            vZLCF <- vZLCF$covmat
        }
    }
     
    if (sum(variance %in% "ZLMB") > 0) {
        if (B > 0) {
            vZLMB <- zlFun(beta = btemp2, Y = Y, delta = delta, X = X,
                           id = id, weights = weights, B = B, vClose = FALSE,
                           rankWeights = rankWeights, gw = gw,
                           method = method, stratify = stratify, sigma = SigmaInit)
            ZLMB.An.inv <- vZLMB$An.inv
            vZLMB <- vZLMB$covmat
            variance[which(variance == "sZLMB")] <- "ZLMB"
        }
    }

    if (sum(variance %in% "ISCF") > 0) {
        if (B > 0) {
            ## ISCF only for gehan
            vISCF <- isFun(beta = btemp2, Y = Y, delta = delta, X = X, id = id,
                           weights = weights, sigma = SigmaInit, B = B, vClose = TRUE,
                           gw = gw, rankWeights = rankWeights, stratify = stratify)
            ISCF.An.inv <- vISCF$An.inv
            vISCF <- vISCF$covmat
        }
    }

    if (sum(variance %in% "ISMB") > 0) {
        if (B > 0) {
            vISMB <- isFun(beta = btemp2, Y = Y, delta = delta, X = X, id = id,
                           weights = weights, sigma = SigmaInit, B = B, vClose = FALSE,
                           gw = gw, rankWeights = rankWeights)
            ISMB.An.inv <- vISMB$An.nv
            vISMB <- vISMB$covmat
        }
    }

    if (sum(variance %in% "sHCF") > 0) {
        if (B > 0) {
            vsHCF <- huangFun(beta = btemp2, Y = Y, delta = delta, X = X, id = id,
                              weights = weights, sigma = SigmaInit, B = B, vClose = TRUE,
                              gw = gw, rankWeights = rankWeights, stratify = stratify)$covmat
        }
    }

    if (sum(variance %in% "sHMB") > 0) {
        if (B > 0) {
            vsHMB <- huangFun(beta = btemp2, Y = Y, delta = delta, X = X, id = id,
                              weights = weights, sigma = SigmaInit, B = B, vClose = FALSE,
                              gw = gw, rankWeights = rankWeights, stratify = stratify)$covmat
        }
    }
    covmat <- list(MB = vMB, ZLCF = vZLCF, ZLMB = vZLMB,
                   sHCF = vsHCF, sHMB = vsHMB,
                   ISCF = vISCF, ISMB = vISMB)
    out <- list(beta = btemp2, covmat = covmat,
                convergence = bconv, bstep = bstep,
                var.meth = variance, bhist = bhist)
}


aftsrr.pe <- function(beta, Y, X, delta, clsize, sigma, weights, Z, 
                      rankWeights, method, pw = rep(1, nrow(X)),
                      control = aftgee.control()) {
    n <- nrow(X)
    bm <- gw <- NULL
    if (!is.numeric(pw)) {
        pw <- rep(1, nrow(X))
    }
    
    if (rankWeights == "gehan") {
        temp <- solveGehan(beta = beta, Y = Y, X = X, delta = delta, clsize = clsize,
                           sigma = sigma, weights = weights, Z = Z)
        bstep <- 1
        btime <- temp$tbeta
        btemp2 <- temp$beta
        bconv <- temp$conv
    }
    
    if (rankWeights != "gehan") {
        if (method == "nonsm") {
            btime <- system.time(temp <- BBsolve(beta, uilogFun, Y = Y, X = X, delta = delta,
                                                 clsize = clsize, sigma = sigma, n = n, Z = Z,
                                                 weights = weights, constant = 0,
                                                 pw = pw,
                                                 rankWeights = rankWeights, rkmethod = "nonsm",
                                                 quiet = TRUE))
            bstep <- 1
            btime <- btime[3]
            btemp2 <- temp$par
            bconv <- temp$convergence
        }
    }

    if (rankWeights != "gehan" & method != "nonsm") {
        if (rankWeights == "logrank" & method == "sm") {
            btime <- system.time(temp <- BBsolve(beta, uilogFun, Y = Y, X = X, delta = delta,
                                                 clsize = clsize, sigma = sigma, n = n, Z = Z,
                                                 weights = weights, constant = 0,
                                                 pw = rep(1, nrow(X)),
                                                 rankWeights = rankWeights, rkmethod = "sm",
                                                 quiet = TRUE))
            bstep <- 1
            btime <- btime[3]
            btemp2 <- temp$par
            bconv <- temp$convergence
        }   
        if (method == "monosm" | rankWeights %in% c("PW", "GP", "userdefined")) {
            btemp1 <- solveGehan(beta = beta, Y = Y, X = X, delta = delta,
                               clsize = clsize, sigma = sigma, weights = weights, Z = Z)$beta
            bstep <- 1
            btime <- 0
            for (i in 1:control$maxiter) {
                temp <- solvePW(beta = btemp1, Y = Y, X = X, delta = delta,
                                clsize = clsize, sigma = sigma, weights = weights,
                                Z = Z, rankWeights = rankWeights,
                                method = method, constant = 0)
                btemp2 <- temp$beta
                ebeta <- abs(btemp1 - btemp2)
                e.rel <- max(ebeta / abs(btemp2))
                bstep <- bstep + 1
                btime <- btime + temp$btime
                if (abs(e.rel) < control$reltol) {
                    bm <- c(bm, btemp2)
                    bconv <- 0
                    break
                }
                bm <- c(bm, btemp2)
                btemp1 <- btemp2
                bconv <- 1
            }            
        }
    }
    pe <- list(btemp2 = as.numeric(btemp2), btime = btime, bconv = bconv, bstep = bstep)
}


