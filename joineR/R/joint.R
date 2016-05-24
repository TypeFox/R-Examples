joint <- function (data, long.formula, surv.formula, model = c("intslope", 
    "int", "quad"), sepassoc = FALSE, longsep = FALSE, survsep = FALSE, 
    gpt, lgpt, max.it, tol) 
{
    id <- data$subj.col
    time.long <- data$time.col
    if (missing(gpt)) {
        gpt <- 3
    }
    if (missing(lgpt)) {
        lgpt <- 10
    }
    if (missing(max.it)) {
        max.it <- 200
    }
    if (missing(tol)) {
        tol <- 0.001
    }
    Call <- match.call()
    if (any(sapply(data$baseline, "class") == "factor")) {
    data$baseline <- drop.levels(data$baseline)
    }
    long.data <- merge(data$longitudinal, data$baseline, by = id, 
        sort = FALSE)
    long.frame <- model.frame(long.formula, data = long.data)
    long.cov <- model.matrix(long.formula, long.frame)
    long.terms <- terms(long.formula, data = long.data)
    long.names <- colnames(long.cov)
    rll <- !is.na(data$longitudinal[[names(long.frame[1])]])
    longdat <- cbind(data$longitudinal[[id]][rll], long.frame[, 
        1], data$longitudinal[[time.long]][rll], long.cov)
    longdat <- as.data.frame(longdat)
    names(longdat) <- c(id, names(long.frame)[1], time.long, 
        long.names)
    surv.frame <- model.frame(surv.formula, data = cbind(data$survival, 
        data$baseline))
    srv <- model.extract(surv.frame, "response")
    surv.terms <- terms(surv.formula, data = cbind(data$survival, 
        data$baseline))
    attr(surv.terms, "intercept") <- 1
    surv.cov <- model.matrix(surv.terms, data = cbind(data$survival, 
        data$baseline))
    surv.cov <- as.matrix(surv.cov[, -1])
    rss <- as.integer(row.names(surv.cov))
    survdat <- cbind(data$survival[[id]][rss], srv[rss, 1], srv[rss, 
        2], surv.cov[rss, ])
    survdat <- as.data.frame(survdat)
    names(survdat) <- c(id, surv.formula[2][[1]][[2]], surv.formula[2][[1]][[3]], 
	   colnames(surv.cov))    
    if (dim(survdat)[2] > 3) {
        survdat[, 4:dim(survdat)[2]] <- scale(survdat[, 4:dim(survdat)[2]], 
            scale = FALSE)
    }
    survdat2 <- data.frame(data$survival[[id]][rss], srv[rss, 1], srv[rss, 
        2], surv.frame[, -1])
    names(survdat2) <- c(id, surv.formula[2][[1]][[2]], surv.formula[2][[1]][[3]], 
          attr(surv.terms, "term.labels"))
    model <- match.arg(model)
    if (model != "intslope" && model != "int" && model != "quad") {
        stop(paste("Unknown model", model))
    }
    ran <- 2
    if (model == "int") {
        ran <- 1
    }
    if (model == "quad") {
        ran <- 3
    }
    lat <- ran
    if (!sepassoc) {
        lat <- 1
    }
    sep <- function(ests, logical) {
        if (logical == FALSE) {
            ests <- "No separate results requested"
        }
        ests
    }
    longst <- function(longdat, long.formula, model, longdat2) {
        if (model == "int") {
            rf <- as.formula(paste("~1", colnames(longdat)[1], 
                sep = "|"))
            long.start <- lme(long.formula, random = rf, method = "ML", 
                data = data.frame(longdat2), na.action = na.omit)
        }
        else if (model == "intslope") {
            rf <- as.formula(paste(paste("~", colnames(longdat)[3], 
                sep = ""), colnames(longdat)[1], sep = "|"))
            long.start <- lme(long.formula, random = rf, method = "ML", 
                data = data.frame(longdat2), na.action = na.omit)
        }
        else {
            tsq <- paste(paste("I(", paste(colnames(longdat)[3], 
                "^2", sep = ""), sep = ""), ")", sep = "")
            rf <- as.formula(paste(paste("~", paste(colnames(longdat)[3], 
                tsq, sep = "+"), sep = ""), colnames(longdat)[1], 
                sep = "|"))
            long.start <- lme(long.formula, random = rf, method = "ML", 
                data = data.frame(longdat2), na.action = na.omit)
        }
        q <- dim(VarCorr(long.start))[1] - 1
        sigma.u <- diag(as.double(VarCorr(long.start)[1:q, 1]), 
            q, q)
        if (q > 1) {
            vv <- tcrossprod(diag(sqrt(sigma.u))) * lower.tri(sigma.u)
            rho <- as.double(VarCorr(long.start, rdig = 8)[-c(1, 
                q + 1), -(1:2)])
            rho <- rho[!is.na(rho)]
            vars <- diag(sigma.u)
            sigma.u[lower.tri(sigma.u)] <- vv[lower.tri(vv)] * 
                rho
            sigma.u <- sigma.u + t(sigma.u) - diag(vars)
        }
        rownames(sigma.u) <- paste("U_", 0:(q - 1), sep = "")
        colnames(sigma.u) <- paste("U_", 0:(q - 1), sep = "")
        sigma.z <- long.start$sigma^2
        ll <- long.start$logLik
        b1 <- fixef(long.start)
        list(b1 = data.frame(b1), sigma.z = sigma.z, sigma.u = sigma.u, 
            log.like = ll)
    }
    survst <- function(survdat, surv.formula, survdat2) {
	survdat2 <- survdat2[order(survdat2[, 2]), ]
        n <- length(survdat[, 2])
        s <- survdat[, 2]
        cen <- survdat[, 3]
        cen[1] <- 1
        survdat[1, 3] <- 1
        survdat2[1, 3] <- 1
        p2 <- dim(survdat)[2] - 3
        surv.start <- coxph(surv.formula, data = survdat2, x = TRUE)
        surv.start.f <- survfit(surv.start)
        sf <- surv.start.f$time[surv.start.f$n.event != 0]
        nf <- length(sf)
        nev <- surv.start.f$n.event[surv.start.f$n.event != 0]
        if (p2 > 0) {
            haz <- coxph.detail(surv.start)$hazard
        }
        else {
            haz <- surv.start.f$n.event/surv.start.f$n.risk
            haz <- haz[surv.start.f$n.event > 0]
        }
        rs <- rep(1:nf, c(diff(match(sf, s)), n + 1 - match(sf, 
            s)[nf]))
        if (cen[1] == 0) {
            rs <- c(0, rs)
        }
        b2 <- coef(surv.start)
        ll <- surv.start$loglik - sum(cen)
        list(b2 = b2, haz = haz, rs = rs, sf = sf, nev = nev, 
            log.like = ll)
    }
    em.alg <- function(longdat, survdat, ran, paraests, gpt, 
        max.it, tol) {
        id <- longdat[, 1]
        Y <- longdat[, 2]
        tt <- longdat[, 3]
        X1 <- as.matrix(longdat[, 4:dim(longdat)[2]])
        n <- length(survdat[, 2])
        s <- survdat[, 2]
        cen <- survdat[, 3]
        p1 <- dim(longdat)[2] - 3
        p2 <- dim(survdat)[2] - 3
        X2 <- 0
        if (p2 > 0) {
            X2 <- as.matrix(survdat[, 4:dim(survdat)[2]])
        }
        else {
            b2x <- matrix(0, n, 1)
        }
        b1 <- paraests$b1[, 1]
        sigma.u <- paraests$sigma.u
        tsigu <- t(sigma.u)
        sigma.z <- paraests$sigma.z
        b2 <- c(paraests$b2, rep(0, lat))
        haz <- paraests$haz
        sf <- paraests$sf
        rs <- paraests$rs
        nev <- paraests$nev
        nn <- diff(match(unique(id), id))
        nn <- c(nn, length(id) - sum(nn))
        N <- sum(nn)
        g <- gauss.quad.prob(gpt, "normal", sigma = sqrt(0.5))
        ab <- g$nodes
        w <- g$weights * sqrt(pi)
        gmat <- matrix(0, gpt^ran, ran)
        gmat[, 1] <- rep(ab, each = gpt^(ran - 1))
        if (model != "int") {
            gmat[, 2] <- rep(ab, gpt)
            w <- as.vector(w %x% w)
        }
        if (model == "quad") {
            gmat[, 3] <- rep(ab, each = gpt)
            w <- as.vector(w %x% g$weights * sqrt(pi))
        }
        EU <- matrix(0, n, ran)
        EUU <- matrix(0, n, sum(1:ran))
        EexpU <- matrix(0, n, length(haz))
        EUexpU <- matrix(0, n, ran)
        EUUexpU <- matrix(0, n, sum(1:ran))
        r <- Y - X1 %*% b1
        Dtt <- getD(ran, tt)
        Dtt2 <- t(Dtt)
        if (model != "int") {
            Dttc <- t(getD(sum(1:ran) - ran, tt)) * tt
        }
        Ds <- getD(ran, s)
        Dst <- t(Ds)
        Dsf <- getD(ran, sf)
        Dsf2 <- Dsf^2
        Dsfc <- t(t(Dsf) * sf)
        Dnsf <- matrix(1, ran, length(sf))
        s1 <- rep(1:(ran - 1), (ran - 1):1)
        s2 <- sequence((ran - 1):1) + rep(1:(ran - 1), (ran - 
            1):1)
        cnn <- c(0, cumsum(nn))
        Inn <- diag(max(nn))
        conv <- FALSE
        for (it in 1 : max.it) {
            if (p2 > 0) {
                b2x <- X2 %*% b2[1:p2]
            }
            eb2x <- exp(b2x)
            sigma.zi <- sigma.z * Inn
            cov <- sigma.u %*% Dtt
            tcov <- Dtt2 %*% sigma.u
            DH <- Dnsf * rep(haz, each = ran)
            for (i in 1:n) {
                rv <- r[(cnn[i] + 1):cnn[i + 1]]
                ttv <- Dtt2[(cnn[i] + 1):cnn[i + 1], ]
                W21 <- cov[, (cnn[i] + 1):cnn[i + 1]]
                W12 <- tcov[(cnn[i] + 1):cnn[i + 1], ]
                if (model == "int") {
                  W11 <- tcrossprod(ttv, W21) + sigma.zi[1:nn[i], 
                    1:nn[i]]
                }
                else {
                  W11 <- ttv %*% W21 + sigma.zi[1:nn[i], 1:nn[i]]
                }
                if (nn[i] == 1) {
                  W3 <- W12/W11
                  if (model == "int") {
                    cvch <- sqrt((sigma.u - tcrossprod(W21, W3)) * 
                      2)
                  }
                  else {
                    cvch <- chol((sigma.u - tcrossprod(W21, W3)) * 
                      2)
                  }
                  cm <- matrix(W3 * rv, gpt^ran, ran, TRUE)
                }
                else {
                  W3 <- solve(W11, W12)
                  if (model == "int") {
                    cvch <- sqrt((sigma.u - W21 %*% W3) * 2)
                  }
                  else {
                    cvch <- chol((sigma.u - W21 %*% W3) * 2)
                  }
                  cm <- matrix(rv %*% W3, gpt^ran, ran, TRUE)
                }
                newu <- gmat %*% cvch + cm
                newu2 <- newu^2
                if (model != "int") {
                  newu2 <- cbind(newu2, newu[, s1] * newu[, s2])
                }
                egDUs <- 1
                if (cen[i] == 1) {
                  egDUs <- exp(newu %*% (Dst[i, ] * b2[(p2 + 
                    1):(p2 + lat)]))
                }
                egDUsf <- exp(newu %*% (Dsf[, 1:rs[i]] * b2[(p2 + 
                  1):(p2 + lat)]))
                ess <- exp(-(eb2x[i, ] * egDUsf) %*% haz[1:rs[i]])
                f <- egDUs * ess * w
                den <- sum(f)
                EU[i, 1:ran] <- f[, 1] %*% newu/den
                EUU[i, 1:sum(1:ran)] <- f[, 1] %*% newu2/den
                C <- egDUsf[, 1:rs[i]]
                EexpU[i, 1:rs[i]] <- f[, 1] %*% C/den
                if (model == "int") {
                  EUexpU[i, 1] <- sum(f[, 1] %*% (newu[, 1] * 
                    C) * haz[1:rs[i]])/den
                  EUUexpU[i, 1] <- sum(f[, 1] %*% (newu[, 1]^2 * 
                    C) * haz[1:rs[i]])/den
                }
                else {
                  EUexpU[i, 1:ran] <- rowSums(crossprod(newu * 
                    f[, 1], C) * Dsf[, 1:rs[i]] * DH[, 1:rs[i]])/den
                  EUUexpU[i, 1:ran] <- rowSums(crossprod(newu2[, 
                    1:ran] * f[, 1], C) * Dsf2[, 1:rs[i]] * DH[, 
                    1:rs[i]])/den
                  if (model == "intslope") {
                    EUUexpU[i, ran + 1] <- 2 * sum(f[, 1] %*% 
                      (newu2[, ran + 1] * C) * haz[1:rs[i]] * 
                      sf[1:rs[i]])/den
                  }
                  else {
                    EUUexpU[i, (ran + 1):sum(1:ran)] <- 2 * rowSums(crossprod(newu2[, 
                      (ran + 1):sum(1:ran)] * f[, 1], C) * Dsfc[, 
                      1:rs[i]] * DH[, 1:rs[i]])/den
                  }
                }
            }
            parac <- data.frame(c(b1, b2, sigma.z, sigma.u))
            EexpUi <- colSums(t(EexpU) * haz)
            haz <- nev/colSums(EexpU * eb2x[, 1])
            EUmat <- apply(EU, 2, rep, nn)
            EUUmat <- apply(EUU, 2, rep, nn)
            Ut <- rowSums(EUmat * Dtt2)
            UUt <- rowSums(EUUmat[, 1:ran] * Dtt2^2)
            UUt2 <- 0
            if (model != "int") {
                UUt2 <- rowSums(EUUmat[, (ran + 1):sum(1:ran)] * 
                  Dttc)
            }
            b1 <- solve(crossprod(X1), crossprod(X1, Y - Ut))
            r <- Y - X1 %*% b1
            sigma.z <- sum(r^2 - 2 * r * Ut + UUt + 2 * UUt2)/N
            diag(sigma.u) <- colMeans(EUU)[1:ran]
            if (model != "int") {
                sigma.u[lower.tri(sigma.u)] <- colMeans(EUU)[-(1:ran)]
                sigma.u[upper.tri(sigma.u)] <- t(sigma.u)[upper.tri(sigma.u)]
            }
            fd <- vector("numeric", p2 + ran)
            sd <- matrix(0, p2 + ran, p2 + ran)
            fd[(p2 + 1):(p2 + ran)] <- colSums(cen * (EU * t(Ds))) - 
                colSums(eb2x[, 1] * EUexpU)
            if (model != "int") {
                sd[(p2 + 1):(p2 + ran), (p2 + 1):(p2 + ran)][upper.tri(sd[(p2 + 
                  1):(p2 + ran), (p2 + 1):(p2 + ran)])] <- -colSums(eb2x[, 
                  1] * 0.5 * EUUexpU)[(ran + 1):sum(1:ran)]
            }
            if (p2 > 0) {
                fd[1:p2] <- c(colSums((cen * X2) - (X2 * eb2x[, 
                  1] * EexpUi)))
                sd[(1:p2), (p2 + 1):(p2 + ran)] <- -t(X2) %*% 
                  (eb2x[, 1] * EUexpU)
                sd <- sd + t(sd)
                for (i in 1:p2) {
                  for (j in 1:p2) {
                    sd[i, j] <- -(sum(X2[, i] * X2[, j] * eb2x[, 
                      1] * EexpUi))
                  }
                }
            }
            if (model == "int") {
                sd[(p2 + 1), (p2 + 1)] <- -colSums(eb2x[, 1] * 
                  EUUexpU)[1:ran]
            }
            else {
                diag(sd[(p2 + 1):(p2 + ran), (p2 + 1):(p2 + ran)]) <- -colSums(eb2x[, 
                  1] * EUUexpU)[1:ran]
            }
            if (!sepassoc) {
                if (model == "int") {
                  fd <- fd
                  sd <- sd
                }
                else {
                  fd[p2 + 1] <- sum(fd[(p2 + 1):(p2 + ran)])
                  fd <- fd[1:(p2 + 1)]
                  if (p2 > 1) {
                    sd[1:p2, p2 + 1] <- rowSums(sd[(1:p2), (p2 + 
                      1):(p2 + ran)])
                  }
                  else {
                    sd[1:p2, p2 + 1] <- sum(sd[(1:p2), (p2 + 
                      1):(p2 + ran)])
                  }
                  sd[p2 + 1, 1:p2] <- sd[1:p2, p2 + 1]
                  sd[p2 + 1, p2 + 1] <- sum(sd[(p2 + 1):(p2 + 
                    ran), (p2 + 1):(p2 + ran)])
                  sd <- sd[1:(p2 + 1), 1:(p2 + 1)]
                }
            }
            b2 <- b2 - solve(sd, fd)
            para <- data.frame(c(b1, b2, sigma.z, sigma.u))
            dd <- abs(parac - para)
            if (max(dd) < tol) {
                conv <- TRUE
                break
            }
        }
        if (conv != TRUE) {
            print("Not converged")
        }
        list(b1 = data.frame(b1), b2 = data.frame(b2), sigma.z = sigma.z, 
            sigma.u = sigma.u, haz = haz, random = EU, conv = conv, 
            iters = it)
    }
    getD <- function(q, arg) {
        D <- matrix(0, q, length(arg))
        for (i in 1:q) {
            D[i, ] <- arg^(i - 1)
        }
        D
    }
    jlike <- function(longdat, survdat, ran, likeests, lgpt) {
        id <- longdat[, 1]
        Y <- longdat[, 2]
        tt <- longdat[, 3]
        X1 <- as.matrix(longdat[, 4:dim(longdat)[2]])
        n <- length(survdat[, 2])
        s <- survdat[, 2]
        cen <- survdat[, 3]
        nn <- diff(match(unique(id), id))
        nn[length(nn) + 1] <- length(id) - sum(nn)
        p1 <- dim(longdat)[2] - 3
        p2 <- dim(survdat)[2] - 3
        X2 <- 0
        if (p2 > 0) {
            X2 <- as.matrix(survdat[, 4:dim(survdat)[2]])
        }
        b1 <- likeests$b1[, 1]
        sigma.u <- likeests$sigma.u
        sigma.z <- likeests$sigma.z
        b2 <- likeests$b2[, 1]
        haz <- likeests$haz
        sf <- likeests$sf
        rs <- likeests$rs
        N <- sum(nn)
        g <- gauss.quad.prob(lgpt, "normal", sigma = sqrt(0.5))
        ab <- g$nodes
        w <- g$weights * sqrt(pi)
        gmat <- matrix(0, lgpt^ran, ran)
        gmat[, 1] <- rep(ab, each = lgpt^(ran - 1))
        if (model != "int") {
            gmat[, 2] <- rep(ab, lgpt)
            w <- as.vector(w %x% w)
        }
        if (model == "quad") {
            gmat[, 3] <- rep(ab, each = lgpt)
            w <- as.vector(w %x% g$weights * sqrt(pi))
        }
        l1 <- 0
        l2 <- 0
        r <- Y - X1 %*% b1
        Dtt <- getD(ran, tt)
        Dtt2 <- t(Dtt)
        Ds <- getD(ran, s)
        Dst <- t(Ds)
        Dsf <- getD(ran, sf)
        cnn <- c(0, cumsum(nn))
        b2x <- X2 %*% b2[1:p2]
        if (p2 == 0) {
            b2x <- matrix(0, n, 1)
        }
        sigma.zi <- sigma.z * diag(max(nn))
        cov <- sigma.u %*% Dtt
        tcov <- Dtt2 %*% sigma.u
        for (i in 1:n) {
            rv <- r[(cnn[i] + 1):cnn[i + 1]]
            ttv <- Dtt2[(cnn[i] + 1):cnn[i + 1], ]
            W21 <- cov[, (cnn[i] + 1):cnn[i + 1]]
            W12 <- tcov[(cnn[i] + 1):cnn[i + 1], ]
            if (model == "int") {
                W11 <- tcrossprod(ttv, W21) + sigma.zi[1:nn[i], 
                  1:nn[i]]
            }
            else {
                W11 <- ttv %*% W21 + sigma.zi[1:nn[i], 1:nn[i]]
            }
            if (nn[i] == 1) {
                W3 <- W12/W11
                if (model == "int") {
                  cvch <- sqrt((sigma.u - tcrossprod(W21, W3)) * 
                    2)
                }
                else {
                  cvch <- chol((sigma.u - tcrossprod(W21, W3)) * 
                    2)
                }
                cm <- matrix(W3 * rv, lgpt^ran, ran, TRUE)
            }
            else {
                W3 <- solve(W11, W12)
                if (model == "int") {
                  cvch <- sqrt((sigma.u - W21 %*% W3) * 2)
                }
                else {
                  cvch <- chol((sigma.u - W21 %*% W3) * 2)
                }
                cm <- matrix(rv %*% W3, lgpt^ran, ran, TRUE)
            }
            newu <- gmat %*% cvch + cm
            DUs <- newu %*% Ds[, i]
            DUsf <- newu %*% Dsf[, 1:rs[i]]
            ss <- exp(b2x[i, ]) * exp(newu %*% (Dsf[, 1:rs[i]] * 
                b2[(p2 + 1):(p2 + lat)])) %*% haz[1:rs[i]]
            den <- sum(exp(cen[i] * (newu %*% (Dst[i, ] * b2[(p2 + 
                1):(p2 + lat)]) + b2x[i, ])) * (haz[rs[i]]^cen[i]) * 
                w * exp(-ss))
            l2 <- l2 + 0
            if (den > 0) {
                l2 <- l2 + log(den)
            }
            l1 <- l1 - nn[i] * 0.5 * log(2 * pi) - 0.5 * log(det(W11)) - 
                0.5 * sum(rv * solve(W11, rv))
        }
        ll <- l1 + l2 - 0.5 * ran * n * log(pi)
        list(log.like = ll, longlog.like = l1, survlog.like = ll - 
            l1)
    }
    sort.dat <- function(longdat, survdat) {
        longid <- longdat[, 1]
        nn <- diff(match(unique(longid), longid))
        nn[length(nn) + 1] <- length(longid) - sum(nn)
        svec <- rep(survdat[, 2], nn)
        sort.long <- longdat[order(svec), ]
        os <- order(survdat[, 2])
        sort.surv <- survdat[os, ]
        list(long.s = data.frame(sort.long), surv.s = data.frame(sort.surv))
    }
    sort <- sort.dat(longdat, survdat)
    longdat <- as.matrix(sort$long.s)
    survdat <- as.matrix(sort$surv.s)
    p2 <- dim(survdat)[2] - 3
    ldaests <- longst(longdat, long.formula, model = model, long.data)
    survests <- survst(survdat, surv.formula, survdat2)
    sep.ll <- ldaests$log.like + survests$log.like[2]
    sep.loglik <- list(seplhood = sep.ll, sepy = ldaests$log.like, 
        sepn = survests$log.like[2])
    paraests <- c(ldaests, survests)
    jointfit <- em.alg(longdat, survdat, ran, paraests, gpt, 
        max.it, tol)
    likeests <- c(jointfit, list(rs = survests$rs, sf = survests$sf))
    b1 <- jointfit$b1
    sigma.u <- jointfit$sigma.u
    rownames(b1) <- rownames(paraests$b1)
    if (p2 > 0) {
        b2 <- jointfit$b2[1:p2, ]
        names(b2) <- names(paraests$b2)
    }
    else {
        b2 <- NULL
    }
    fixed <- list(longitudinal = b1, survival = b2)
    latent <- jointfit$b2[(p2 + 1):(p2 + lat), ]
    names(latent) <- paste("gamma_", 0:(lat - 1), sep = "")
    random <- jointfit$random
    colnames(random) <- paste("U_", 0:(ran - 1), sep = "")
    rownames(random) <- survdat[, 1]
    coefficients <- list(fixed = fixed, random = random, latent = latent)
    jointll <- jlike(longdat, survdat, ran, likeests, lgpt)
    loglik <- list(jointlhood = jointll$log.like, jointy = jointll$longlog.like, 
        jointn = jointll$survlog.like)
    sepests <- list(longests = sep(ldaests, longsep), survests = sep(survests, 
        survsep))
    results <- list(coefficients = coefficients, sigma.z = jointfit$sigma.z, 
        sigma.u = jointfit$sigma.u, hazard = jointfit$haz, loglik = loglik, 
        numIter = jointfit$iters, convergence = jointfit$conv, 
        model = model, sepassoc = sepassoc, sepests = sepests, 
        sep.loglik = sep.loglik, formulae = list(lformula = long.formula, 
            sformula = surv.formula), data = data, call = Call)
    class(results) <- "joint"
    results
}
