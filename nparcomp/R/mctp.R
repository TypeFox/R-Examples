mctp<-function (formula, data, type = c("Tukey", "Dunnett", "Sequen", 
    "Williams", "Changepoint", "AVE", "McDermott", "Marcus", 
    "UmbrellaWilliams", "UserDefined"), conf.level = 0.95, alternative = c("two.sided", 
    "less", "greater"), asy.method = c("fisher", "mult.t", "normal"), 
    plot.simci = FALSE, control = NULL, info = TRUE, rounds = 3, 
    contrast.matrix = NULL, correlation = FALSE, effect = c("unweighted", 
        "weighted")) 
{
    input.list <- list(formula = formula, data = data, type = type[1], 
        conf.level = conf.level, alternative = alternative, asy.method = asy.method, 
        plot.simci = plot.simci, control = control, info = info, 
        rounds = rounds, contrast.matrix = contrast.matrix, correlation = correlation, 
        effect = effect)
    conflevel <- conf.level

    if (conflevel >= 1 || conflevel <= 0) {
        stop("The confidence level must be between 0 and 1!")
        if (is.null(alternative)) {
            stop("Please declare the alternative! (two.sided, lower, greater)")
        }
    }
    type <- match.arg(type)
    alternative <- match.arg(alternative)
    asy.method <- match.arg(asy.method)
    effect <- match.arg(effect)
    if (length(formula) != 3) {
        stop("You can only analyse one-way layouts!")
    }
    if (length(formula) != 3) {
        stop("You can only analyse one-way layouts!")
    }
    dat <- model.frame(formula, data)
    if (ncol(dat) != 2) {
        stop("Specify one response and only one class variable in the formula")
    }
    if (is.numeric(dat[, 1]) == FALSE) {
        stop("Response variable must be numeric")
    }
    response <- dat[, 1]
    factorx <- as.factor(dat[, 2])
    samples <- split(response, factorx)
    fl <- levels(factorx)
    a <- nlevels(factorx)
    if (a <= 2) {
        stop("You want to perform a two-sample test. Please use the function npar.t.test")
    }
    n <- sapply(samples, length)
    if (any(n <= 1)) {
        warn <- paste("The factor level", fl[n <= 1], "has got only one observation!")
        stop(warn)
    }
    N <- sum(n)
    tmp1 <- sort(rep(1:a, a))
    tmp2 <- rep(1:a, a)
    pairRanks <- lapply(1:(a^2), function(arg) rank(c(samples[[tmp1[arg]]], 
        samples[[tmp2[arg]]])))
    p <- sapply(1:(a^2), function(arg) {
        x1 <- samples[[tmp1[arg]]]
        x2 <- samples[[tmp2[arg]]]
        rx1x2 <- rank(c(x1, x2))
        l1 <- length(x1)
        l2 <- length(x2)
        1/(l1 + l2) * (mean(rx1x2[(l1 + 1):(l1 + l2)]) - mean(rx1x2[1:l1])) + 
            0.5
    })
    intRanks <- lapply(samples, rank)
    placements <- lapply(1:(a^2), function(arg) 1/n[tmp1[arg]] * 
        (pairRanks[[arg]][(n[tmp1[arg]] + 1):(n[tmp1[arg]] + 
            n[tmp2[arg]])] - intRanks[[tmp2[arg]]]))
    V <- rep(0, a^4)
    help <- expand.grid(1:a, 1:a, 1:a, 1:a)
    h1 <- help[, 4]
    h2 <- help[, 3]
    h3 <- help[, 2]
    h4 <- help[, 1]
    for (u in 1:(a^4)) {
        i <- h1[u]
        j <- h2[u]
        r <- h3[u]
        s <- h4[u]
        if (i == r && j == s && i != j && r != s) {
            xi <- samples[[i]]
            xj <- samples[[j]]
            ni <- length(xi)
            nj <- length(xj)
            ri <- rank(xi)
            rj <- rank(xj)
            rij <- rank(c(xi, xj))
            pj <- 1/ni * (rij[(ni + 1):(ni + nj)] - rj)
            pi <- 1/nj * (rij[1:ni] - ri)
            vi <- var(pi)/ni
            vj <- var(pj)/nj
            V[u] <- N * (vi + vj)
        }
        if (i == s && j == r && i != j && r != s) {
            xi <- samples[[i]]
            xj <- samples[[j]]
            ni <- length(xi)
            nj <- length(xj)
            ri <- rank(xi)
            rj <- rank(xj)
            rij <- rank(c(xi, xj))
            pj <- 1/ni * (rij[(ni + 1):(ni + nj)] - rj)
            pi <- 1/nj * (rij[1:ni] - ri)
            vi <- var(pi)/ni
            vj <- var(pj)/nj
            V[u] <- -N * (vi + vj)
        }
        if (i == r && j != s && i != j && r != s) {
            xi <- samples[[i]]
            xj <- samples[[j]]
            xs <- samples[[s]]
            ni <- length(xi)
            nj <- length(xj)
            ns <- length(xs)
            ri <- rank(xi)
            rj <- rank(xj)
            rs <- rank(xs)
            rij <- rank(c(xi, xj))
            ris <- rank(c(xi, xs))
            pij <- 1/nj * (rij[1:ni] - ri)
            pis <- 1/ns * (ris[1:ni] - ri)
            V[u] <- N * (cov(pij, pis)/ni)
        }
        if (i != r && j == s && i != j && r != s) {
            xi <- samples[[i]]
            xj <- samples[[j]]
            xr <- samples[[r]]
            ni <- length(xi)
            nj <- length(xj)
            nr <- length(xr)
            ri <- rank(xi)
            rj <- rank(xj)
            rr <- rank(xr)
            rji <- rank(c(xj, xi))
            rjr <- rank(c(xj, xr))
            pji <- 1/ni * (rji[1:nj] - rj)
            prj <- 1/nr * (rjr[1:nj] - rj)
            V[u] <- N * (cov(pji, prj)/nj)
        }
        if (i == s && j != r && i != j && r != s) {
            xi <- samples[[i]]
            xj <- samples[[j]]
            xr <- samples[[r]]
            ni <- length(xi)
            nj <- length(xj)
            nr <- length(xr)
            ri <- rank(xi)
            rj <- rank(xj)
            rr <- rank(xr)
            rij <- rank(c(xi, xj))
            rir <- rank(c(xi, xr))
            pij <- 1/nj * (rij[1:ni] - ri)
            pir <- 1/nr * (rir[1:ni] - ri)
            V[u] <- -N * (cov(pij, pir)/ni)
        }
        if (i != s && j == r && i != j && r != s) {
            xi <- samples[[i]]
            xj <- samples[[j]]
            xs <- samples[[s]]
            ni <- length(xi)
            nj <- length(xj)
            ns <- length(xs)
            ri <- rank(xi)
            rj <- rank(xj)
            rs <- rank(xs)
            rji <- rank(c(xj, xi))
            rjs <- rank(c(xj, xs))
            pji <- 1/ni * (rji[1:nj] - rj)
            pjs <- 1/ns * (rjs[1:nj] - rj)
            V[u] <- -N * (cov(pji, pjs)/nj)
        }
    }
    V1 <- matrix(V, ncol = a^2, nrow = a^2)
    switch(effect, weighted = {
        W <- kronecker(t(n/N), diag(a))
text.output.W <- paste("Global Ranks")
    }, unweighted = {
        W <- kronecker(t(rep(1/a, a)), diag(a))
text.output.W <- paste("Global Pseudo Ranks")
    })
    pd <- W %*% p
    VV <- W %*% V1 %*% t(W)
    logit.pd<-log(c(pd/(1-pd)))
    logit.pd.dev<-diag(1/c((pd*(1-pd))))

Lower.logit1 <-logit.pd-qnorm(conf.level)/sqrt(N)*sqrt(c(diag(logit.pd.dev%*%VV%*%t(logit.pd.dev))))
Upper.logit1 <-logit.pd+qnorm(conf.level)/sqrt(N)*sqrt(c(diag(logit.pd.dev%*%VV%*%t(logit.pd.dev))))
Lower.logit <- exp(Lower.logit1)/(1+exp(Lower.logit1))
Upper.logit <- exp(Upper.logit1)/(1+exp(Upper.logit1))
#Lower.logit <- pd-qnorm(1-(1-conf.level)/2)/sqrt(N)*sqrt(c(diag(VV)))
#Upper.logit <- pd+qnorm(1-(1-conf.level)/2)/sqrt(N)*sqrt(c(diag(VV)))
    if (type == "UserDefined") {
        if (is.null(contrast.matrix)) {
            stop("Please eanter a contrast matrix!")
        }
        ch <- contrast.matrix
        rownames(ch) <- paste("C", 1:nrow(ch))
        colnames(ch) <- fl
    }
    if (type != "UserDefined") {
        if (is.null(control)) {
            icon <- 1
        }
        if (!is.null(control)) {
            icon <- which(fl == control)
        }
        ch <- contrMat(n = n, type, base = icon)
    }
    nc <- nrow(ch)
    connames <- rownames(ch)
    Con <- matrix(ch, ncol = a)
    rownames(Con) <- connames
    colnames(Con) <- colnames(ch)
    degrees <- function(CC) {
        nc <- nrow(CC)
        dfs <- c()
        for (hhh in 1:nc) {
            cc <- CC[hhh, ]
            Yk <- list()
            for (l in 1:a) {
                Yk[[l]] <- 0
                for (i in 1:(a^2)) {
                  r <- tmp1[i]
                  s <- tmp2[i]
                  if (s == l && r != l) {
                    Ykhelp <- placements[[i]]
                    Yk[[l]] <- Yk[[l]] + Ykhelp
                  }
                }
            }
            Ykstern <- list()
            for (l in 1:a) {
                Ykstern[[l]] <- cc[l] * Yk[[l]]
                for (i in 1:(a^2)) {
                  r <- tmp1[i]
                  s <- tmp2[i]
                  if (s == l && r != l) {
                    Yksternhelp <- -cc[r] * placements[[i]]
                    Ykstern[[l]] <- Ykstern[[l]] + Yksternhelp
                  }
                }
            }
            variances <- sapply(Ykstern, var)/(a * n)
            varii2 <- (variances == 0)
            variances[varii2] <- 1/N
            dfs[hhh] <- (sum(variances))^2/sum(variances^2/(n - 
                1))
        }
        dfs
    }
    dfT <- round(max(4, min(degrees(Con))))
    Cpd <- Con %*% pd
    CV <- Con %*% VV %*% t(Con)
    rhobf <- cov2cor(CV)
    p.adj <- c()
    switch(asy.method, mult.t = {
        T <- sqrt(N) * (Cpd)/sqrt(c(diag(CV)))
        AsyMethod <- paste("Multi - T with", round(dfT, rounds), 
            "DF")
        switch(alternative, two.sided = {
            text.Output <- paste("True differences of relative effects are less or equal than 0")
            for (pp in 1:nc) {
                p.adj[pp] <- 1 - pmvt(lower = -abs(T[pp]), upper = abs(T[pp]), 
                  delta = rep(0, nc), df = dfT, corr = rhobf)[1]
            }
            crit <- qmvt(conflevel, corr = rhobf, tail = "both", 
                df = dfT)$quantile
            Lower <- Cpd - crit/sqrt(N) * sqrt(c(diag(CV)))
            Upper <- Cpd + crit/sqrt(N) * sqrt(c(diag(CV)))
        }, less = {
            text.Output <- paste("True differences of relative effects are less than 0")
            for (pp in 1:nc) {
                p.adj[pp] <- 1-pmvt(lower = T[pp], upper = Inf, 
                  df = dfT, delta = rep(0, nc), corr = rhobf)
            }
            crit <- qmvt(conflevel, df = dfT, corr = rhobf, tail = "lower")$quantile
            Lower <- rep(-1, nc)
            Upper <- Cpd + crit/sqrt(N) * sqrt(c(diag(CV)))
        }, greater = {
            text.Output <- paste("True differences of relative effects are greater than 0")
            for (pp in 1:nc) {
                p.adj[pp] <- 1 - pmvt(lower = -Inf, upper = T[pp], 
                  df = dfT, delta = rep(0, nc), corr = rhobf)[1]
            }
            crit <- qmvt(conflevel, corr = rhobf, df = dfT, tail = "lower")$quantile
            Lower <- Cpd - crit/sqrt(N) * sqrt(c(diag(CV)))
            Upper <- rep(1, nc)
        })
    }, normal = {
        AsyMethod <- "Normal - Approximation"
        T <- sqrt(N) * (Cpd)/sqrt(c(diag(CV)))
        switch(alternative, two.sided = {
            text.Output <- paste("True differences of relative effects are less or equal than 0")
            for (pp in 1:nc) {
                p.adj[pp] <- 1 - pmvnorm(lower = -abs(T[pp]), 
                  upper = abs(T[pp]), mean = rep(0, nc), corr = rhobf)[1]
            }
            crit <- qmvnorm(conflevel, corr = rhobf, tail = "both")$quantile
            Lower <- Cpd - crit/sqrt(N) * sqrt(c(diag(CV)))
            Upper <- Cpd + crit/sqrt(N) * sqrt(c(diag(CV)))
        }, less = {
            text.Output <- paste("True differences of relative effects are less than 0")
            for (pp in 1:nc) {
                p.adj[pp] <- 1-pmvnorm(lower = T[pp], upper = Inf, 
                  mean = rep(0, nc), corr = rhobf)
            }
            crit <- qmvnorm(conflevel, corr = rhobf, tail = "lower")$quantile
            Lower <- rep(-1, nc)
            Upper <- Cpd + crit/sqrt(N) * sqrt(c(diag(CV)))
        }, greater = {
            text.Output <- paste("True differences of relative effects are greater than 0")
            for (pp in 1:nc) {
                p.adj[pp] <- 1 - pmvnorm(lower = -Inf, upper = T[pp], 
                  mean = rep(0, nc), corr = rhobf)
            }
            crit <- qmvnorm(conflevel, corr = rhobf, tail = "lower")$quantile
            Lower <- Cpd - crit/sqrt(N) * sqrt(c(diag(CV)))
            Upper <- rep(1, nc)
        })
    }, fisher = {
        AsyMethod <- paste("Fisher with", round(dfT, rounds), 
            "DF")
        Cfisher <- 1/2 * log((1 + Cpd)/(1 - Cpd))
        Vfisherdev <- diag(c(1/(1 - Cpd^2)))

        Vfisher <- Vfisherdev %*% CV %*% t(Vfisherdev)
        T <- sqrt(N) * Cfisher/sqrt(c(diag(Vfisher)))
        switch(alternative, two.sided = {
            text.Output <- paste("True differences of relative effects are less or equal than 0")
            for (pp in 1:nc) {
                p.adj[pp] <- 1 - pmvt(lower = -abs(T[pp]), upper = abs(T[pp]), 
                  delta = rep(0, nc), corr = rhobf, df = dfT)[1]
            }
            crit <- qmvt(conflevel, corr = rhobf, tail = "both", 
                df = dfT)$quantile
            Lower1 <- Cfisher - crit/sqrt(N) * sqrt(c(diag(Vfisher)))
            Upper1 <- Cfisher + crit/sqrt(N) * sqrt(c(diag(Vfisher)))
            Lower <- (exp(2 * Lower1) - 1)/(exp(2 * Lower1) + 1)
            Upper <- (exp(2 * Upper1) - 1)/(exp(2 * Upper1) + 1)
        }, less = {
            text.Output <- paste("True differences of relative effects are less than 0")
            for (pp in 1:nc) {
                p.adj[pp] <- 1-pmvt(lower = T[pp], upper = Inf, 
                  delta = rep(0, nc), df = dfT, corr = rhobf)
            }
            crit <- qmvt(conflevel, corr = rhobf, tail = "lower", 
                df = dfT)$quantile
            Lower <- rep(-1, nc)
            Upper1 <- Cfisher + crit/sqrt(N) * sqrt(c(diag(Vfisher)))
            Upper <- (exp(2 * Upper1) - 1)/(exp(2 * Upper1) + 
                1)
        }, greater = {
            text.Output <- paste("True differences of relative effects are greater than 0")
            for (pp in 1:nc) {
                p.adj[pp] <- 1 - pmvt(lower = -Inf, upper = T[pp], 
                  delta = rep(0, nc), corr = rhobf, df = dfT)
            }
            crit <- qmvnorm(conflevel, corr = rhobf, tail = "lower")$quantile
            Lower1 <- Cfisher - crit/sqrt(N) * sqrt(c(diag(Vfisher)))
            Lower <- (exp(2 * Lower1) - 1)/(exp(2 * Lower1) + 
                1)
            Upper <- rep(1, nc)
        })
    })
    if (plot.simci == TRUE) {
        text.Ci <- paste(conflevel * 100, "%", "Simultaneous Confidence Intervals")
        Lowerp <- "|"
        plot(Cpd, 1:nc, xlim = c(-1, 1), pch = 15, axes = FALSE, 
            xlab = "", ylab = "")
        points(Lower, 1:nc, pch = Lowerp, font = 2, cex = 2)
        points(Upper, 1:nc, pch = Lowerp, font = 2, cex = 2)
        abline(v = 0, lty = 3, lwd = 2)
        for (ss in 1:nc) {
            polygon(x = c(Lower[ss], Upper[ss]), y = c(ss, ss), 
                lwd = 2)
        }
        axis(1, at = seq(-1, 1, 0.1))
        axis(2, at = 1:nc, labels = connames)
        box()
        title(main = c(text.Ci, paste("Type of Contrast:", type), 
            paste("Method:", AsyMethod)))
    }
    data.info <- data.frame(Sample = fl, Size = n, Effect = pd, 
Lower=Lower.logit, Upper=Upper.logit)
    Analysis.of.Relative.Effects <- data.frame(Estimator = round(Cpd, 
        rounds), Lower = round(Lower, rounds), Upper = round(Upper, 
        rounds), Statistic = round(T, rounds), p.Value = p.adj)
    Overall <- data.frame(Quantile = crit, p.Value = min(p.adj))
    result <- list(Data.Info = data.info, Contrast = Con, Analysis = Analysis.of.Relative.Effects, 
        Overall = Overall)
    if (info == TRUE) {
        cat("\n", "#----------------Nonparametric Multiple Comparisons for relative effects---------------#", 
            "\n", "\n", "-", "Alternative Hypothesis: ", text.Output, 
            "\n", "-", "Estimation Method: ", text.output.W,
            "\n", "-", "Type of Contrast", ":", type, "\n", "-", 
            "Confidence Level:", conflevel * 100, "%", "\n", 
            "-", "Method", "=", AsyMethod, "\n", "\n", "#--------------------------------------------------------------------------------------#", 
            "\n", "\n")
    }
    if (correlation == TRUE) {
        result$Covariance <- CV
        result$Correlation <- rhobf
    }
    result$input <- input.list
    result$text.Output <- text.Output
    result$text.output.W <- text.output.W
    result$connames <- connames
    result$AsyMethod <- AsyMethod
    class(result) <- "mctp"

    return(result)
}

