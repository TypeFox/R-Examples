gao <-
function (formula, data, alpha = 0.05, control = NULL, silent = FALSE)
{

    dat <- model.frame(formula, data)
    if (ncol(dat) != 2) {
        stop("Specify one response and only one class variable in the formula !")
    }
    if (is.numeric(dat[, 1]) == FALSE) {
        stop("Response variable must be numeric !")
    }
    response <- dat[, 1]
    group <- as.factor(dat[, 2])
    fl <- levels(group)
    a <- nlevels(group)
    N <- length(response)
    n <- aggregate(response, list(group), FUN = "length")$x
    if (any(n <= 2)) {
        warn <- paste("The factor level", fl[n <= 2], "has got less than two observations!")
        stop(warn)
    }
    if (is.null(control)) {
        cont <- 1
    }
    if (!is.null(control)) {
        if (!any(fl == control)) {
            stop("The dataset doesn't contain this control group!")
        }
        cont <- which(fl == control)
    }
    C <- contrMat(n=c(rep(10,a)), "Dunnett", base = cont)
    rx <- c()
    for (i in 1:N) {
        help <- expand.grid(response[i], response)
        help1 <- (help[, 1] > help[, 2]) + 1/2 * (help[, 1] ==
            help[, 2])
        help2 <- data.frame(h1 = help1, h2 = group)
        samples2 <- split(help2$h1, help2$h2)
        pseudo <- sapply(1:a, function(arg) {
            mean(samples2[[arg]])
        })
        rx[i] <- N * mean(pseudo)
    }
    new.data <- data.frame(res = rx, group = group)
    pd <- 1/N * aggregate(new.data$res, list(group), FUN = "mean")$x
    Cpd <- C %*% pd
    v1 <- 1/N^2 * aggregate(new.data$res, list(group), FUN = "var")$x
    lambda <- N/n
    v11 <- c(v1 * lambda)
    v2 <- diag(v1 * lambda)
    Cv <- C %*% v2 %*% t(C)
    T <- sqrt(N) * Cpd/sqrt(c(diag(Cv)))
    ncont <- which((1:a) != cont)
    numerator <- c(diag(Cv))^2
    denu1 <- v1[cont]^2/(n[cont]^2 * (n[cont] - 1))
    denu2 <- v1[ncont]^2/(n[ncont]^2 * (n[ncont] - 1))
    denu <- N^2 * (denu1 + denu2)
    df <- numerator/denu
    pv <- c()
    for (h in 1:(a - 1)) {
        pv[h] <- min(2 * pt(T[h], df[h]), 2 - 2 * pt(T[h], df[h]))
    }
    adj.p <- p.adjust(pv, "hochberg")
    p.bonf <- p.adjust(pv, "bonferroni")
    p.holm <- p.adjust(pv, "holm")
    Rejected <- (adj.p <= alpha)
    vj <- which((1:a) != cont)
    vi <- rep(cont, a - 1)
    cmpid <- sapply(1:(a - 1), function(arg) {
        i <- vi[arg]
        j <- vj[arg]
        paste("F", "(", fl[j], ")", "-", "F", "(", fl[i], ")",
            sep = "")
    })
    result <- data.frame(Comparison = cmpid,
    Estimator = round(Cpd,4),
        df = round(df,4),
        Statistic = round(T,4),
        P.Raw = round(pv,4),
        P.Hochberg = round(adj.p,4),
        Rejected = Rejected,
        P.Bonf = round(p.bonf,4),
        P.Holm = round(p.holm,4))

    rownames(result) <- 1:(a - 1)
    output = list(Info = data.frame(Sample = fl,
    Size = n,
    Effect = round(pd,4),
    Variance = round(v1,4)),
        Analysis = result)
    if (!silent) {
        cat("\n","#----Xin Gao et al's (2008) Non-Parametric Multiple Test Procedure \n",
        "#----Type of Adjustment: Hochberg \n",
       "#----Level of significance", "=", alpha, "\n",
       "#----The procedure compares if the distribution functions F() are equal. The FWER is strongly controlled \n",
       "#---- This function uses pseudo ranks of the data!\n",
       "#----Reference: Gao, X. et al. (2008). Nonparametric Multiple Comparison Procedures for Unbalanced One-Way Factorial Designs. JSPI 138, 2574 - 2591. \n\n")

    }
    return(output)
}
