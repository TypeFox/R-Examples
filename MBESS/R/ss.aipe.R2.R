ss.aipe.R2 <- function (Population.R2 = NULL, conf.level = 0.95, width = NULL, Random.Predictors = TRUE, Random.Regressors, which.width = "Full", p = NULL, K, degree.of.certainty = NULL, assurance = NULL, certainty = NULL, verify.ss = FALSE, Tol = 1e-09, ...) 
{
if(!requireNamespace("gsl", quietly = TRUE)) stop("The package 'gsl' is needed; please install the package and try again.")


    if (!is.null(certainty) & is.null(degree.of.certainty) & 
        is.null(assurance)) 
        degree.of.certainty <- certainty
    if (is.null(assurance) && !is.null(degree.of.certainty) & 
        is.null(certainty)) 
        assurance <- degree.of.certainty
    if (!is.null(assurance) && is.null(degree.of.certainty) & 
        is.null(certainty)) 
        degree.of.certainty <- assurance
    if (!is.null(assurance) && !is.null(degree.of.certainty) && 
        assurance != degree.of.certainty) 
        stop("The arguments 'assurance' and 'degree.of.certainty' must have the same value.")
    if (!is.null(assurance) && !is.null(certainty) && assurance != 
        certainty) 
        stop("The arguments 'assurance' and 'certainty' must have the same value.")
    if (!is.null(degree.of.certainty) && !is.null(certainty) && 
        degree.of.certainty != certainty) 
        stop("The arguments 'degree.of.certainty' and 'certainty' must have the same value.")
    if (!missing(Random.Regressors)) 
        Random.Predictors <- Random.Regressors
    if (!missing(K)) 
        p <- K
    char.expand(which.width, c("Full", "Lower", "Upper"), nomatch = stop("Problems with 'which.width' specification. You must choose either 'Full', 'Lower', or 'Upper'.", 
        call. = FALSE))
    if (is.null(p)) 
        stop("You need to specify 'p', the number of predictors.")
    if (is.null(Population.R2)) 
        stop("You need to specify the population squared multiple correlation coefficient, 'Population.R2'.")
    if (!is.null(degree.of.certainty)) {
        if (degree.of.certainty <= 0.49 | degree.of.certainty > 
            1) 
            stop("The 'degree.of.certainty' must be between .50 and 1.")
    }
    Expected.R2 <- function(Population.R2, N, p) {
        Value <- 1 - ((N - p - 1)/(N - 1)) * (1 - Population.R2) * 
            gsl::hyperg_2F1(1, 1, 0.5 * (N + 1), Population.R2)
        Value <- max(0, Value)
        return(Value)
    }
    To.Find.Pop.R2.Given.Expectation <- function(E.R2.Low = Conf.Limit.Desired.Certainty.Lower, 
        E.R2.Up = Conf.Limit.Desired.Certainty.Upper, N = N, 
        p = p) {
        True.vals <- seq(0.001, 0.999, 0.001)
        Exp.vals <- rep(NA, length(True.vals))
        for (i in 1:length(True.vals)) {
            Exp.vals[i] <- Expected.R2(True.vals[i], N, p)
        }
        match(round(E.R2.Low, 2), round(Exp.vals, 3))
        return(c(mean(True.vals[match(round(E.R2.Low, 2), round(Exp.vals, 
            3))]), mean(True.vals[match(round(E.R2.Up, 2), round(Exp.vals, 
            3))])))
    }
    alpha.lower <- alpha.upper <- (1 - conf.level)/2
    N.0 <- p + 1 + p
    Continue <- TRUE
    while (Continue == TRUE) {
        CI.0 <- ci.R2(R2 = Expected.R2(Population.R2 = Population.R2, 
            N = N.0, p = p), conf.level = NULL, alpha.lower = alpha.lower, 
            alpha.upper = alpha.upper, N = N.0, p = p, Random.Predictors = FALSE)
        Continue <- sum(is.na(CI.0)) > 0
        N.0 <- N.0 + 1
    }
    if (which.width == "Full") {
        N.1 <- N.0 + 1
        CI.1 <- ci.R2(R2 = Expected.R2(Population.R2 = Population.R2, 
            N = N.1, p = p), conf.level = NULL, alpha.lower = alpha.lower, 
            alpha.upper = alpha.upper, N = N.1, p = p, Random.Predictors = FALSE)
        w.F <- CI.1$Upper.Conf.Limit.R2 - CI.1$Lower.Conf.Limit.R2
        Diff <- w.F - width
        while (Diff > Tol | is.na(Diff)) {
            N.1 <- N.1 + 1
            CI.1 <- ci.R2(R2 = Expected.R2(Population.R2 = Population.R2, 
                N = N.1, p = p), alpha.lower = alpha.lower, alpha.upper = alpha.upper, 
                N = N.1, p = p, Random.Predictors = FALSE)
            w.F <- CI.1$Upper.Conf.Limit.R2 - CI.1$Lower.Conf.Limit.R2
            Diff <- w.F - width
        }
        if (Random.Predictors == TRUE) {
            CI.1 <- ci.R2(R2 = Expected.R2(Population.R2 = Population.R2, 
                N = N.1, p = p), alpha.lower = alpha.lower, alpha.upper = alpha.upper, 
                N = N.1, p = p, Random.Predictors = TRUE)
            w.F <- CI.1$Upper.Conf.Limit.R2 - CI.1$Lower.Conf.Limit.R2
            Diff <- w.F - width
            while (Diff > Tol | is.na(Diff)) {
                N.1 <- N.1 + 1
                CI.1 <- ci.R2(R2 = Expected.R2(Population.R2 = Population.R2, 
                  N = N.1, p = p), alpha.lower = alpha.lower, 
                  alpha.upper = alpha.upper, N = N.1, p = p, 
                  Random.Predictors = TRUE)
                w.F <- CI.1$Upper.Conf.Limit.R2 - CI.1$Lower.Conf.Limit.R2
                Diff <- w.F - width
            }
        }
        Result.Full <- list(Required.Sample.Size = N.1)
        if (verify.ss == FALSE) {
            if (is.null(degree.of.certainty) == TRUE) {
                print("The approximate sample size is given below; you should consider using the additional")
                print("argument 'verify.ss=TRUE' to ensure the exact sample size value is obtained.")
                return(Result.Full)
            }
        }
        if (verify.ss == TRUE) {
            Result.Full <- list(Required.Sample.Size = verify.ss.aipe.R2(Population.R2 = Population.R2, 
                conf.level = conf.level, width = width, Random.Predictors = Random.Predictors, 
                which.width = "Full", p = p, n = Result.Full$Required.Sample.Size, 
                degree.of.certainty = degree.of.certainty, ...))
            if (is.null(degree.of.certainty) == TRUE) 
                return(Result.Full)
        }
        if (!is.null(degree.of.certainty)) {
            Conf.Limit.Desired.Certainty.Lower <- ci.R2(R2 = Expected.R2(Population.R2 = Population.R2, 
                N = N.1, p = p), alpha.lower = 1 - degree.of.certainty, 
                alpha.upper = 0, N = N.1, p = p, Random.Predictors = Random.Predictors)$Lower.Conf.Limit.R2
            Conf.Limit.Desired.Certainty.Upper <- ci.R2(R2 = Expected.R2(Population.R2 = Population.R2, 
                N = N.1, p = p), alpha.lower = 0, alpha.upper = 1 - 
                degree.of.certainty, N = N.1, p = p, Random.Predictors = Random.Predictors)$Upper.Conf.Limit.R2
            N.2 <- N.1
            CI.2 <- ci.R2(R2 = Conf.Limit.Desired.Certainty.Upper, 
                alpha.lower = alpha.lower, alpha.upper = alpha.upper, 
                N = N.2, p = p, Random.Predictors = Random.Predictors)
            w.F <- CI.2$Upper.Conf.Limit.R2 - CI.2$Lower.Conf.Limit.R2
            Diff <- w.F - width
            while (Diff > Tol | is.na(Diff)) {
                N.2 <- N.2 + 1
                Conf.Limit.Desired.Certainty.Upper.i <- ci.R2(R2 = Expected.R2(Population.R2 = Population.R2, 
                  N = N.2, p = p), alpha.lower = 0, alpha.upper = 1 - 
                  degree.of.certainty, N = N.2, p = p, Random.Predictors = Random.Predictors)$Upper.Conf.Limit.R2
                CI.2 <- ci.R2(R2 = Conf.Limit.Desired.Certainty.Upper.i, 
                  alpha.lower = alpha.lower, alpha.upper = alpha.upper, 
                  N = N.2, p = p, Random.Predictors = Random.Predictors)
                w.F <- CI.2$Upper.Conf.Limit.R2 - CI.2$Lower.Conf.Limit.R2
                Diff <- w.F - width
            }
            N.Upper.Conf.Lim <- N.2
            Ex.Width.Upper <- w.F
            N.3 <- N.1
            CI.3 <- ci.R2(R2 = Conf.Limit.Desired.Certainty.Lower, 
                alpha.lower = alpha.lower, alpha.upper = alpha.upper, 
                N = N.3, p = p, Random.Predictors = Random.Predictors)
            w.F <- CI.3$Upper.Conf.Limit.R2 - CI.3$Lower.Conf.Limit.R2
            Diff <- w.F - width
            while (Diff > Tol | is.na(Diff)) {
                N.3 <- N.3 + 1
                Conf.Limit.Desired.Certainty.Lower.i <- ci.R2(R2 = Expected.R2(Population.R2 = Population.R2, 
                  N = N.3, p = p), alpha.lower = 1 - degree.of.certainty, 
                  alpha.upper = 0, N = N.3, p = p, Random.Predictors = Random.Predictors)$Lower.Conf.Limit.R2
                CI.3 <- ci.R2(R2 = Conf.Limit.Desired.Certainty.Lower.i, 
                  alpha.lower = alpha.lower, alpha.upper = alpha.upper, 
                  N = N.3, p = p, Random.Predictors = Random.Predictors)
                w.F <- CI.3$Upper.Conf.Limit.R2 - CI.3$Lower.Conf.Limit.R2
                Diff <- w.F - width
            }
            N.Lower.Conf.Lim <- N.3
            Ex.Width.Lower <- w.F
            if (N.Upper.Conf.Lim >= N.Lower.Conf.Lim) {
                Result.Full <- list(Required.Sample.Size = N.Upper.Conf.Lim)
            }
            if (N.Lower.Conf.Lim > N.Upper.Conf.Lim) {
                Result.Full <- list(Required.Sample.Size = N.Lower.Conf.Lim)
            }
            CI.WIDTH.R2 <- function(R2, conf.level, N, p, Random.Predictors) {
                Lims <- ci.R2(R2 = Expected.R2(Population.R2 = R2, 
                  N, p), conf.level = conf.level, N = N, p = p, 
                  Random.Predictors = Random.Predictors)
                Lims$Upper - Lims$Lower
            }
            Pop.Values <- To.Find.Pop.R2.Given.Expectation(E.R2.Low = Conf.Limit.Desired.Certainty.Lower, 
                E.R2.Up = Conf.Limit.Desired.Certainty.Upper, 
                N = N.1, p = p)
            Optimize.Result <- optimize(f = CI.WIDTH.R2, interval = c(max(Pop.Values[1] * 
                0.98, 1e-06), min(Pop.Values[2] * 1.02, 0.999999)), 
                maximum = TRUE, tol = .Machine$double.eps^0.5, 
                N = N.1, p = p, conf.level = conf.level, Random.Predictors = Random.Predictors)$maximum
            Optimize.Result <- Expected.R2(Population.R2 = Optimize.Result, 
                N.1, p)
            N.Op <- NULL
            N.Op.Up <- NULL
            N.Op.Low <- NULL
            if ((round(Conf.Limit.Desired.Certainty.Lower, 4) < 
                round(Optimize.Result, 4)) & (round(Optimize.Result, 
                4) < round(Conf.Limit.Desired.Certainty.Upper, 
                4))) {
                if (Random.Predictors == TRUE) {
                  Conf.Limit.Desired.Certainty.Lower <- ci.R2(R2 = Optimize.Result, 
                    alpha.lower = 1 - degree.of.certainty, alpha.upper = 0, 
                    N = N.1, p = p, Random.Predictors = Random.Predictors)$Lower.Conf.Limit.R2
                  Conf.Limit.Desired.Certainty.Upper <- ci.R2(R2 = Optimize.Result, 
                    alpha.lower = 0, alpha.upper = 1 - degree.of.certainty, 
                    N = N.1, p = p, Random.Predictors = Random.Predictors)$Upper.Conf.Limit.R2
                  if (Conf.Limit.Desired.Certainty.Lower > Optimize.Result) 
                    Conf.Limit.Desired.Certainty.Lower <- Optimize.Result
                  N.Op.Low <- N.1
                  CI.Op.Low <- ci.R2(R2 = Conf.Limit.Desired.Certainty.Lower, 
                    alpha.lower = alpha.lower, alpha.upper = alpha.upper, 
                    N = N.Op.Low, p = p, Random.Predictors = Random.Predictors)
                  w.F <- CI.Op.Low$Upper.Conf.Limit.R2 - CI.Op.Low$Lower.Conf.Limit.R2
                  Diff <- w.F - width
                  while (Diff > Tol | is.na(Diff)) {
                    N.Op.Low <- N.Op.Low + 1
                    CI.Op.Low <- ci.R2(R2 = Conf.Limit.Desired.Certainty.Lower, 
                      alpha.lower = alpha.lower, alpha.upper = alpha.upper, 
                      N = N.Op.Low, p = p, Random.Predictors = Random.Predictors)
                    w.F <- CI.Op.Low$Upper.Conf.Limit.R2 - CI.Op.Low$Lower.Conf.Limit.R2
                    Diff <- w.F - width
                  }
                  w.Op.Width.Low <- w.F
                  if (Conf.Limit.Desired.Certainty.Upper < Optimize.Result) 
                    Conf.Limit.Desired.Certainty.Upper <- Optimize.Result
                  N.Op.Up <- N.1
                  CI.Op.Up <- ci.R2(R2 = Conf.Limit.Desired.Certainty.Upper, 
                    alpha.lower = alpha.lower, alpha.upper = alpha.upper, 
                    N = N.Op.Up, p = p, Random.Predictors = Random.Predictors)
                  w.F <- CI.Op.Low$Upper.Conf.Limit.R2 - CI.Op.Low$Lower.Conf.Limit.R2
                  Diff <- w.F - width
                  while (Diff > Tol | is.na(Diff)) {
                    N.Op.Up <- N.Op.Up + 1
                    CI.Op.Up <- ci.R2(R2 = Conf.Limit.Desired.Certainty.Upper, 
                      alpha.lower = alpha.lower, alpha.upper = alpha.upper, 
                      N = N.Op.Up, p = p, Random.Predictors = Random.Predictors)
                    w.F <- CI.Op.Up$Upper.Conf.Limit.R2 - CI.Op.Up$Lower.Conf.Limit.R2
                    Diff <- w.F - width
                  }
                  w.Op.Width.Up <- w.F
                  if (N.1 == N.Lower.Conf.Lim & N.1 == N.Upper.Conf.Lim & 
                    N.1 == N.Op.Up & N.1 == N.Op.Low) {
                    N.Op <- N.1
                    CI.Op <- ci.R2(R2 = Optimize.Result, alpha.lower = alpha.lower, 
                      alpha.upper = alpha.upper, N = N.Op, p = p, 
                      Random.Predictors = Random.Predictors)
                    w.F <- CI.Op$Upper.Conf.Limit.R2 - CI.Op$Lower.Conf.Limit.R2
                    Diff <- w.F - width
                    while (Diff > Tol | is.na(Diff)) {
                      N.Op <- N.Op + 1
                      CI.Op <- ci.R2(R2 = Optimize.Result, alpha.lower = alpha.lower, 
                        alpha.upper = alpha.upper, N = N.Op, 
                        p = p, Random.Predictors = Random.Predictors)
                      w.F <- CI.Op$Upper.Conf.Limit.R2 - CI.Op$Lower.Conf.Limit.R2
                      Diff <- w.F - width
                    }
                    w.Op <- w.F
                  }
                  Necessary.N <- max(N.1, N.Lower.Conf.Lim, N.Upper.Conf.Lim, 
                    N.Op.Low, N.Op.Up, N.Op, na.rm = TRUE)
                  For.Ex.Width <- ci.R2(R2 = Expected.R2(Population.R2 = Population.R2, 
                    N = Necessary.N, p = p), alpha.lower = alpha.lower, 
                    alpha.upper = alpha.upper, N = Necessary.N, 
                    p = p, Random.Predictors = Random.Predictors)
                  Result.Full <- list(Required.Sample.Size = Necessary.N)
                }
                if (Random.Predictors == FALSE) {
                  Prob.F.Max.Width <- pf(q = Rsquare2F(R2 = Optimize.Result, 
                    p = p, N = N.1), df1 = p, df2 = N.1 - p - 
                    1, ncp = Rsquare2Lambda(R2 = Population.R2, 
                    N = N.1))
                  low.lim.F <- qf(p = max(0, (Prob.F.Max.Width - 
                    (1 - degree.of.certainty)/2)), df1 = p, df2 = N.1 - 
                    p - 1, ncp = Rsquare2Lambda(R2 = Population.R2, 
                    N = N.1))
                  up.lim.F <- qf(p = min((Prob.F.Max.Width + 
                    (1 - degree.of.certainty)/2), 1), df1 = p, 
                    df2 = N.1 - p - 1, ncp = Rsquare2Lambda(R2 = Population.R2, 
                      N = N.1))
                  low.R2 <- F2Rsquare(F.value = low.lim.F, df.1 = p, 
                    df.2 = N.1 - p - 1)
                  up.R2 <- F2Rsquare(F.value = up.lim.F, df.1 = p, 
                    df.2 = N.1 - p - 1)
                  n.1 <- N.1 - 1
                  if (up.R2 > 0 & up.R2 < 1) {
                    N.4 <- n.1 + 1
                    CI.4 <- ci.R2(R2 = Expected.R2(Population.R2 = up.R2, 
                      N = N.4, p = p), alpha.lower = alpha.lower, 
                      alpha.upper = alpha.upper, N = N.4, p = p, 
                      Random.Predictors = Random.Predictors)
                    w.F <- CI.4$Upper.Conf.Limit.R2 - CI.4$Lower.Conf.Limit.R2
                    Diff <- w.F - width
                    while (Diff > Tol | is.na(Diff)) {
                      N.4 <- N.4 + 1
                      CI.4 <- ci.R2(R2 = Expected.R2(Population.R2 = up.R2, 
                        N = N.4, p = p), alpha.lower = alpha.lower, 
                        alpha.upper = alpha.upper, N = N.4, p = p, 
                        Random.Predictors = Random.Predictors)
                      w.F <- CI.4$Upper.Conf.Limit.R2 - CI.4$Lower.Conf.Limit.R2
                      Diff <- w.F - width
                    }
                    n.up.R2 <- N.4
                    Ex.Width.Upper.after.optim <- w.F
                  }
                  if (low.R2 > 0 & low.R2 < 1) {
                    N.5 <- n.1 + 1
                    CI.5 <- ci.R2(R2 = Expected.R2(Population.R2 = low.R2, 
                      N = N.5, p = p), alpha.lower = alpha.lower, 
                      alpha.upper = alpha.upper, N = N.5, p = p, 
                      Random.Predictors = Random.Predictors)
                    w.F <- CI.5$Upper.Conf.Limit.R2 - CI.5$Lower.Conf.Limit.R2
                    Diff <- w.F - width
                    while (Diff > Tol | is.na(Diff)) {
                      N.5 <- N.5 + 1
                      CI.5 <- ci.R2(R2 = Expected.R2(Population.R2 = low.R2, 
                        N = N.5, p = p), alpha.lower = alpha.lower, 
                        alpha.upper = alpha.upper, N = N.5, p = p, 
                        Random.Predictors = Random.Predictors)
                      w.F <- CI.5$Upper.Conf.Limit.R2 - CI.5$Lower.Conf.Limit.R2
                      Diff <- w.F - width
                    }
                    n.low.R2 <- N.5
                    Ex.Width.Lower.after.optim <- w.F
                  }
                  if (low.R2 == 0 | low.R2 == 1) 
                    n.low.R2 <- N.1
                  if (up.R2 == 0 | up.R2 == 1) 
                    n.up.R2 <- N.1
                  if ((n.low.R2 >= n.up.R2) & (n.low.R2 > Result.Full$Required.Sample.Size)) {
                    Result.Full <- list(Required.Sample.Size = n.low.R2)
                  }
                  if ((n.up.R2 > n.low.R2) & (n.up.R2 > Result.Full$Required.Sample.Size)) {
                    Result.Full <- list(Required.Sample.Size = n.up.R2)
                  }
                }
            }
            if (verify.ss == FALSE) {
                print("The approximate sample size is given below; you should consider using the additional")
                print("argument 'verify.ss=TRUE' to ensure the exact sample size value is obtained.")
                return(Result.Full)
            }
            if (verify.ss == TRUE) {
                Result.Full <- list(Required.Sample.Size = verify.ss.aipe.R2(Population.R2 = Population.R2, 
                  conf.level = conf.level, width = width, Random.Predictors = Random.Predictors, 
                  which.width = "Full", p = p, n = Result.Full$Required.Sample.Size, 
                  degree.of.certainty = degree.of.certainty, 
                  ...))
                return(Result.Full)
            }
        }
    }
    if (which.width == "Lower") {
        if (verify.ss == TRUE) 
            stop("verify.ss can only be used with 'width=FULL'.")
        N.1 <- N.0 + 1
        CI.1 <- ci.R2(R2 = Expected.R2(Population.R2 = Population.R2, 
            N = N.1, p = p), alpha.lower = alpha.lower, alpha.upper = alpha.upper, 
            N = N.1, p = p, Random.Predictors = FALSE)
        w.L <- Expected.R2(Population.R2 = Population.R2, N = N.1, 
            p = p) - CI.1$Lower.Conf.Limit.R2
        Diff <- w.L - width
        while (Diff > Tol | is.na(Diff)) {
            N.1 <- N.1 + 1
            CI.1 <- ci.R2(R2 = Expected.R2(Population.R2 = Population.R2, 
                N = N.1, p = p), alpha.lower = alpha.lower, alpha.upper = alpha.upper, 
                N = N.1, p = p, Random.Predictors = FALSE)
            w.L <- Expected.R2(Population.R2 = Population.R2, 
                N = N.1, p = p) - CI.1$Lower.Conf.Limit.R2
            Diff <- w.L - width
        }
        if (Random.Predictors == TRUE) {
            CI.1 <- ci.R2(R2 = Expected.R2(Population.R2 = Population.R2, 
                N = N.1, p = p), alpha.lower = alpha.lower, alpha.upper = alpha.upper, 
                N = N.1, p = p, Random.Predictors = TRUE)
            w.L <- Expected.R2(Population.R2 = Population.R2, 
                N = N.1, p = p) - CI.1$Lower.Conf.Limit.R2
            Diff <- w.L - width
            while (Diff > Tol | is.na(Diff)) {
                N.1 <- N.1 + 1
                CI.1 <- ci.R2(R2 = Expected.R2(Population.R2 = Population.R2, 
                  N = N.1, p = p), alpha.lower = alpha.lower, 
                  alpha.upper = alpha.upper, N = N.1, p = p, 
                  Random.Predictors = TRUE)
                w.L <- Expected.R2(Population.R2 = Population.R2, 
                  N = N.1, p = p) - CI.1$Lower.Conf.Limit.R2
                Diff <- w.L - width
            }
        }
        if (is.null(degree.of.certainty)) {
            Result.Low <- list(Required.Sample.Size = N.1)
            print("This sample size should be regarded as VERY APPROXIMATE.")
            print("The methods used have only been throughly evaluated for Full widths.")
            return(Result.Low)
        }
        if (!is.null(degree.of.certainty)) {
            Conf.Limit.Desired.Certainty <- ci.R2(R2 = Expected.R2(Population.R2 = Population.R2, 
                N = N.1, p = p), alpha.lower = 0, alpha.upper = 1 - 
                degree.of.certainty, N = N.1, p = p, Random.Predictors = Random.Predictors)$Upper.Conf.Limit.R2
            N.2 <- N.1
            CI.2 <- ci.R2(R2 = Expected.R2(Population.R2 = Conf.Limit.Desired.Certainty, 
                N = N.2, p = p), alpha.lower = alpha.lower, alpha.upper = alpha.upper, 
                N = N.2, p = p, Random.Predictors = Random.Predictors)
            w.L <- Conf.Limit.Desired.Certainty - CI.2$Lower.Conf.Limit.R2
            Diff <- w.L - width
            while (Diff > Tol) {
                N.2 <- N.2 + 1
                CI.2 <- ci.R2(R2 = Expected.R2(Population.R2 = Conf.Limit.Desired.Certainty, 
                  N = N.2, p = p), alpha.lower = alpha.lower, 
                  alpha.upper = alpha.upper, N = N.2, p = p, 
                  Random.Predictors = Random.Predictors)
                w.L <- Conf.Limit.Desired.Certainty - CI.2$Lower.Conf.Limit.R2
                Diff <- w.L - width
            }
            w.L.M.2 <- Expected.R2(Population.R2 = Population.R2, 
                N = N.2, p = p) - ci.R2(R2 = Expected.R2(Population.R2 = Population.R2, 
                N = N.2, p = p), alpha.lower = alpha.lower, alpha.upper = alpha.upper, 
                N = N.2, p = p, Random.Predictors = Random.Predictors)$Lower.Conf.Limit.R2
            Result.Low.2 <- list(Required.Sample.Size = N.2)
            Conf.Limit.Desired.Certainty <- ci.R2(R2 = Expected.R2(Population.R2 = Population.R2, 
                N = N.1, p = p), alpha.lower = 1 - degree.of.certainty, 
                alpha.upper = 0, N = N.1, p = p, Random.Predictors = Random.Predictors)$Lower.Conf.Limit.R2
            N.3 <- N.1
            CI.2 <- ci.R2(R2 = Expected.R2(Population.R2 = Conf.Limit.Desired.Certainty, 
                N = N.3, p = p), alpha.lower = alpha.lower, alpha.upper = alpha.upper, 
                N = N.3, p = p, Random.Predictors = Random.Predictors)
            w.L <- Conf.Limit.Desired.Certainty - CI.2$Lower.Conf.Limit.R2
            Diff <- w.L - width
            while (Diff > Tol) {
                N.3 <- N.3 + 1
                CI.2 <- ci.R2(R2 = Expected.R2(Population.R2 = Conf.Limit.Desired.Certainty, 
                  N = N.3, p = p), alpha.lower = alpha.lower, 
                  alpha.upper = alpha.upper, N = N.3, p = p, 
                  Random.Predictors = Random.Predictors)
                w.L <- Conf.Limit.Desired.Certainty - CI.2$Lower.Conf.Limit.R2
                Diff <- w.L - width
            }
            w.L.M.3 <- Expected.R2(Population.R2 = Population.R2, 
                N = N.3, p = p) - ci.R2(R2 = Expected.R2(Population.R2 = Population.R2, 
                N = N.3, p = p), alpha.lower = alpha.lower, alpha.upper = alpha.upper, 
                N = N.3, p = p, Random.Predictors = Random.Predictors)$Lower.Conf.Limit.R2
            Result.Low.3 <- list(Required.Sample.Size = N.3)
            if (N.2 >= N.3) 
                Result.Low <- Result.Low.2
            if (N.3 > N.2) 
                Result.Low <- Result.Low.3
            print("This sample size should be regarded as VERY APPROXIMATE.")
            print("The methods used have only been throughly evaluated for Full widths.")
            return(Result.Low)
        }
    }
    if (which.width == "Upper") {
        if (verify.ss == TRUE) 
            stop("verify.ss can only be used with 'width=FULL'.")
        N.1 <- N.0 + 1
        CI.1 <- ci.R2(R2 = Expected.R2(Population.R2 = Population.R2, 
            N = N.1, p = p), alpha.lower = alpha.lower, alpha.upper = alpha.upper, 
            N = N.1, p = p, Random.Predictors = FALSE)
        w.U <- CI.1$Upper.Conf.Limit.R2 - Expected.R2(Population.R2 = Population.R2, 
            N = N.1, p = p)
        Diff <- w.U - width
        while (Diff > Tol | is.na(Diff)) {
            N.1 <- N.1 + 1
            CI.1 <- ci.R2(R2 = Expected.R2(Population.R2 = Population.R2, 
                N = N.1, p = p), alpha.lower = alpha.lower, alpha.upper = alpha.upper, 
                N = N.1, p = p, Random.Predictors = FALSE)
            w.U <- CI.1$Upper.Conf.Limit.R2 - Expected.R2(Population.R2 = Population.R2, 
                N = N.1, p = p)
            Diff <- w.U - width
        }
        if (Random.Predictors == TRUE) {
            CI.1 <- ci.R2(R2 = Expected.R2(Population.R2 = Population.R2, 
                N = N.1, p = p), alpha.lower = alpha.lower, alpha.upper = alpha.upper, 
                N = N.1, p = p, Random.Predictors = TRUE)
            w.U <- CI.1$Upper.Conf.Limit.R2 - Expected.R2(Population.R2 = Population.R2, 
                N = N.1, p = p)
            Diff <- w.U - width
            while (Diff > Tol | is.na(Diff)) {
                N.1 <- N.1 + 1
                CI.1 <- ci.R2(R2 = Expected.R2(Population.R2 = Population.R2, 
                  N = N.1, p = p), alpha.lower = alpha.lower, 
                  alpha.upper = alpha.upper, N = N.1, p = p, 
                  Random.Predictors = TRUE)
                w.U <- CI.1$Upper.Conf.Limit.R2 - Expected.R2(Population.R2 = Population.R2, 
                  N = N.1, p = p)
                Diff <- w.U - width
            }
        }
        if (is.null(degree.of.certainty)) {
            Result.Up <- list(Required.Sample.Size = N.1)
            print("This sample size should be regarded as VERY APPROXIMATE.")
            print("The methods used have only been throughly evaluated for Full widths.")
            return(Result.Up)
        }
        if (!is.null(degree.of.certainty)) {
            Conf.Limit.Desired.Certainty <- ci.R2(R2 = Expected.R2(Population.R2 = Population.R2, 
                N = N.1, p = p), alpha.lower = 0, alpha.upper = 1 - 
                degree.of.certainty, N = N.1, p = p, Random.Predictors = Random.Predictors)$Upper.Conf.Limit.R2
            N.2 <- N.1
            CI.2 <- ci.R2(R2 = Expected.R2(Population.R2 = Conf.Limit.Desired.Certainty, 
                N = N.2, p = p), alpha.lower = alpha.lower, alpha.upper = alpha.upper, 
                N = N.2, p = p, Random.Predictors = Random.Predictors)
            w.U <- CI.2$Upper.Conf.Limit.R2 - Conf.Limit.Desired.Certainty
            Diff <- w.U - width
            while (Diff > Tol) {
                N.2 <- N.2 + 1
                CI.2 <- ci.R2(R2 = Expected.R2(Population.R2 = Conf.Limit.Desired.Certainty, 
                  N = N.2, p = p), alpha.lower = alpha.lower, 
                  alpha.upper = alpha.upper, N = N.2, p = p, 
                  Random.Predictors = Random.Predictors)
                w.U <- CI.2$Upper.Conf.Limit.R2 - Conf.Limit.Desired.Certainty
                Diff <- w.U - width
            }
            w.U.M.2 <- ci.R2(R2 = Expected.R2(Population.R2 = Population.R2, 
                N = N.2, p = p, Random.Predictors = Random.Predictors), 
                alpha.lower = alpha.lower, alpha.upper = alpha.upper, 
                N = N.2, p = p)$Upper.Conf.Limit.R2 - Expected.R2(Population.R2 = Population.R2, 
                N = N.2, p = p)
            Result.Up.2 <- list(Required.Sample.Size = N.2)
            Conf.Limit.Desired.Certainty <- ci.R2(R2 = Expected.R2(Population.R2 = Population.R2, 
                N = N.1, p = p), alpha.lower = 1 - degree.of.certainty, 
                alpha.upper = 0, N = N.1, p = p, Random.Predictors = Random.Predictors)$Lower.Conf.Limit.R2
            N.3 <- N.1
            CI.2 <- ci.R2(R2 = Expected.R2(Population.R2 = Conf.Limit.Desired.Certainty, 
                N = N.3, p = p), alpha.lower = alpha.lower, alpha.upper = alpha.upper, 
                N = N.3, p = p, Random.Predictors = Random.Predictors)
            w.U <- CI.2$Upper.Conf.Limit.R2 - Conf.Limit.Desired.Certainty
            Diff <- w.U - width
            while (Diff > Tol) {
                N.3 <- N.3 + 1
                CI.2 <- ci.R2(R2 = Expected.R2(Population.R2 = Conf.Limit.Desired.Certainty, 
                  N = N.3, p = p), alpha.lower = alpha.lower, 
                  alpha.upper = alpha.upper, N = N.3, p = p, 
                  Random.Predictors = Random.Predictors)
                w.U <- CI.2$Upper.Conf.Limit.R2 - Conf.Limit.Desired.Certainty
                Diff <- w.U - width
            }
            w.U.M.3 <- ci.R2(R2 = Expected.R2(Population.R2 = Population.R2, 
                N = N.3, p = p), alpha.lower = alpha.lower, alpha.upper = alpha.upper, 
                N = N.3, p = p)$Upper.Conf.Limit.R2 - Expected.R2(Population.R2 = Population.R2, 
                N = N.3, p = p, Random.Predictors = Random.Predictors)
            Result.Up.3 <- list(Required.Sample.Size = N.3)
            if (N.2 >= N.3) 
                Result.Up <- Result.Up.2
            if (N.3 > N.2) 
                Result.Up <- Result.Up.3
            print("This sample size should be regarded as VERY APPROXIMATE.")
            print("The methods used have only been throughly evaluated for Full widths.")
            return(Result.Up)
        }
    }
}
