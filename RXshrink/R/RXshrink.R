"plot.RXlarlso" <-
function (x, trace = "all", trkey = FALSE, ...) 
{
    mcal <- x$sext[, 2]
    mcalp <- rep(mcal, times = x$p)
    if (trace != "coef" && trace != "rmse" && trace != "exev" && 
        trace != "infd" && trace != "spat" && trace != "seq") 
        trace <- "all"
    opar <- par(no.readonly = TRUE)
    on.exit(par(opar))
    if (trace == "all")
        par(mfrow=c(3,2))
    else
        par(mfrow=c(1,1))
    if (trace == "all" || trace == "seq" || trace == "coef") {
        plot(mcalp, x$coef, ann = FALSE, type = "n")
        abline(h = 0, col = gray(0.9))
        for (i in 1:x$p) lines(mcal, x$coef[, i], col = i, lty = i, 
            lwd = 2)
        title(main = paste("COEFFICIENT TRACE:", x$lars$type), 
            xlab = "m = Multicollinearity Allowance", ylab = "Fitted Coefficients")
        if( trkey )
            legend("bottom",all.vars(x$form)[2:(x$p+1)], col=1:(x$p), lty=1:(x$p), lwd=2)
    }
    if (trace == "seq") {
        cat("\nPress the Enter key to view the RMSE trace...")
        scan()
    }
    if (trace == "all" || trace == "seq" || trace == "rmse") {
        plot(mcalp, x$rmse, ann = FALSE, type = "n")
        abline(h = 0, col = gray(0.9))
        for (i in 1:x$p) lines(mcal, x$rmse[, i], col = i, lty = i, 
            lwd = 2)
        title(main = paste("RELATIVE MEAN SQ. ERROR:", x$lars$type), 
            xlab = "m = Multicollinearity Allowance", ylab = "Scaled MSE Risk")
        if( trkey )
            legend("bottom",all.vars(x$form)[2:(x$p+1)], col=1:(x$p), lty=1:(x$p), lwd=2)
    }
    if (trace == "seq") {
        cat("\nPress the Enter key to view the EXEV trace...")
        scan()
    }
    if (trace == "all" || trace == "seq" || trace == "exev") {
        plot(mcalp, x$exev, ann = FALSE, type = "n")
        abline(h = 0, col = gray(0.9))
        for (i in 1:x$p) lines(mcal, x$exev[, i], col = i, lty = i, 
            lwd = 2)
        title(main = paste("EXCESS EIGENVALUES:", x$lars$type), 
            xlab = "m = Multicollinearity Allowance", ylab = "Least Squares minus larlso")
        if( trkey )
            legend("bottom",paste("Component",1:(x$p)), col=1:(x$p), lty=1:(x$p), lwd=2)
    }
    if (trace == "seq") {
        cat("\nPress the Enter key to view the INFD trace...")
        scan()
    }
    if (trace == "all" || trace == "seq" || trace == "infd") {
        plot(mcalp, x$infd, ann = FALSE, type = "n")
        abline(h = 0, col = gray(0.9))
        for (i in 1:x$p) lines(mcal, x$infd[, i], col = i, lty = i, 
            lwd = 2)
        title(main = paste("INFERIOR DIRECTION:", x$lars$type), 
            xlab = "m = Multicollinearity Allowance", ylab = "Direction Cosines")
        if( trkey )
            legend("bottom",all.vars(x$form)[2:(x$p+1)], col=1:(x$p), lty=1:(x$p), lwd=2)
    }
    if (trace == "seq") {
        cat("\nPress the Enter key to view the SPAT trace...")
        scan()
    }
    if (trace == "all" || trace == "seq" || trace == "spat") {
        plot(mcalp, x$spat, ann = FALSE, type = "n")
        abline(h = 0, col = gray(0.9))
        for (i in 1:x$p) lines(mcal, x$spat[, i], col = i, lty = i, 
            lwd = 2)
        title(main = paste("SHRINKAGE PATTERN:", x$lars$type), 
            xlab = "m = Multicollinearity Allowance", ylab = "larlso Delta Factors")
        if( trkey )
            legend("bottom",paste("Component",1:(x$p)), col=1:(x$p), lty=1:(x$p), lwd=2)
    }
}

"plot.RXridge" <-
function (x, trace = "all", trkey = FALSE, ...) 
{
    mcal <- x$sext[, 3]
    mcalp <- rep(mcal, times = x$p)
    if (trace != "coef" && trace != "rmse" && trace != "exev" && 
        trace != "infd" && trace != "spat" && trace != "seq")
        trace <- "all"
    opar <- par(no.readonly = TRUE)
    on.exit(par(opar))
    if (trace == "all")
        par(mfrow=c(3,2))
    else
        par(mfrow=c(1,1))
    if (trace == "all" || trace == "seq" || trace == "coef") {
        plot(mcalp, x$coef, ann = FALSE, type = "n")
        abline(h = 0, col = gray(0.9))
        for (i in 1:x$p) lines(mcal, x$coef[, i], col = i, lty = i, 
            lwd = 2)
        title(main = paste("COEFFICIENT TRACE: Q-shape =", x$qp), 
            xlab = "m = Multicollinearity Allowance", ylab = "Fitted Coefficients")
        if( trkey )
            legend("bottom",all.vars(x$form)[2:(x$p+1)], col=1:(x$p), lty=1:(x$p), lwd=2)
    }
    if (trace == "seq") {
        cat("\nPress the Enter key to view the RMSE trace...")
        scan()
    }
    if (trace == "all" || trace == "seq" || trace == "rmse") {
        plot(mcalp, x$rmse, ann = FALSE, type = "n")
        abline(h = 0, col = gray(0.9))
        for (i in 1:x$p) lines(mcal, x$rmse[, i], col = i, lty = i, 
            lwd = 2)
        title(main = paste("RELATIVE MEAN SQ. ERROR: Q-shape =", 
            x$qp), xlab = "m = Multicollinearity Allowance", 
            ylab = "Scaled MSE Risk")
        if( trkey )
            legend("bottom",all.vars(x$form)[2:(x$p+1)], col=1:(x$p), lty=1:(x$p), lwd=2)
    }
    if (trace == "seq") {
        cat("\nPress the Enter key to view the EXEV trace...")
        scan()
    }
    if (trace == "all" || trace == "seq" || trace == "exev") {
        plot(mcalp, x$exev, ann = FALSE, type = "n")
        abline(h = 0, col = gray(0.9))
        for (i in 1:x$p) lines(mcal, x$exev[, i], col = i, lty = i, 
            lwd = 2)
        title(main = paste("EXCESS EIGENVALUES: Q-shape =", x$qp), 
            xlab = "m = Multicollinearity Allowance", ylab = "Least Squares minus Ridge")
        if( trkey )
            legend("bottom",paste("Component",1:(x$p)), col=1:(x$p), lty=1:(x$p), lwd=2)
    }
    if (trace == "seq") {
        cat("\nPress the Enter key to view the INFD trace...")
        scan()
    }
    if (trace == "all" || trace == "seq" || trace == "infd") {
        plot(mcalp, x$infd, ann = FALSE, type = "n")
        abline(h = 0, col = gray(0.9))
        for (i in 1:x$p) lines(mcal, x$infd[, i], col = i, lty = i, 
            lwd = 2)
        title(main = paste("INFERIOR DIRECTION: Q-shape =", x$qp), 
            xlab = "m = Multicollinearity Allowance", ylab = "Direction Cosines")
        if( trkey )
            legend("bottom",all.vars(x$form)[2:(x$p+1)], col=1:(x$p), lty=1:(x$p), lwd=2)
    }
    if (trace == "seq") {
        cat("\nPress the Enter key to view the SPAT trace...")
        scan()
    }
    if (trace == "all" || trace == "seq" || trace == "spat") {
        plot(mcalp, x$spat, ann = FALSE, type = "n")
        abline(h = 0, col = gray(0.9))
        for (i in 1:x$p) lines(mcal, x$spat[, i], col = i, lty = i, 
            lwd = 2)
        title(main = paste("SHRINKAGE PATTERN: Q-shape =", x$qp), 
            xlab = "m = Multicollinearity Allowance", ylab = "Ridge Delta Factors")
        if( trkey )
            legend("bottom",paste("Component",1:(x$p)), col=1:(x$p), lty=1:(x$p), lwd=2)
    }
}

"plot.RXtrisk" <-
function (x, trace = "all", trkey = FALSE, ...) 
{
    mcal <- x$sext[, 3]
    mcalp <- rep(mcal, times = x$p)
    if (trace != "coef" && trace != "rmse" && trace != "exev" && 
        trace != "infd" && trace != "spat" && trace != "seq")
        trace <- "all"
    opar <- par(no.readonly = TRUE)
    on.exit(par(opar))
    if (trace == "all")
        par(mfrow=c(3,2))
    else
        par(mfrow=c(1,1))
    if (trace == "all" || trace == "seq" || trace == "coef") {
        plot(mcalp, x$coef, ann = FALSE, type = "n")
        abline(h = 0, col = gray(0.9))
        for (i in 1:x$p) lines(mcal, x$coef[, i], col = i, lty = i, 
            lwd = 2)
        title(main = paste("EXPECTED COEFFICIENTS: Q-shape =", 
            x$qp), xlab = "m = Multicollinearity Allowance", 
            ylab = "Expected Coefficients")
        if( trkey )
            legend("bottom",all.vars(x$form)[2:(x$p+1)], col=1:(x$p), lty=1:(x$p), lwd=2)
    }
    if (trace == "seq") {
        cat("\nPress the Enter key to view the RMSE trace...")
        scan()
    }
    if (trace == "all" || trace == "seq" || trace == "rmse") {
        plot(mcalp, x$rmse, ann = FALSE, type = "n")
        abline(h = 0, col = gray(0.9))
        for (i in 1:x$p) lines(mcal, x$rmse[, i], col = i, lty = i, 
            lwd = 2)
        title(main = paste("RELATIVE TRUE MeanSqError: Q-shape =", 
            x$qp), xlab = "m = Multicollinearity Allowance", 
            ylab = "Scaled True MSE Risk")
        if( trkey )
            legend("bottom",all.vars(x$form)[2:(x$p+1)], col=1:(x$p), lty=1:(x$p), lwd=2)
    }
    if (trace == "seq") {
        cat("\nPress the Enter key to view the EXEV trace...")
        scan()
    }
    if (trace == "all" || trace == "seq" || trace == "exev") {
        plot(mcalp, x$exev, ann = FALSE, type = "n")
        abline(h = 0, col = gray(0.9))
        for (i in 1:x$p) lines(mcal, x$exev[, i], col = i, lty = i, 
            lwd = 2)
        title(main = paste("TRUE EXCESS EIGENVALUES: Q-shape =", 
            x$qp), xlab = "m = Multicollinearity Allowance", 
            ylab = "Least Squares minus Ridge")
        if( trkey )
            legend("bottom",paste("Component",1:(x$p)), col=1:(x$p), lty=1:(x$p), lwd=2)
    }
    if (trace == "seq") {
        cat("\nPress the Enter key to view the INFD trace...")
        scan()
    }
    if (trace == "all" || trace == "seq" || trace == "infd") {
        plot(mcalp, x$infd, ann = FALSE, type = "n")
        abline(h = 0, col = gray(0.9))
        for (i in 1:x$p) lines(mcal, x$infd[, i], col = i, lty = i, 
            lwd = 2)
        title(main = paste("TRUE INFERIOR DIRECTION: Q-shape =", 
            x$qp), xlab = "m = Multicollinearity Allowance", 
            ylab = "Inferior Direction Cosines")
        if( trkey )
            legend("bottom",all.vars(x$form)[2:(x$p+1)], col=1:(x$p), lty=1:(x$p), lwd=2)
    }
    if (trace == "seq") {
        cat("\nPress the Enter key to view the SPAT trace...")
        scan()
    }
    if (trace == "all" || trace == "seq" || trace == "spat") {
        plot(mcalp, x$spat, ann = FALSE, type = "n")
        abline(h = 0, col = gray(0.9))
        for (i in 1:x$p) lines(mcal, x$spat[, i], col = i, lty = i, 
            lwd = 2)
        title(main = paste("SHRINKAGE PATTERN: Q-shape =", x$qp), 
            xlab = "m = Multicollinearity Allowance", ylab = "Ridge Delta Factors")
        if( trkey )
            legend("bottom",paste("Component",1:(x$p)), col=1:(x$p), lty=1:(x$p), lwd=2)
    }
}

"plot.RXtsimu" <-
function (x, trace = "all", trkey = FALSE, ...) 
{
    mcal <- x$sext[, 3]
    mcalp <- rep(mcal, times = x$p)
    if (trace != "coef" && trace != "rsel" &&
        trace != "spat" && trace != "seq")
        trace <- "all"
    opar <- par(no.readonly = TRUE)
    on.exit(par(opar))
    if (trace == "all")
        par(mfrow=c(2,2))
    else
        par(mfrow=c(1,1))
    if (trace == "all" || trace == "seq" || trace == "coef") {
        plot(mcalp, x$coef, ann = FALSE, type = "n")
        abline(h = 0, col = gray(0.9))
        for (i in 1:x$p) lines(mcal, x$coef[, i], col = i, lty = i, 
            lwd = 2)
        title(main = paste("SIMULATED COEFFICIENTS: Q-shape =", 
            x$qp), xlab = "m = Multicollinearity Allowance", 
            ylab = "Simulated Coefficients")
        if( trkey )
            legend("bottom",all.vars(x$form)[2:(x$p+1)], col=1:(x$p), lty=1:(x$p), lwd=2)
    }
    if (trace == "seq") {
        cat("\nPress the Enter key to view the RSE Loss trace...")
        scan()
    }
    if (trace == "all" || trace == "seq" || trace == "rsel") {
        plot(mcalp, x$rsel, ann = FALSE, type = "n")
        abline(h = 0, col = gray(0.9))
        for (i in 1:x$p) lines(mcal, x$rsel[, i], col = i, lty = i, 
            lwd = 2)
        title(main = paste("RELATIVE SqError LOSS: Q-shape =", 
            x$qp), xlab = "m = Multicollinearity Allowance", 
            ylab = "Scaled True SE Loss")
        if( trkey )
            legend("bottom",all.vars(x$form)[2:(x$p+1)], col=1:(x$p), lty=1:(x$p), lwd=2)
    }
    if (trace == "seq") {
        cat("\nPress the Enter key to view the SPAT trace...")
        scan()
    }
    if (trace == "all" || trace == "seq" || trace == "spat") {
        plot(mcalp, x$spat, ann = FALSE, type = "n")
        abline(h = 0, col = gray(0.9))
        for (i in 1:x$p) lines(mcal, x$spat[, i], col = i, lty = i, 
            lwd = 2)
        title(main = paste("SHRINKAGE PATTERN: Q-shape =", x$qp), 
            xlab = "m = Multicollinearity Allowance", ylab = "Ridge Delta Factors")
        if( trkey )
            legend("bottom",paste("Component",1:(x$p)), col=1:(x$p), lty=1:(x$p), lwd=2)
    }
}

"plot.RXuclars" <-
function (x, trace = "all", trkey = FALSE, ...) 
{
    mcal <- x$sext[, 2]
    mcalp <- rep(mcal, times = x$p)
    if (trace != "coef" && trace != "rmse" && trace != "exev" && 
        trace != "infd" && trace != "spat" && trace != "seq")
        trace <- "all"
    opar <- par(no.readonly = TRUE)
    on.exit(par(opar))
    if (trace == "all")
        par(mfrow=c(3,2))
    else
        par(mfrow=c(1,1))
    if (trace == "all" || trace == "seq" || trace == "coef") {
        plot(mcalp, x$coef, ann = FALSE, type = "n")
        abline(h = 0, col = gray(0.9))
        for (i in 1:x$p) lines(mcal, x$coef[, i], col = i, lty = i, 
            lwd = 2)
        title(main = paste("COEFFICIENT TRACE:", x$lars$type), 
            xlab = "m = Multicollinearity Allowance", ylab = "Fitted Coefficients")
        if( trkey )
            legend("bottom",all.vars(x$form)[2:(x$p+1)], col=1:(x$p), lty=1:(x$p), lwd=2)
    }
    if (trace == "seq") {
        cat("\nPress the Enter key to view the RMSE trace...")
        scan()
    }
    if (trace == "all" || trace == "seq" || trace == "rmse") {
        plot(mcalp, x$rmse, ann = FALSE, type = "n")
        abline(h = 0, col = gray(0.9))
        for (i in 1:x$p) lines(mcal, x$rmse[, i], col = i, lty = i, 
            lwd = 2)
        title(main = paste("RELATIVE MEAN SQ. ERROR:", x$lars$type), 
            xlab = "m = Multicollinearity Allowance", ylab = "Scaled MSE Risk")
        if( trkey )
            legend("bottom",all.vars(x$form)[2:(x$p+1)], col=1:(x$p), lty=1:(x$p), lwd=2)
    }
    if (trace == "seq") {
        cat("\nPress the Enter key to view the EXEV trace...")
        scan()
    }
    if (trace == "all" || trace == "seq" || trace == "exev") {
        plot(mcalp, x$exev, ann = FALSE, type = "n")
        abline(h = 0, col = gray(0.9))
        for (i in 1:x$p) lines(mcal, x$exev[, i], col = i, lty = i, 
            lwd = 2)
        title(main = paste("EXCESS EIGENVALUES:", x$lars$type), 
            xlab = "m = Multicollinearity Allowance", ylab = "Least Squares minus uclars")
        if( trkey )
            legend("bottom",paste("Component",1:(x$p)), col=1:(x$p), lty=1:(x$p), lwd=2)
    }
    if (trace == "seq") {
        cat("\nPress the Enter key to view the INFD trace...")
        scan()
    }
    if (trace == "all" || trace == "seq" || trace == "infd") {
        plot(mcalp, x$infd, ann = FALSE, type = "n")
        abline(h = 0, col = gray(0.9))
        for (i in 1:x$p) lines(mcal, x$infd[, i], col = i, lty = i, 
            lwd = 2)
        title(main = paste("INFERIOR DIRECTION:", x$lars$type), 
            xlab = "m = Multicollinearity Allowance", ylab = "Direction Cosines")
        if( trkey )
            legend("bottom",all.vars(x$form)[2:(x$p+1)], col=1:(x$p), lty=1:(x$p), lwd=2)
    }
    if (trace == "seq") {
        cat("\nPress the Enter key to view the SPAT trace...")
        scan()
    }
    if (trace == "all" || trace == "seq" || trace == "spat") {
        plot(mcalp, x$spat, ann = FALSE, type = "n")
        abline(h = 0, col = gray(0.9))
        for (i in 1:x$p) lines(mcal, x$spat[, i], col = i, lty = i, 
            lwd = 2)
        title(main = paste("SHRINKAGE PATTERN:", x$lars$type), 
            xlab = "m = Multicollinearity Allowance", ylab = "uclars Delta Factors")
        if( trkey )
            legend("bottom",paste("Component",1:(x$p)), col=1:(x$p), lty=1:(x$p), lwd=2)
    }
}

"print.RXlarlso" <-
function (x, ...) 
{
    cat("\nRXlarlso Object: LARS Maximum Likelihood Shrinkage\n")
    cat("Data Frame:", x$data, "\n")
    cat("Regression Equation:\n")
    print(x$form)
    cat("\n    Number of Regressor Variables, p =", x$p, "\n")
    cat("    Number of Observations, n =", x$n, "\n")
    cat("\nPrincipal Axis Summary Statistics of Ill-Conditioning...\n")
    print.default(x$prinstat, quote = FALSE)
    cat("\n    Residual Mean Square for Error =", x$s2, "\n")
    cat("    Estimate of Residual Std. Error =", sqrt(x$s2), 
        "\n\n")
    cat("\nThe extent of shrinkage (M value) most likely to be optimal\n")
    cat("depends upon whether one uses the Classical, Empirical Bayes, or\n")
    cat("Random Coefficient criterion.  In each case, the objective is to\n")
    cat("minimize the minus-two-log-likelihood statistics listed below:\n")
    print.default(x$mlik, quote = FALSE)
    cat("\nExtent of shrinkage statistics...\n")
    print.default(x$sext, quote = FALSE)
    cat("\nOutput from LARS invocation...\n")
    print(x$lars)
}

"print.RXridge" <-
function (x, ...) 
{
    cat("\nRXridge Object: Shrinkage-Ridge Regression Model Specification\n")
    cat("Data Frame:", x$data, "\n")
    cat("Regression Equation:\n")
    print(x$form)
    cat("\n    Number of Regressor Variables, p =", x$p, "\n")
    cat("    Number of Observations, n =", x$n, "\n")
    cat("\nPrincipal Axis Summary Statistics of Ill-Conditioning...\n")
    print.default(x$prinstat, quote = FALSE)
    cat("\n    Residual Mean Square for Error =", x$s2, "\n")
    cat("    Estimate of Residual Std. Error =", sqrt(x$s2), 
        "\n\n")
    cat("Classical Maximum Likelihood choice of SHAPE(Q) and EXTENT(M) of\n")
    cat("shrinkage in the 2-parameter generalized ridge family...\n")
    print.default(x$crlqstat, quote = FALSE)
    cat("\n Q =", x$qmse, " is the path shape most likely to lead to minimum\n")
    cat("MSE risk because this shape maximizes CRLQ and minimizes CHISQ.\n")
    cat("\n\nRXridge: Shrinkage PATH Shape =", x$qp, "\n")
    cat("\nThe extent of shrinkage (M value) most likely to be optimal\n")
    cat("in the Q-shape =", x$qp, " 2-parameter ridge family can depend\n")
    cat("upon whether one uses the Classical, Empirical Bayes, or Random\n")
    cat("Coefficient criterion.  In each case, the objective is to\n")
    cat("minimize the minus-two-log-likelihood statistics listed below:\n")
    print.default(x$mlik, quote = FALSE)
    cat("\nExtent of shrinkage statistics...\n")
    print.default(x$sext, quote = FALSE)
}

"print.RXtrisk" <-
function (x, ...) 
{
    cat("\nRXtrisk Object: True Risk of Shrinkage in Regression\n")
    cat("Data Frame:", x$data, "\n")
    cat("Regression Equation:\n")
    print(x$form)
    cat("\n    Number of Regressor Variables, p =", x$p, "\n")
    cat("    Number of Observations, n =", x$n, "\n")
    cat("\nPrincipal Axis Summary Statistics of Ill-Conditioning...\n")
    print.default(x$prinstat, quote = FALSE)
    cat("\n\nRXtrisk: Shrinkage PATH Shape =", x$qp, "\n")
    cat("\nExtent of shrinkage statistics...\n")
    print.default(x$sext, quote = FALSE)
}

"print.RXtsimu" <-
function (x, ...) 
{
    cat("\nRXtsimu Object: Simulated Loss of Shrinkage in Regression\n")
    cat("Data Frame:", x$data, "\n")
    cat("Regression Equation:\n")
    print(x$form)
    cat("\n    Number of Regressor Variables, p =", x$p, "\n")
    cat("    Number of Observations, n =", x$n, "\n")
    cat("\nPrincipal Axis Summary Statistics of Ill-Conditioning...\n")
    print.default(x$prinstat, quote = FALSE)
    cat("\n\nRXtsimu: Shrinkage PATH Shape =", x$qp, "\n")
    cat("\nExtent of shrinkage statistics...\n")
    print.default(x$sext, quote = FALSE)
    cat("\nSimulated Response and Expected Values...\n")
    print.default(x$ydat, quote = FALSE)
}

"print.RXuclars" <-
function (x, ...) 
{
    cat("\nRXuclars Object: Uncorrelated Component LARS Shrinkage\n")
    cat("Data Frame:", x$data, "\n")
    cat("Regression Equation:\n")
    print(x$form)
    cat("\n    Number of Regressor Variables, p =", x$p, "\n")
    cat("    Number of Observations, n =", x$n, "\n")
    cat("\nPrincipal Axis Summary Statistics of Ill-Conditioning...\n")
    print.default(x$prinstat, quote = FALSE)
    cat("\n    Residual Mean Square for Error =", x$s2, "\n")
    cat("    Estimate of Residual Std. Error =", sqrt(x$s2), 
        "\n\n")
    cat("\nThe extent of shrinkage (M value) most likely to be optimal\n")
    cat("depends upon whether one uses the Classical, Empirical Bayes, or\n")
    cat("Random Coefficient criterion.  In each case, the objective is to\n")
    cat("minimize the minus-two-log-likelihood statistics listed below:\n")
    print.default(x$mlik, quote = FALSE)
    cat("\nExtent of shrinkage statistics...\n")
    print.default(x$sext, quote = FALSE)
    cat("\nOutput from LARS invocation...\n")
    print(x$lars)
}

"RXlarlso" <-
function (form, data, rscale = 1, type = "lar", trace = FALSE,
          eps = .Machine$double.eps, omdmin = 9.9e-13, ...) 
{
    if (missing(form) || class(form) != "formula") 
        stop("First argument to RXlarlso must be a valid linear regression formula.")
    yvar <- deparse(form[[2]])
    if (missing(data) || !inherits(data, "data.frame")) 
        stop("Second argument to RXlarlso must be an existing Data Frame.")
    dfname <- deparse(substitute(data))
    if (!is.element(yvar, dimnames(data)[[2]])) 
        stop("Response variable in the RXlarlso formula must be an existing variable.")
    lmobj <- lm(form, data)
    yvec <- as.matrix(lmobj$model[, 1])
    xmat <- as.matrix(lmobj$model[, 2:length(lmobj$model)])
    if (rscale < 1) 
        rscale <- 0
    if (rscale > 1) 
        rscale <- 2
    p <- ncol(xmat)
    if (p < 2) 
        stop("Number on non-constant regressor variables must be at least 2.")
    n <- nrow(xmat)
    if (n != nrow(yvec)) 
        stop("Numbers of observations in XMAT and YVEC must match.")
    if (n < p + 4) 
        stop("Number of observations must exceed number of regressors by at least 4.")
    mx <- matrix(apply(xmat, 2, "mean"), nrow = 1)
    crx <- xmat - matrix(1, n, 1) %*% mx
    xscale <- diag(p)
    if (rscale >= 1) {
        xscale <- diag(sqrt(diag(var(crx))))
        crx <- crx %*% solve(xscale)
    }
    sx <- svd(crx)
    eigval <- matrix(sx$d^2, ncol = 1)
    eiginv <- solve(diag(sx$d^2, ncol = p))
    cry <- matrix(yvec - mean(yvec), ncol = 1)
    yscale <- 1
    if (rscale >= 1) {
        yscale <- sqrt(var(cry))
        cry <- cry/yscale[1, 1]
    }
    smse <- sx$v %*% eiginv %*% t(sx$v)
    risk <- diag(smse)
    tsmse <- sum(risk)
    comp <- solve(diag(sx$d, ncol = p)) %*% t(sx$u) %*% cry
    bstar <- sx$v %*% comp
    exev <- matrix(0, p, 1)
    infd <- matrix(0, p, 1)
    delta <- matrix(1, p, 1)
    idty <- diag(p)
    d <- idty
    cold <- delta
    sv <- matrix(sx$d, p, 1)
    ssy <- t(cry) %*% cry
    rho <- (sv * comp)/sqrt(ssy[1, 1])
    arho <- matrix(abs(rho), nrow = 1)
    r2 <- sum(arho^2)
    if (r2 >= 1) 
        stop(" Maximum likelihood shrinkage cannot be performed when RSQUARE=1.")
    res <- cry - crx %*% bstar
    s2 <- t(res) %*% res/(n - p - 1)
    varrho <- s2[1, 1]/ssy[1, 1]
    tstat <- rho/sqrt(varrho)
    frat <- rho^2/varrho
    stat <- cbind(eigval, sv, comp, rho, tstat)
    dimnames(stat) <- list(1:p, c("LAMBDA", "SV", "COMP", "RHO", 
        "TRAT"))
    RXolist <- list(data = dfname, form = form, p = p, n = n, 
        r2 = r2, s2 = s2, prinstat = stat, gmat = sx$v)
    larsobj <- lars(crx, cry, type, trace, eps, ...)
    bhat <- as.matrix(larsobj$beta)
    if (p != length(bhat[1, ])) 
        stop("Number of coefficients for LARS and shrinkage must match.")
    steps <- length(bhat[, 1])
    binc <- bstar - bhat[steps, ]
    if (sum(binc^2) > omdmin) 
        stop("OLS coefficients for LARS and shrinkage must match.")
    mcal <- 0
    const <- (n - p - 3)/(n - p - 1)
    srat <- solve(diag(as.vector(sv), ncol = p)) %*% tstat
    MCAL <- 0
    C <- Inf
    E <- Inf
    R <- Inf
    for (inc in 2:steps) {
        binc <- bhat[steps + 1 - inc, ]
        dgen <- (t(sx$v) %*% binc)/comp
        dinc <- pmax(pmin(dgen, 1), 0)
        d <- matrix(dinc, p, p)
        minc <- p - sum(dinc)
        omd <- pmax(1 - dinc, omdmin)
        ddomd <- dinc/omd
        rxi <- sum(t(arho) * sqrt(ddomd))
        slik <- 2/(rxi + sqrt(4 * n + rxi^2))
        clik <- 2 * n * log(slik) + sum(ddomd) - (rxi/slik) - 
            n * log((1 - r2)/n)
        ebay <- sum(frat * omd - log(omd))
        sr2d <- sum(dinc * rho^2)
        rcof <- -sum(log(omd)) + n * log((1 - sr2d)/(1 - r2))
        C <- rbind(C, clik)
        E <- rbind(E, ebay)
        R <- rbind(R, rcof)
        MCAL <- rbind(MCAL, minc)
        vecr <- (idty - d) %*% srat
        compr <- const * vecr %*% t(vecr) + (2 * d - idty) * 
            eiginv
        diagc <- diag(diag(compr), ncol = p)
        lowr <- eiginv * d^2
        maxd <- matrix(pmax(diagc, lowr), p, p)
        compr <- compr - diagc + maxd
        smse <- sx$v %*% compr %*% t(sx$v)
        rinc <- diag(smse)
        lowb <- sx$v %*% lowr %*% t(sx$v)
        rinc <- pmax(rinc, diag(lowb))
        tinc <- sum(rinc)
        emse <- eiginv - compr
        sfac <- min(abs(diag(emse)))/100
        if (sfac < 1e-05) 
            sfac <- 1e-05
        eign <- eigen(emse/sfac)
        einc <- rev(eign$values) * sfac
        cinc <- matrix(0, p, 1)
        if (is.na(einc[1])) 
            einc[1] <- 0
        if (einc[1] < 0) {
            eign$vectors <- sx$v %*% eign$vectors
            cinc <- eign$vectors[, p]
            if (rscale == 2) {
                cinc <- cinc %*% xscale
                cinc <- cinc/sqrt(sum(cinc^2))
            }
            if (t(cold) %*% cinc < 0) 
                cinc <- -1 * cinc
            cold <- cinc
        }
        bstar <- cbind(bstar, binc)
        risk <- cbind(risk, rinc)
        exev <- cbind(exev, einc)
        infd <- cbind(infd, cinc)
        delta <- cbind(delta, dgen)
        tsmse <- rbind(tsmse, tinc)
        mcal <- rbind(mcal, minc)
    }
    if (rscale == 2) {
        bstar <- yscale * solve(xscale) %*% bstar
        risk <- yscale^2 * solve(xscale^2) %*% risk
    }
    mlik <- cbind(MCAL, C, E, R)
    dimnames(mlik) <- list(0:(steps - 1), c("M", "CLIK", "EBAY", 
        "RCOF"))
    bstar <- t(bstar)
    risk <- t(risk)
    exev <- t(exev)
    infd <- t(infd)
    delta <- t(delta)
    sext <- cbind(tsmse, mcal)
    dimnames(sext) <- list(0:(steps - 1), c("TSMSE", "MCAL"))
    RXolist <- c(RXolist, list(lars = larsobj, coef = bstar, 
        rmse = risk, exev = exev, infd = infd, spat = delta, 
        mlik = mlik, sext = sext))
    class(RXolist) <- "RXlarlso"
    RXolist
}

"RXmstep" <-
function (mobj, kinc, p, qp, eqm1) 
{
    if (mobj <= 0) {
        d <- diag(p)
        kinc <- 0
        return(list(kinc = kinc, d = d))
    }
    if (mobj >= p) {
        d <- matrix(0, p, p)
        kinc <- Inf
        return(list(kinc = kinc, d = d))
    }
    funs <- mobj - p
    if (qp == 1) 
        kinc <- (-1 * mobj)/funs
    else {
        while (abs(funs) > 1e-05) {
            funs <- mobj - p + sum((1 + kinc * eqm1)^-1)
            derivs <- sum(eqm1/(1 + kinc * eqm1)^2)
            kinc <- kinc + funs/derivs
        }
    }
    d <- diag(as.vector(1/(1 + kinc * eqm1)), p)
    iter <- list(kinc = kinc, d = d)
    iter
}

"RXridge" <-
function (form, data, rscale = 1, Q = "qmse", steps = 8, nq = 21, 
    qmax = 5, qmin = -5, omdmin = 9.9e-13) 
{
    if (missing(form) || class(form) != "formula") 
        stop("First argument to RXridge must be a valid linear regression formula.")
    yvar <- deparse(form[[2]])
    if (missing(data) || !inherits(data, "data.frame")) 
        stop("Second argument to RXridge must be an existing Data Frame.")
    dfname <- deparse(substitute(data))
    if (!is.element(yvar, dimnames(data)[[2]])) 
        stop("Response variable in the RXridge formula must be an existing variable.")
    lmobj <- lm(form, data)
    yvec <- as.matrix(lmobj$model[, 1])
    xmat <- as.matrix(lmobj$model[, 2:length(lmobj$model)])
    nqval <- round(nq)
    if (nqval < 3) 
        nqval <- 3
    if (rscale < 1) 
        rscale <- 0
    if (rscale > 1) 
        rscale <- 2
    p <- ncol(xmat)
    if (p < 2) 
        stop("Number on non-constant regressor variables must be at least 2.")
    n <- nrow(xmat)
    if (n != nrow(yvec)) 
        stop("Numbers of observations in XMAT and YVEC must match.")
    if (n < p + 4) 
        stop("Number of observations must exceed number of regressors by at least 4.")
    mx <- matrix(apply(xmat, 2, "mean"), nrow = 1)
    crx <- xmat - matrix(1, n, 1) %*% mx
    xscale <- diag(p)
    if (rscale >= 1) {
        xscale <- diag(sqrt(diag(var(crx))))
        crx <- crx %*% solve(xscale)
    }
    sx <- svd(crx)
    eigval <- matrix(sx$d^2, ncol = 1)
    eiginv <- solve(diag(sx$d^2, ncol = p))
    cry <- matrix(yvec - mean(yvec), ncol = 1)
    yscale <- 1
    if (rscale >= 1) {
        yscale <- sqrt(var(cry))
        cry <- cry/yscale[1, 1]
    }
    smse <- sx$v %*% eiginv %*% t(sx$v)
    risk <- diag(smse)
    tsmse <- sum(risk)
    comp <- solve(diag(sx$d, ncol = p)) %*% t(sx$u) %*% cry
    bstar <- sx$v %*% comp
    exev <- matrix(0, p, 1)
    infd <- matrix(0, p, 1)
    delta <- matrix(1, p, 1)
    idty <- diag(p)
    d <- idty
    cold <- delta
    sv <- matrix(sx$d, p, 1)
    ssy <- t(cry) %*% cry
    rho <- (sv * comp)/sqrt(ssy[1, 1])
    arho <- matrix(abs(rho), nrow = 1)
    r2 <- sum(arho^2)
    if (r2 >= 1) 
        stop(" Maximum likelihood ridge shrinkage cannot be applied when RSQUARE=1.")
    res <- cry - crx %*% bstar
    s2 <- t(res) %*% res/(n - p - 1)
    varrho <- s2[1, 1]/ssy[1, 1]
    tstat <- rho/sqrt(varrho)
    frat <- rho^2/varrho
    stat <- cbind(eigval, sv, comp, rho, tstat)
    dimnames(stat) <- list(1:p, c("LAMBDA", "SV", "COMP", "RHO", 
        "TRAT"))
    RXolist <- list(data = dfname, form = form, p = p, n = n, 
        r2 = r2, s2 = s2, prinstat = stat)
    qvec <- matrix(0, nqval, 1)
    crlq <- matrix(0, nqval, 1)
    mvec <- matrix(0, nqval, 1)
    kvec <- matrix(0, nqval, 1)
    chisq <- matrix(0, nqval, 1)
    for (it in 1:nqval) {
        qnow <- ((qmin - qmax) * it + (nqval * qmax) - qmin)/(nqval - 
            1)
        qvec[it] <- qnow
        s1mq <- matrix(1, p, 1)
        if (qnow != 1) 
            s1mq <- exp((1 - qnow) * log(sv))
        sq2 <- sum(s1mq^2)
        crlq[it] <- (arho %*% s1mq)/sqrt(r2 * sq2)
        r2c2 <- r2 * crlq[it]^2
        kvec[it] <- (sq2 * (1 - r2c2))/(n * r2c2)
        s1mq <- kvec[it] * s1mq^-2
        mvec[it] <- sum(s1mq/(1 + s1mq))
        chisq[it] <- n * log((1 - r2c2)/(1 - r2))
    }
    c2min <- chisq[1]
    qmse <- qvec[1]
    for (it in 2:nqval) {
        if (chisq[it] < c2min) {
            c2min <- chisq[it]
            qmse <- qvec[it]
        }
    }
    stat <- cbind(qvec, crlq, mvec, kvec, chisq)
    dimnames(stat) <- list(1:nqval, c("Q", "CRLQ", "M", "K", 
        "CHISQ"))
    RXolist <- c(RXolist, list(crlqstat = stat, qmse = qmse))
    if (missing(Q) || Q == "qmse") 
        qp <- qmse
    else {
        Q <- as.numeric(Q)
        if (is.na(Q)) 
            Q <- 0
        qp <- min(max(Q, qmin), qmax)
    }
    qm1 <- qp - 1
    if (qp == 1) 
        eqm1 <- matrix(1, p, 1)
    else eqm1 <- exp(qm1 * log(eigval))
    mcal <- 0
    konst <- 0
    kinc <- 0
    const <- (n - p - 3)/(n - p - 1)
    srat <- solve(diag(as.vector(sv), ncol = p)) %*% tstat
    MCAL <- 0
    KONST <- 0
    C <- Inf
    E <- Inf
    R <- Inf
    maxinc <- p * steps
    for (inc in 1:maxinc) {
        mobj <- inc/steps
        iter <- RXmstep(mobj, kinc, p, qp, eqm1)
        kinc <- iter$kinc
        d <- matrix(iter$d, p, p)
        dinc <- diag(d)
        omd <- pmax(1 - dinc, omdmin)
        ddomd <- dinc/omd
        rxi <- sum(arho * sqrt(ddomd))
        slik <- 2/(rxi + sqrt(4 * n + rxi^2))
        clik <- 2 * n * log(slik) + sum(ddomd) - (rxi/slik) - 
            n * log((1 - r2)/n)
        ebay <- sum(frat * omd - log(omd))
        sr2d <- sum(dinc * rho^2)
        rcof <- -sum(log(omd)) + n * log((1 - sr2d)/(1 - r2))
        minc <- p - sum(dinc)
        C <- rbind(C, clik)
        E <- rbind(E, ebay)
        R <- rbind(R, rcof)
        KONST <- rbind(KONST, kinc)
        MCAL <- rbind(MCAL, minc)
        binc <- sx$v %*% d %*% comp
        vecr <- (idty - d) %*% srat
        compr <- const * vecr %*% t(vecr) + (2 * d - idty) * 
            eiginv
        diagc <- diag(diag(compr), ncol = p)
        lowr <- eiginv * d^2
        maxd <- matrix(pmax(diagc, lowr), p, p)
        compr <- compr - diagc + maxd
        smse <- sx$v %*% compr %*% t(sx$v)
        rinc <- diag(smse)
        lowb <- sx$v %*% lowr %*% t(sx$v)
        rinc <- pmax(rinc, diag(lowb))
        tinc <- sum(rinc)
        emse <- eiginv - compr
        sfac <- min(abs(diag(emse)))/100
        if (sfac < 1e-05) 
            sfac <- 1e-05
        eign <- eigen(emse/sfac)
        einc <- rev(eign$values) * sfac
        cinc <- matrix(0, p, 1)
        if (is.na(einc[1])) 
            einc[1] <- 0
        if (einc[1] < 0) {
            eign$vectors <- sx$v %*% eign$vectors
            cinc <- eign$vectors[, p]
            if (rscale == 2) {
                cinc <- cinc %*% xscale
                cinc <- cinc/sqrt(sum(cinc^2))
            }
            if (t(cold) %*% cinc < 0) 
                cinc <- -1 * cinc
            cold <- cinc
        }
        bstar <- cbind(bstar, binc)
        risk <- cbind(risk, rinc)
        exev <- cbind(exev, einc)
        infd <- cbind(infd, cinc)
        delta <- cbind(delta, dinc)
        tsmse <- rbind(tsmse, tinc)
        konst <- rbind(konst, kinc)
        mcal <- rbind(mcal, minc)
    }
    if (rscale == 2) {
        bstar <- yscale * solve(xscale) %*% bstar
        risk <- yscale^2 * solve(xscale^2) %*% risk
    }
    mlik <- cbind(MCAL, KONST, C, E, R)
    dimnames(mlik) <- list(0:maxinc, c("M", "K", "CLIK", "EBAY", 
        "RCOF"))
    bstar <- t(bstar)
    risk <- t(risk)
    exev <- t(exev)
    infd <- t(infd)
    delta <- t(delta)
    sext <- cbind(tsmse, konst, mcal)
    dimnames(sext) <- list(0:maxinc, c("TSMSE", "KONST", "MCAL"))
    RXolist <- c(RXolist, list(qp = qp, coef = bstar, rmse = risk, 
        exev = exev, infd = infd, spat = delta, mlik = mlik, 
        sext = sext))
    class(RXolist) <- "RXridge"
    RXolist
}

"RXtrisk" <-
function (form, data, trug, trus, Q = 0, rscale = 1, steps = 8, 
    qmax = 5, qmin = -5) 
{
    if (missing(form) || class(form) != "formula") 
        stop("First argument to RXtrisk must be a valid linear regression formula.")
    yvar <- deparse(form[[2]])
    if (missing(data) || !inherits(data, "data.frame")) 
        stop("Second argument to RXtrisk must be an existing Data Frame.")
    dfname <- deparse(substitute(data))
    lmobj <- lm(form, data)
    xmat <- as.matrix(lmobj$model[, 2:length(lmobj$model)])
    if (rscale < 1) 
        rscale <- 0
    if (rscale > 1) 
        rscale <- 2
    qp <- min(max(Q, qmin), qmax)
    n <- nrow(xmat)
    p <- ncol(xmat)
    if (p < 2) 
        stop("Number on non-constant regressor variables must be at least 2.")
    if (p != nrow(trug)) 
        stop("Number of true components must equal number of columns in XMAT.")
    mx <- matrix(apply(xmat, 2, "mean"), nrow = 1)
    crx <- xmat - matrix(1, n, 1) %*% mx
    xscale <- diag(p)
    if (rscale >= 1) {
        xscale <- diag(sqrt(diag(var(crx))))
        crx <- crx %*% solve(xscale)
    }
    sx <- svd(crx)
    eigval <- matrix(sx$d^2, ncol = 1)
    if (eigval[p] <= 0) {
        cat("\n\nRXtrisk requires XMAT to have Full Column Rank.")
        print(sx$d)
        print(sx$v)
        stop()
    }
    eiginv <- solve(diag(sx$d^2, ncol = p))
    smse <- sx$v %*% eiginv %*% t(sx$v)
    risk <- diag(smse)
    tsmse <- sum(risk)
    fact <- trus^2
    for (i in 1:p) fact <- fact + trug[i]^2 * eigval[i]/(n - 
        1)
    trug <- trug/sqrt(fact)
    trus <- trus/sqrt(fact)
    beta <- sx$v %*% trug
    maxinc <- p * steps
    qm1 <- qp - 1
    exev <- matrix(0, p, 1)
    infd <- matrix(0, p, 1)
    delta <- matrix(1, p, 1)
    d <- diag(delta)
    cold <- delta
    sv <- matrix(sx$d, p, 1)
    idty <- diag(p)
    if (qp == 1) 
        eqm1 <- matrix(1, p, 1)
    else eqm1 <- exp(qm1 * log(eigval))
    stat <- cbind(eigval, sv, trug, beta)
    dimnames(stat) <- list(1:p, c("LAMBDA", "SV", "GAMMA", "BETA"))
    RXolist <- list(data = dfname, form = form, trug = trug, 
        trus = trus, qp = qp, p = p, n = n, prinstat = stat)
    mcal <- 0
    konst <- 0
    kinc <- 0
    for (inc in 1:maxinc) {
        mobj <- inc/steps
        iter <- RXmstep(mobj, kinc, p, qp, eqm1)
        kinc <- iter$kinc
        d <- matrix(iter$d, p, p)
        binc <- sx$v %*% d %*% trug
        rbias <- (idty - d) %*% trug/trus
        compr <- rbias %*% t(rbias) + eiginv * d^2
        smse <- sx$v %*% compr %*% t(sx$v)
        rinc <- diag(smse)
        tinc <- sum(rinc)
        emse <- eiginv - compr
        sfac <- min(abs(diag(emse)))/100
        if (sfac < 1e-05) 
            sfac <- 1e-05
        eign <- eigen(emse/sfac)
        einc <- rev(eign$values) * sfac
        cinc <- matrix(0, p, 1)
        if (is.na(einc[1])) 
            einc[1] <- 0
        if (einc[1] < 0) {
            eign$vectors <- sx$v %*% eign$vectors
            cinc <- eign$vectors[, p]
            if (rscale == 2) {
                cinc <- cinc %*% xscale
                cinc <- cinc/sqrt(sum(cinc^2))
            }
            if (t(cold) %*% cinc < 0) 
                cinc <- -1 * cinc
            cold <- cinc
        }
        dinc <- diag(d)
        minc <- p - sum(dinc)
        beta <- cbind(beta, binc)
        risk <- cbind(risk, rinc)
        exev <- cbind(exev, einc)
        infd <- cbind(infd, cinc)
        delta <- cbind(delta, dinc)
        tsmse <- rbind(tsmse, tinc)
        konst <- rbind(konst, kinc)
        mcal <- rbind(mcal, minc)
    }
    beta <- t(beta)
    risk <- t(risk)
    exev <- t(exev)
    infd <- t(infd)
    delta <- t(delta)
    sext <- cbind(tsmse, konst, mcal)
    dimnames(sext) <- list(0:maxinc, c("TSMSE", "KONST", "MCAL"))
    RXolist <- c(RXolist, list(coef = beta, rmse = risk, exev = exev, 
        infd = infd, spat = delta, sext = sext))
    class(RXolist) <- "RXtrisk"
    RXolist
}

"RXtsimu" <-
function (form, data, trug, trus, Q = 0, rscale = 1, steps = 8, 
    seed = 123, qmax = 5, qmin = -5) 
{
    if (missing(form) || class(form) != "formula") 
        stop("First argument to RXtsimu must be a valid linear regression formula.")
    yvar <- deparse(form[[2]])
    if (missing(data) || !inherits(data, "data.frame")) 
        stop("Second argument to RXtsimu must be an existing Data Frame.")
    dfname <- deparse(substitute(data))
    lmobj <- lm(form, data)
    xmat <- as.matrix(lmobj$model[, 2:length(lmobj$model)])
    if (rscale < 1) 
        rscale <- 0
    if (rscale > 1) 
        rscale <- 2
    qp <- min(max(Q, qmin), qmax)
    n <- nrow(xmat)
    p <- ncol(xmat)
    if (p < 2) 
        stop("Number on non-constant regressor variables must be at least 2.")
    if (p != nrow(trug)) 
        stop("Number of true components must equal number of columns in XMAT.")
    mx <- matrix(apply(xmat, 2, "mean"), nrow = 1)
    crx <- xmat - matrix(1, n, 1) %*% mx
    xscale <- diag(p)
    if (rscale >= 1) {
        xscale <- diag(sqrt(diag(var(crx))))
        crx <- crx %*% solve(xscale)
    }
    sx <- svd(crx)
    eigval <- matrix(sx$d^2, ncol = 1)
    if (eigval[p] <= 0) {
        cat("\n\nRXtsimu requires XMAT to have Full Column Rank.")
        print(sx$d)
        print(sx$v)
        stop()
    }
    eiginv <- solve(diag(sx$d^2, ncol = p))
    fact <- trus^2
    for (i in 1:p) fact <- fact + trug[i]^2 * eigval[i]/(n - 
        1)
    trug <- trug/sqrt(fact)
    trus <- trus/sqrt(fact)
    beta <- sx$v %*% trug
    expy <- crx %*% beta
    if (missing(seed)) 
        seed <- sample(1001, 1) - 1
    set.seed(seed)
    error <- matrix(rnorm(n, 0, trus), ncol = 1)
    yvec <- matrix(expy + error - mean(error), ncol = 1)
    ydat <- cbind(yvec, expy)
    dimnames(ydat) <- list(1:n, c("YVEC", "EXPY"))
    ssy <- t(yvec) %*% yvec
    comp <- solve(diag(sx$d)) %*% t(sx$u) %*% yvec
    bstar <- sx$v %*% comp
    loss <- (bstar - beta)^2/trus^2
    totlos <- sum(loss)
    maxinc <- p * steps
    qm1 <- qp - 1
    delta <- matrix(1, p, 1)
    d <- diag(p)
    sv <- matrix(sx$d, p, 1)
    idty <- diag(p)
    if (qp == 1) 
        eqm1 <- matrix(1, p, 1)
    else eqm1 <- exp(qm1 * log(eigval))
    stat <- cbind(eigval, sv, trug, beta)
    dimnames(stat) <- list(1:p, c("LAMBDA", "SV", "GAMMA", "BETA"))
    RXolist <- list(data = dfname, form = form, trug = trug, 
        trus = trus, qp = qp, p = p, n = n, ydat = ydat, prinstat = stat)
    rho <- (sv * comp)/sqrt(ssy[1, 1])
    res <- yvec - crx %*% bstar
    s2 <- t(res) %*% res/(n - p - 1)
    varrho <- s2[1, 1]/ssy[1, 1]
    trat <- rho/sqrt(varrho)
    mcal <- 0
    konst <- 0
    kinc <- 0
    for (inc in 1:maxinc) {
        mobj <- inc/steps
        iter <- RXmstep(mobj, kinc, p, qp, eqm1)
        kinc <- iter$kinc
        d <- matrix(iter$d, p, p)
        binc <- sx$v %*% d %*% comp
        linc <- (binc - beta)^2/trus^2
        tinc <- sum(linc)
        dinc <- diag(d)
        minc <- p - sum(dinc)
        bstar <- cbind(bstar, binc)
        loss <- cbind(loss, linc)
        delta <- cbind(delta, dinc)
        totlos <- rbind(totlos, tinc)
        konst <- rbind(konst, kinc)
        mcal <- rbind(mcal, minc)
    }
    bstar <- t(bstar)
    loss <- t(loss)
    delta <- t(delta)
    sext <- cbind(totlos, konst, mcal)
    dimnames(sext) <- list(0:maxinc, c("TLOSS", "KONST", "MCAL"))
    RXolist <- c(RXolist, list(coef = bstar, rsel = loss, spat = delta, 
        sext = sext))
    class(RXolist) <- "RXtsimu"
    RXolist
}

"RXuclars" <-
function (form, data, rscale = 1, type = "lar", trace = FALSE, 
    eps = .Machine$double.eps, omdmin = 9.9e-13, ...) 
{
    if (missing(form) || class(form) != "formula") 
        stop("First argument to RXuclars must be a valid linear regression formula.")
    yvar <- deparse(form[[2]])
    if (missing(data) || !inherits(data, "data.frame")) 
        stop("Second argument to RXuclars must be an existing Data Frame.")
    dfname <- deparse(substitute(data))
    if (!is.element(yvar, dimnames(data)[[2]])) 
        stop("Response variable in the RXuclars formula must be an existing variable.")
    lmobj <- lm(form, data)
    yvec <- as.matrix(lmobj$model[, 1])
    xmat <- as.matrix(lmobj$model[, 2:length(lmobj$model)])
    if (rscale < 1) 
        rscale <- 0
    if (rscale > 1) 
        rscale <- 2
    p <- ncol(xmat)
    if (p < 2) 
        stop("Number on non-constant regressor variables must be at least 2.")
    n <- nrow(xmat)
    if (n != nrow(yvec)) 
        stop("Numbers of observations in XMAT and YVEC must match.")
    if (n < p + 4) 
        stop("Number of observations must exceed number of regressors by at least 4.")
    mx <- matrix(apply(xmat, 2, "mean"), nrow = 1)
    crx <- xmat - matrix(1, n, 1) %*% mx
    xscale <- diag(p)
    if (rscale >= 1) {
        xscale <- diag(sqrt(diag(var(crx))))
        crx <- crx %*% solve(xscale)
    }
    sx <- svd(crx)
    eigval <- matrix(sx$d^2, ncol = 1)
    eiginv <- solve(diag(sx$d^2, ncol = p))
    cry <- matrix(yvec - mean(yvec), ncol = 1)
    yscale <- 1
    if (rscale >= 1) {
        yscale <- sqrt(var(cry))
        cry <- cry/yscale[1, 1]
    }
    smse <- sx$v %*% eiginv %*% t(sx$v)
    risk <- diag(smse)
    tsmse <- sum(risk)
    comp <- solve(diag(sx$d, ncol = p)) %*% t(sx$u) %*% cry
    bstar <- sx$v %*% comp
    exev <- matrix(0, p, 1)
    infd <- matrix(0, p, 1)
    delta <- matrix(1, p, 1)
    idty <- diag(p)
    d <- idty
    cold <- delta
    sv <- matrix(sx$d, p, 1)
    ssy <- t(cry) %*% cry
    rho <- (sv * comp)/sqrt(ssy[1, 1])
    arho <- matrix(abs(rho), nrow = 1)
    r2 <- sum(arho^2)
    if (r2 >= 1) 
        stop(" Maximum likelihood shrinkage cannot be performed when RSQUARE=1.")
    res <- cry - crx %*% bstar
    s2 <- t(res) %*% res/(n - p - 1)
    varrho <- s2[1, 1]/ssy[1, 1]
    tstat <- rho/sqrt(varrho)
    frat <- rho^2/varrho
    stat <- cbind(eigval, sv, comp, rho, tstat)
    dimnames(stat) <- list(1:p, c("LAMBDA", "SV", "COMP", "RHO", 
        "TRAT"))
    RXolist <- list(data = dfname, form = form, p = p, n = n, 
        r2 = r2, s2 = s2, prinstat = stat, gmat = sx$v)
    larsobj <- lars(sx$u, cry, type, trace, eps, ...) 
    bhat <- as.matrix(larsobj$beta) %*% solve(diag(sx$d, ncol = p)) %*% 
        t(sx$v)
    if (p != length(bhat[1, ])) 
        stop("Number of coefficients for LARS and shrinkage must match.")
    steps <- length(bhat[, 1])
    binc <- bstar - bhat[steps, ]
    if (sum(binc^2) > omdmin) 
        stop("OLS coefficients for LARS and shrinkage must match.")
    mcal <- 0
    const <- (n - p - 3)/(n - p - 1)
    srat <- solve(diag(as.vector(sv), ncol = p)) %*% tstat
    MCAL <- 0
    C <- Inf
    E <- Inf
    R <- Inf
    for (inc in 2:steps) {
        binc <- bhat[steps + 1 - inc, ]
        dinc <- (t(sx$v) %*% binc)/comp
        for (j in 1:p) {
            if (abs(dinc[j]) < omdmin) 
                dinc[j] <- 0
            if (dinc[j] < 0) 
                stop("A negative Uncorrelated Component Shrinkage Factor should not occur.")
        }
        d <- matrix(dinc, p, p)
        minc <- p - sum(dinc)
        omd <- pmax(1 - dinc, omdmin)
        ddomd <- dinc/omd
        rxi <- sum(t(arho) * sqrt(ddomd))
        slik <- 2/(rxi + sqrt(4 * n + rxi^2))
        clik <- 2 * n * log(slik) + sum(ddomd) - (rxi/slik) - 
            n * log((1 - r2)/n)
        ebay <- sum(frat * omd - log(omd))
        sr2d <- sum(dinc * rho^2)
        rcof <- -sum(log(omd)) + n * log((1 - sr2d)/(1 - r2))
        C <- rbind(C, clik)
        E <- rbind(E, ebay)
        R <- rbind(R, rcof)
        MCAL <- rbind(MCAL, minc)
        vecr <- (idty - d) %*% srat
        compr <- const * vecr %*% t(vecr) + (2 * d - idty) * 
            eiginv
        diagc <- diag(diag(compr), ncol = p)
        lowr <- eiginv * d^2
        maxd <- matrix(pmax(diagc, lowr), p, p)
        compr <- compr - diagc + maxd
        smse <- sx$v %*% compr %*% t(sx$v)
        rinc <- diag(smse)
        lowb <- sx$v %*% lowr %*% t(sx$v)
        rinc <- pmax(rinc, diag(lowb))
        tinc <- sum(rinc)
        emse <- eiginv - compr
        sfac <- min(abs(diag(emse)))/100
        if (sfac < 1e-05) 
            sfac <- 1e-05
        eign <- eigen(emse/sfac)
        einc <- rev(eign$values) * sfac
        cinc <- matrix(0, p, 1)
        if (is.na(einc[1])) 
            einc[1] <- 0
        if (einc[1] < 0) {
            eign$vectors <- sx$v %*% eign$vectors
            cinc <- eign$vectors[, p]
            if (rscale == 2) {
                cinc <- cinc %*% xscale
                cinc <- cinc/sqrt(sum(cinc^2))
            }
            if (t(cold) %*% cinc < 0) 
                cinc <- -1 * cinc
            cold <- cinc
        }
        bstar <- cbind(bstar, binc)
        risk <- cbind(risk, rinc)
        exev <- cbind(exev, einc)
        infd <- cbind(infd, cinc)
        delta <- cbind(delta, dinc)
        tsmse <- rbind(tsmse, tinc)
        mcal <- rbind(mcal, minc)
    }
    if (rscale == 2) {
        bstar <- yscale * solve(xscale) %*% bstar
        risk <- yscale^2 * solve(xscale^2) %*% risk
    }
    mlik <- cbind(MCAL, C, E, R)
    dimnames(mlik) <- list(0:(steps - 1), c("M", "CLIK", "EBAY", 
        "RCOF"))
    bstar <- t(bstar)
    risk <- t(risk)
    exev <- t(exev)
    infd <- t(infd)
    delta <- t(delta)
    sext <- cbind(tsmse, mcal)
    dimnames(sext) <- list(0:(steps - 1), c("TSMSE", "MCAL"))
    RXolist <- c(RXolist, list(lars = larsobj, coef = bstar, 
        rmse = risk, exev = exev, infd = infd, spat = delta, 
        mlik = mlik, sext = sext))
    class(RXolist) <- "RXuclars"
    RXolist
}
