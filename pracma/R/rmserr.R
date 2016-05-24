##
##  r m s e r r o r . R  RMS Error
##


rmserr <- function(x, y, summary = FALSE) {
    if (!is.numeric(x) || !is.numeric(y))
        stop("Arguments 'x' and 'y' must be numeric vectors.")
    if (length(x) != length(y))
        stop("Vectors 'x' and ' y' must have the same length.")

    n <- length(x);

    mae    <- sum(abs(y - x))/n
    mae_f  <- formatC(mae, digits=4, format="f")
    mse    <- sum((y - x)^2)/n
    mse_f  <- formatC(mse, digits=4, format="f")
    rmse   <- sqrt(sum((y - x)^2)/n)
    rmse_f <- formatC(rmse, digits=4, format="f")
    mape   <- sum( abs((y - x)/x) )/n
    mape_f <- formatC(mape, digits=4, format="f")
    nmse   <- sum((y - x)^2)/sum((x - mean(x))^2)
    nmse_f <- formatC(nmse, digits=4, format="f")
    rstd   <- sqrt(sum((y - x)^2)/n) / mean(x)
    rstd_f <- formatC(rstd, digits=4, format="f")

    if (summary) {
cat("-- Error Terms --------------------------------------------------\n");
cat(" MAE:  ", mae_f,  "  \t- mean absolute error (in range [", range(x), "])\n")
cat(" MSE:  ", mse_f,  "  \t- mean squared error (the variance?!)\n");
cat(" RMSE: ", rmse_f, "  \t- root mean squared error (std. dev.)\n");
cat(" MAPE: ", mape_f, "  \t- mean absolute percentage error\n");
cat(" LMSE: ", nmse_f, "  \t- normalized mean squared error\n");
cat(" rSTD: ", rstd_f, "  \t- relative standard deviation (", mean(x), ")\n");
cat("-----------------------------------------------------------------\n");
    }

    R <- list(mae  = mae,  mse  = mse,  rmse = rmse,
              mape = mape, nmse = nmse, rstd = rstd)

    if (summary) {
        invisible(R)
    } else {
        return(R)
    }
}
