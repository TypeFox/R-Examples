#' Mean, Standard deviation, Confident Intervals
#'
#' Get the Mean, Standard deviation and Confidence Intervals of a sampling trace
#'
#' @param df a dafa frame, each column is a sampling trace.
#' @param sig significance level, defaul is 0.05
#' @return a data frame with MSCI of each trace
#' @export
#' @author Li Qinglong <liqinglong0830@@163.com>

msci <- function(df, sig = 0.05)
{
    m = apply(df, 2, mean)
    s = apply(df, 2, sd)
    ci = apply(df, 2, quantile, probs = c(sig, 1 - sig))
    return(data.frame(rbind(m, s, ci), row.names = c("Mean", "Std.Dev", "5%", "95%")))
}