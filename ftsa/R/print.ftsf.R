print.ftsf = function (x, ...) 
{
    st <- ifelse(frequency(x$mean$time) != 1, deparse(start(x$mean$time)), 
        format(tsp(x$mean$time)[1]))
    ed <- ifelse(frequency(x$mean$time) != 1, deparse(end(x$mean$time)), 
        format(tsp(x$mean$time)[2]))
    cat("Functional time series forecasts")
    cat("\n  for", st, "to", ed, "\n  using", x$coeff[[1]]$method, 
        "method\n  based on\n")
    print(x$model)
}
