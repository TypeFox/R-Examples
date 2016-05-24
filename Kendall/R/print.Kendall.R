"print.Kendall" <-
function(x, ...)
{
    cat(paste(
   "tau = ", format(x$tau, digits = 3),    
    ", 2-sided pvalue =", format.pval( x$sl), sep = ""), 
    fill = TRUE)
}

