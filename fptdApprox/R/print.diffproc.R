print.diffproc <-
function (x, ...) 
{
    if (!is.diffproc(x)) 
        stop(paste(sQuote("x"), "is not of class", shQuote("diffproc")))
    cat("\nAn object of class", shQuote("diffproc"), "which represents a diffusion process, X(t), with:")
    cat("\n   Infinitesimal mean, A1(x,t) =", x$mean)
    cat("\n   Infinitesimal variance, A2(x,t) =", x$var)
    cat("\n   Transition probability density function, f(x,t|y,s) =", 
        x$tpdf)
    cat("\n   Transition probability distribution function, F(x,t|y,s) =", 
        x$tpdF)
    cat("\n\n")
}
