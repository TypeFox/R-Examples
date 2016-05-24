Const <- function(name = c("pi", "gamma", "catalan", "log2"), prec = 120L,
                  rnd.mode = c('N','D','U','Z','A'))
{
    stopifnot(is.numeric(prec))
    if(is.na(i <- pmatch(name, eval(formals()$name))))
        stop("'name' must be one of ",
             paste(paste("'",eval(formals()$name),"'",sep=""),
                   collapse=", "))
    new("mpfr", list(.Call(const_asMpfr, i, prec, match.arg(rnd.mode))))
}
## fails here; must happen *after* dyn.load, i.e. in
## ./zzz.R : Pi <- Const("pi")
