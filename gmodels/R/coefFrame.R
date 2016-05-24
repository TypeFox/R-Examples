# $Id: coefFrame.R 625 2005-06-09 14:20:30Z nj7w $

coefFrame <- 
function (mod, data, by = NULL, fit.on = TRUE, fitfun, keep.unused.levels = TRUE, 
    byvar.sep = "\001" , ...) 
{
    fit.on <- eval(substitute(fit.on), data, parent.frame())
    out <- frameApply(data, on = intersect(all.vars(mod), names(data)), 
        by = by, subset = fit.on, byvar.sep = byvar.sep, fun = function(sub.dat, 
            ...) {
            fit <- try(fitfun(mod, data = sub.dat, ...), silent = TRUE)
            if (inherits(fit, "try-error")) 
                return(fit)
            outi <- coef(fit)
            outi
        }, ...)
    if (keep.unused.levels) 
        out <- unique(merge(data[by], out, all.x = TRUE))
    out
}


