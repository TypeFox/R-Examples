print.glmpath.cr <-
function (x, ...) 
{
    cat("Call:\n")
    dput(x$call)
    actions <- x$actions
    xn <- x$xnames[-1]
    k <- length(actions)
    for (i in 1:k) {
        if (length(actions[[i]]) > 0) {
            cat("Step", i, ":")
            for (ii in actions[[i]]) {
                if (ii > 0) 
                  cat(" ", xn[ii])
                else cat(" -", xn[-ii])
            }
            cat("\n")
        }
    }
}
