fits.compare <- function (...)
{
    tmp.models <- list(...)
    j <- 0
    models <- NULL
    the.call <- sys.call()
    n <- length(the.call) - 1
    nams <- rep("", n)
    nm <- names(the.call)
    if (is.null(nm))
        nm <- rep("", n + 1)
    for (i in 1:length(tmp.models)) {
        xx <- class(tmp.models[[i]])[1]
        if (xx == "aov" || xx == "lm" || xx == "lm.robust" ||
            xx == "lm.Huber" || xx == "glm" || xx == "cubinf" || xx == "TML") {
            j <- j + 1
            if (j == 1)
                models <- list(tmp.models[[i]])
            else models[[j]] <- tmp.models[[i]]
            if (nm[i + 1] != "")
                models[[j]]$name <- nm[i + 1]
            else models[[j]]$name <- deparse(the.call[[i + 1]])
        }
        else warning("Only class lm, aov, glm, TML and robust objects are allowed as input to fits.compare().")
    }
    class(models) <- "fits.compare"
    models
}
