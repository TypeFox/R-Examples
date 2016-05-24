multiStart <- function(par, fn, gr = NULL, action = c("solve", "optimize"), 
    method = c(2, 3, 1),   lower=-Inf, upper=Inf,
    project=NULL, projectArgs=NULL, 
    control = list(), quiet=FALSE, details = FALSE, ...) 
    {
    if (is.null(dim(par))) 
        par <- matrix(par, nrow = 1, ncol = length(par))

    dtls <- list()
    cvg <- vector("logical", length = nrow(par))
    values <- rep(NA, length = nrow(par))  
    action <- match.arg(action)
    feval <- iter <- 0
    pmat <- matrix(NA, nrow(par), ncol(par))

    for (k in 1:nrow(par)) {
        if (!quiet) cat("Parameter set : ", k, "... \n")
        if (action == "solve") 
            ans <- try(BBsolve(par[k, ], fn = fn, method = method, 
                control = control, quiet=quiet, ...), silent=TRUE)
        if (action == "optimize") 
            ans <- try(BBoptim(par[k, ], fn = fn, gr = gr, 
	        method = method, lower=lower, upper=upper,
                project=project, projectArgs=projectArgs,
		control = control,  quiet=quiet, ...), silent=TRUE)
        if (inherits(ans, "try-error")) next
        cvg[k] <- (ans$convergence == 0)
        values[k] <- if (action == "solve") ans$residual else ans$value
        pmat[k, ] <- ans$par
        if (details) dtls <- append(dtls, ans)
        }
    ans.ret <- list(par = pmat, fvalue = values, converged = cvg)
    if (details) attr(ans.ret, "details") <- dtls
    ans.ret
    }

