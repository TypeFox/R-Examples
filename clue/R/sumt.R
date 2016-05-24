sumt <-
function(x0, L, P, grad_L = NULL, grad_P = NULL, method = NULL,
         eps = NULL, q = NULL, verbose = NULL, control = list())
{
    ## Default values: make it nice for others to call us.

    if(is.null(eps)) eps <- sqrt(.Machine$double.eps)
    if(is.null(method)) method <- "CG"
    if(is.null(q)) q <- 10
    if(is.null(verbose)) verbose <- getOption("verbose")
    
    Phi <- function(rho, x) L(x) + rho * P(x)

    if(is.null(grad_L) || is.null(grad_P)) {
        make_Phi <- function(rho) { function(x) Phi(rho, x) }
        make_grad_Phi <- function(rho) NULL
    }
    else {
        grad_Phi <- function(rho, x) grad_L(x) + rho * grad_P(x)
        make_Phi <- if(method == "nlm") {
            function(rho) {
                function(x)
                .structure(Phi(rho, x), gradient = grad_Phi(rho, x))
            }
        }
        else
            function(rho) {
                function(x)
                Phi(rho, x)
            }
        make_grad_Phi <- function(rho) { function(x) grad_Phi(rho, x) }
    }

    ## <NOTE>
    ## For the penalized minimization, the Newton-type nlm() may be
    ## computationally infeasible (although it works much faster for
    ## fitting ultrametrics to the Phonemes data).
    ## De Soete recommends using Conjugate Gradients.
    ## We provide a simple choice: by default, optim(method = "CG") is
    ## used.  If method is non-null and not "nlm", we use optim() with
    ## this method.  In both cases, control gives the control parameters
    ## for optim().
    ## If method is "nlm", nlm() is used, in which case control is
    ## ignored.  Note that we call nlm() with checking analyticals
    ## turned off, as in some cases (e.g. when fitting ultrametrics) the
    ## penalty function is not even continuous ...
    optimize_with_penalty <- if(method == "nlm")
        function(rho, x)
            nlm(make_Phi(rho), x, check.analyticals = FALSE) $ estimate
    else {
        function(rho, x)
            optim(x, make_Phi(rho), gr = make_grad_Phi(rho),
                  method = method, control = control) $ par
    }
    ## Note also that currently we do not check whether optimization was
    ## "successful" ...
    ## </NOTE>

    ## We currently require that x0 be a *list* of start values, the
    ## length of which gives the number of SUMT runs.  But as always,
    ## let's be nice to users and developers, just in case ...
    if(!is.list(x0))
        x0 <- list(x0)

    v_opt <- Inf
    x_opt <- NULL
    rho_opt <- NULL
    for(run in seq_along(x0)) {
        if(verbose)
            message(gettextf("SUMT run: %d", run))            
        x <- x0[[run]]
        ## <TODO>
        ## Better upper/lower bounds for rho?
        rho <- max(L(x), 0.00001) / max(P(x), 0.00001)
        ## </TODO>
        if(verbose)
            message(gettextf("Iteration: 0 Rho: %g P: %g", rho, P(x)))
        iter <- 1L
        repeat {
            ## <TODO>
            ## Shouldnt't we also have maxiter, just in case ...?
            ## </TODO>
            if(verbose)
                message(gettextf("Iteration: %d Rho: %g P: %g",
                                 iter, rho, P(x)))
            x_old <- x
            x <- optimize_with_penalty(rho, x)
            if(max(abs(x_old - x)) < eps)
                break
            iter <- iter + 1L
            rho <- q * rho
        }
        v <- Phi(rho, x)
        if(v < v_opt) {
            v_opt <- v
            x_opt <- x
            rho_opt <- rho
        }
        if(verbose)
            message(gettextf("Minimum: %g", v_opt))
    }

    .structure(list(x = x_opt, L = L(x_opt), P = P(x_opt), rho = rho_opt,
                    call = match.call()),
               class = "sumt")
}
