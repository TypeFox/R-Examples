## Try to estimate movMF via EM.

movMF <-
function(x, k, control = list(), ...)
{
    ## Be nice a la David.
    control <- c(control, list(...))

    if(missing(k)) k <- NULL
    ## Normalize data just in case.
    x <- skmeans:::row_normalize(x)

    n <- nrow(x)
    d <- ncol(x)
    s <- d / 2 - 1

    ## Control parameters.

    maxiter <- control$maxiter
    if(is.null(maxiter))
        maxiter <- 100L

    reltol <- control$reltol
    if(is.null(reltol))
        reltol <- sqrt(.Machine$double.eps)

    verbose <- control$verbose
    if(is.null(verbose))
        verbose <- getOption("verbose")

    E_methods <-
        c("softmax",
          "hardmax",
          "stochmax")
    E <- control$E
    if(is.null(E))
        E <- "softmax"
    else {
        pos <- pmatch(tolower(E), tolower(E_methods))
        if(is.na(pos))
            stop("Invalid E-step method.")
        E <- E_methods[pos]
    }
    do_P <- switch(EXPR = E,
                   "softmax" =
                   function(G, g) exp(G  - g),
                   "hardmax" =
                   function(G, g)
                   .posterior_from_ids(max.col(G), ncol(G)),
                   function(G, g)
                   t(apply(exp(G - g), 1,
                           function(prob) rmultinom(1, 1, prob))))

    nu <- if (is.null(control$nu)) 0 else control$nu
    
    kappa_solvers <-
        c("Banerjee_et_al_2005",
          "Tanabe_et_al_2007",
          "Sra_2012",
          "Song_et_al_2012",
          "uniroot",
          "Newton",
          "Halley",
          "hybrid",
          "Newton_Fourier")
    kappa <- control$kappa
    if(is.numeric(kappa)) {
        ## CASE A: Use a common given value of kappa.
        if (length(kappa) > 1)
          warning("only the first element of 'kappa' is used for a common given kappa")
        kappa_given <- kappa[1]
        use_common_kappa <- TRUE
        do_kappa <- function(norms, alpha)
            rep.int(kappa_given, length(alpha))
        df_kappa <- function(k) 0L
    } else {
        kappa <- as.list(kappa)
        ## Should a common value be used or not?
        pos <- match("common", names(kappa), nomatch = 0L)
        if(pos > 0L) {
            use_common_kappa <- identical(kappa[[pos]], TRUE)
            kappa <- kappa[-pos]
            if (length(nu) > 1)
              warning("only the first element of 'nu' is used for common kappa")
            nu <- nu[1]
        } else {
            use_common_kappa <- FALSE
        }
        if(length(kappa)) {
            ## Solver specifications.
            kname <- kappa[[1L]]
            kargs <- kappa[-1L]
            pos <- pmatch(tolower(kname), tolower(kappa_solvers))
            if(is.na(pos))
                stop("Invalid kappa solver.")
            kappa <- kappa_solvers[pos]
            kfun <- get(sprintf("solve_kappa_%s", kappa))
            solve_kappa <- function(Rbar)
                do.call(kfun, c(list(Rbar, d), kargs))
        } else {
            ## Default solver.
            solve_kappa <- function(Rbar)
                solve_kappa_Newton_Fourier(Rbar, d)
        }
        if(use_common_kappa) {
            do_kappa <- function(norms, alpha) 
                rep.int(solve_kappa(sum(norms) / (n + nu)), length(alpha))
            df_kappa <- function(k) 1L
        } else {
            do_kappa <- function(norms, alpha)
                solve_kappa(norms / (n * alpha + nu))
            df_kappa <- function(k) k
        }
    }

    ## Allow to specify "known" ids for fitting the respective movMF
    ## model by maximum likelihood.  This can be accomplished within our
    ## framework as follows:
    ## * Compute P as the corresponding membership matrix.
    ## * Perform one M step to estimate the parameters.
    ## * Perform one E step to compute the likelihood, but keep P.
    ids <- control$ids
    if(!is.null(ids)) {
        if(identical(ids, TRUE))
            ids <- attr(x, "z")         # Be nice for the case of data
                                        # simulated by rmovMF().
        if (nrow(x) != length(ids))
          stop("length of 'ids' needs to match the number of observations")
        P0 <- .posterior_from_ids(ids, k)
        do_P <- function(G, g) G
        start  <- list(P0)
        maxiter <- 1L
        if(!is.null(control$nruns))
            warning("control argument 'nruns' ignored because ids are specified")
        if(!is.null(control$start))
            warning("control argument 'start' ignored because ids are specified")
    } else {
        ## Initialization.
        start <- control$start
        nruns <- control$nruns
        if(is.null(start)) {
            if(is.null(nruns))
                nruns <- 1L
            start <- as.list(rep.int("p", nruns))
        }
        else {
            if(!is.list(start))
                start <- list(start)
            if(!is.null(nruns))
                warning("control argument 'nruns' ignored because 'start' is specified")
        }
    }
    
    nruns <- length(start)

    minalpha <- control$minalpha
    if(is.null(minalpha)) minalpha <- 0
    if(minalpha >= 1)
        minalpha <- minalpha / n

    converge <- control$converge
    if(is.null(converge)) {
        converge <- if(E == "stochmax") FALSE else TRUE
    }
    
    L_opt <- -Inf
    opt_old <- opt <- NULL

    run <- 1L
    logLiks <- vector(length = nruns)
    
    if(verbose && (nruns > 1L))
        message(gettextf("Run: %d", run),
                domain = NA)

    repeat {
        G <- NULL
        P <- movMF_init(x, k, start[[run]])
        if (!use_common_kappa) {
          if (length(nu) > 1 & length(nu) != ncol(P))
            warning("nu is changed to have length", ncol(P))
          nu <- rep_len(nu, ncol(P))
        }
        L_old <- -Inf
        iter <- 0L
        logLiks[run] <- tryCatch({
            while(iter < maxiter) {
            ## M step.
                alpha <- colMeans(P)
                while(any(alpha < minalpha)) {
                    if(verbose) 
                        message("*** Removing one component ***")
                    nok <- which.min(alpha)
                    P <- P[, -nok, drop = FALSE]
                    if (!use_common_kappa)
                      nu <- nu[-nok]
                    P <- do_P(P, log_row_sums(P))
                    alpha <- colMeans(P)
                    if(!is.null(G))
                        L_old <- sum(log_row_sums(G[, -nok, drop = FALSE]))
                }
                if(any(alpha == 0 & nu <= 0))
                    stop("Cannot handle empty components")
                M <- skmeans:::g_crossprod(P, x)
                norms <- row_norms(M)
                M <- M / ifelse(norms > 0, norms, 1)
                ## If a cluster contains only identical observations,
                ## Rbar = 1.
                kappa <- do_kappa(norms, alpha)
            
                ## E step.
                G <- cadd(skmeans:::g_tcrossprod(x, kappa * M),
                          log(alpha) -  lH(kappa, s))
                g <- log_row_sums(G)
                L_new <- sum(g)
                if(verbose && (iter %% verbose == 0))
                    message(gettextf("Iteration: %d *** L: %g",
                                     iter, L_new),
                            domain = NA)
                if(converge) {
                    if(abs(L_old - L_new) < reltol * (abs(L_old) + reltol)) {
                        L_old <- L_new     
                        break
                    }
                    L_old <- L_new     
                } else if(L_new > L_old) {
                    L_old <- L_new
                    opt_old <- .movMF_object(kappa * M, alpha, L_old, P, iter)
                }
            
                P <- do_P(G, g)
                iter <- iter + 1L
            }
            if(L_old > L_opt) {
                opt <- if(converge)
                    .movMF_object(kappa * M, alpha, L_old, P, iter)
                else opt_old
                L_opt <- L_old
            }
            L_old
        }, error = function(e) {
            if(verbose) {
                ## <NOTE>
                ## Reporting problems is a bit of a mess in cases the
                ## above explicitly called stop(), in which case the
                ## condition call is
                ##   doTryCatch(return(expr), name, parentenv, handler)
                ## Ideally, in these cases we would throw suitably
                ## classed conditions, and provide call/message methods
                ## for them.
                call <- conditionCall(e)
                msg <- conditionMessage(e)
                s <- if(!is.null(call) &&
                        (substring(s <- deparse(call)[1L], 1L, 10L) !=
                         "doTryCatch"))
                    sprintf("Error in %s: %s", s, msg)
                else
                    sprintf("Error: %s", msg)
                message(gettextf("EM algorithm did not converge:\n%s",
                                 s),
                        domain = NA)
                ## </NOTE>
            }
            NA
        })
        
        if(run >= nruns) break
        
        run <- run + 1L
        if(verbose)
            message(gettextf("Run: %d", run),
                    domain = NA)
    }

    ## Compute log-likelihood.
    if(is.null(opt))
        stop("EM algorithm did not converge for any run")
    k <- length(alpha)
    dimnames(opt$theta)<- list(seq_len(k), colnames(x))
    dimnames(opt$P)<- list(rownames(x), seq_len(k))
    ll <- L(x, opt$theta, opt$alpha)
    ## Add the "degrees of freedom" (number of (estimated) parameters in
    ## the model): with k the number of classes actually used,
    ##   \mu: k * (d - 1) (as constrained to unit length),
    ##   \kappa: 0, 1 or k depending on whether we use a common given,
    ##        a common (but not given) kappa, or individual kappas.
    ##   \alpha: k - 1    (as constrained to unit sum),
    ## for a total of
    ##   k d - 1 + df_kappa(k)
    attr(ll, "df") <- d * k - 1L + df_kappa(k)
    attr(ll, "nobs") <- n
    class(ll) <- "logLik"
    opt$ll <- ll
    opt$details <- list(reltol = reltol,
                        iter = c(iter = opt$iter, maxiter = maxiter),
                        logLiks = logLiks,
                        E = E,
                        kappa = control$kappa,
                        minalpha = minalpha,
                        converge = converge)
    opt
}

## Generator.

.movMF_object <-
function(theta, alpha, L, P, iter)
{
    o <- list(theta = theta, alpha = alpha, L = L, P = P, iter = iter)
    class(o) <- "movMF"
    o
}

## Methods.

print.movMF <-
function(x, ...)
{
    cat("theta:\n")
    print(x$theta)
    cat("alpha:\n")
    print(x$alpha)
    cat("L:\n")
    print(x$L)
    invisible(x)
}

coef.movMF <-
function(object, ...)
    object[c("theta", "alpha")]


logLik.movMF <-
function(object, newdata, ...)
{
  if (missing(newdata))
    return(object$ll)
  else {
    newdata <- skmeans:::row_normalize(newdata)
    return(L(newdata, object$theta, object$alpha))
  }
}

predict.movMF <-
function(object, newdata = NULL,
         type = c("class_ids", "memberships"), ...)
{
    type <- match.arg(type)
    P <- if(is.null(newdata)) object$P else {
        x <- skmeans:::row_normalize(newdata)
        theta <- object$theta
        alpha <- object$alpha
        kappa <- row_norms(theta)
        ## Could maybe check dimensions.
        ## Same as for E step in movMF().
        d <- nrow(theta)
        s <- d / 2 - 1 
        G <- cadd(skmeans:::g_tcrossprod(x, theta),
                  log(alpha) -  lH(kappa, s))
        g <- log_row_sums(G)
        exp(G - g)
    }
    if(type == "class_ids")
        max.col(P)
    else
        P
}

## Initializer for normalized x.

## Try something similar to what we do for skmeans, but note that
## eventually we need a matrix of posterior probabilities P.

movMF_init <-
function(x, k, start)
{
    if(is.character(start)) {
        if(start %in% c("p", "S", "s")) {
            ## <FIXME>
            ## How should we turn initial centroids into a posterior?
            ## For now, do fuzzy classification with m = 2 and cosine
            ## dissimilarity.
            M <- skmeans:::.skmeans_init_for_normalized_x(x, k, start)
            D <- pmax(1 - skmeans:::g_tcrossprod(x, M), 0)
            clue:::.memberships_from_cross_dissimilarities(D, 2)
            ## </FIXME>
        }
        else if(start == "i") {
            ## Initialize by choosing random class ids, and building a
            ## binary posterior from these.
            ids <- sample.int(k, nrow(x), replace = TRUE)
            .posterior_from_ids(ids, k)
        }
        else 
            stop("Invalid control option 'start'")
    }
    else if(!is.null(dim(start)))
        start
    else {
        ## A vector of class ids, hopefully.
        ids <- match(start, sort(unique(start)))
        .posterior_from_ids(ids, k)
    }
}

## Normalizing constant for the von Mises Fisher distribution.
## (In case the reference density on the unit sphere is 1.)

C <-
function(kappa, d)
{
    s <- d / 2 - 1
    kappa^s / ((2 * pi)^(s + 1) * besselI(kappa, s))
}

## Utility function for scaling the columns of a matrix.

cmult <-
function(A, x)
    A * rep.int(x, rep.int(nrow(A), ncol(A)))


cadd <-
function(A, x)
    A + rep.int(x, rep.int(nrow(A), ncol(A)))

## Utility for computing sums avoiding *flows.

log_row_sums <-
function(x)
{
    M <- x[cbind(seq_len(nrow(x)), max.col(x, "first"))]
    M + log(rowSums(exp(x - M)))
}

## Utility functions for estimating kappa.

## <NOTE>
## Work only for 0 <= Rbar < 1.
## We could try to additionally implement solve_kappa(1, ...) = Inf.
## But the utilities are internal only and the other functions cannot
## gracefully handle kappa = Inf anyways ...
## </NOTE>

solve_kappa_Banerjee_et_al_2005 <-
function(Rbar, d)
{
    if(any(Rbar >= 1))
        stop("Cannot handle infinite concentration parameters")
    
    Rbar * (d - Rbar ^ 2) / (1 - Rbar ^ 2)
}

solve_kappa_Tanabe_et_al_2007 <-
function(Rbar, d, c = 1, tol = 1e-6)
{
    if(any(Rbar >= 1))
        stop("Cannot handle infinite concentration parameters")

    old <- Rbar * (d - c) / (1 - Rbar ^ 2)
    repeat {
      kappa <- ifelse(old == 0, 0, old * (Rbar / A(old, d)))
      if(max(abs(kappa - old)) < tol) break
      old <- kappa
    }
    kappa
}

solve_kappa_Sra_2012 <-
function(Rbar, d)
{
    ## Initialize using the Banerjee et al approximation.
    kappa <- solve_kappa_Banerjee_et_al_2005(Rbar, d)
    ## Perform two Newton steps.
    A <- A(kappa, d)
    kappa <- kappa - (A - Rbar) / Aprime(kappa, d, A = A)
    A <- A(kappa, d)
    kappa <- kappa - (A - Rbar) / Aprime(kappa, d, A = A)
    kappa
}

solve_kappa_Halley <-
function(Rbar, d, tol = 1e-6, maxiter = 100L)
{
    if(any(Rbar >= 1))
        stop("Cannot handle infinite concentration parameters")
    n <- max(length(Rbar), length(d))
    d <- rep_len(d, n)
    Rbar <- rep_len(Rbar, n)

    sapply(seq_along(Rbar),
           function(i) {
             r <- Rbar[i]
             r2 <- r^2
             r3 <- r2 * r
             r4 <- r3 * r
             D <- d[i]
             Dm1 <- D - 1
             nu <- D / 2 - 1
             kappa_0 <- Inf
             kappa <- Rinv_lower_Amos_bound(r, nu)
             A <- A(kappa, D)
             Adiff <- A - r
             for (n in seq_len(maxiter)) {
               Adiff2 <- Adiff^2
               Adiff3 <- Adiff^3
               Dkappa <- Dm1 / kappa
               numerator <-
                 2 * (Adiff * (1 - r * Dkappa - r2) - Adiff2 * (2 * r + Dkappa) - Adiff3)
               denominator <-
                 2 * (r4 + r3 * 2 * Dkappa + r2 * (-2 + Dkappa^2) -  r * 2 * Dkappa + 1) + 
                   Adiff * (6 * r3 + 9 * r2 * Dkappa + r * (-6 + Dkappa * (3 * D - 4) / kappa) - 3 * Dkappa) + 
                     Adiff2 * (6 * r2 + r * 6 * Dkappa + Dkappa * (D - 2)/kappa - 2) + 
                       Adiff3 * (2 * r + Dkappa)
               kappa <- kappa - numerator / denominator
               A <- A(kappa, D)
               Adiff <- A - r
               if ((abs(Adiff) < tol) ||
                   (abs(kappa_0 - kappa) < tol * (kappa_0 + tol))) {
                 break
               }
               kappa_0 <- kappa
             }
             kappa
           })
}

solve_kappa_Song_et_al_2012 <-
function(Rbar, d)
{
    ## Initialize using the Banerjee et al approximation.
    kappa <- solve_kappa_Banerjee_et_al_2005(Rbar, d)
    ## Perform two Halley steps.
    A <- A(kappa, d)
    Adiff <- A - Rbar
    Aprime <- Aprime(kappa, d, A = A)
    kappa <- kappa - 2 * Adiff * Aprime /
        (2 * Aprime^2 - Adiff * Adoubleprime(kappa, d, A = A))
    A <- A(kappa, d)
    Adiff <- A - Rbar
    Aprime <- Aprime(kappa, d, A = A)
    kappa <- kappa - 2 * Adiff * Aprime /
        (2 * Aprime^2 - Adiff * Adoubleprime(kappa, d, A = A))
    kappa
}

solve_kappa_uniroot <-
function(Rbar, d, tol = 1e-6)
{
    if(any(Rbar >= 1))
        stop("Cannot handle infinite concentration parameters")
    n <- max(length(Rbar), length(d))
    d <- rep_len(d, n)
    Rbar <- rep_len(Rbar, n)
    
    sapply(seq_along(Rbar),
           function(i) {
             r <- Rbar[i]
             nu <- d[i] / 2 - 1
             interval <- c(Rinv_lower_Amos_bound(r, nu),
                           Rinv_upper_Amos_bound(r, nu))
             if (abs(diff(interval)) < tol)
               mean(interval)
             else 
               uniroot(function(kappa) A(kappa, d) - r,
                       interval = interval,
                       tol = tol)$root})
}

solve_kappa_Newton <-
function(Rbar, d, tol = 1e-6, maxiter = 100)
{
    if(any(Rbar >= 1))
        stop("Cannot handle infinite concentration parameters")
    n <- max(length(Rbar), length(d))
    d <- rep_len(d, n)
    Rbar <- rep_len(Rbar, n)

    sapply(seq_along(Rbar),
           function(i) {
             r <- Rbar[i]
             D <- d[i]
             nu <- D / 2 - 1
             kappa_0 <- Inf
             kappa <- Rinv_lower_Amos_bound(r, nu)
             for (n in seq_len(maxiter)) {
               A <- A(kappa, D)
               Adiff <- A - r
               kappa <- kappa - Adiff / Aprime(kappa, D, A = A)
               A <- A(kappa, D)
               if ((abs(Adiff) < tol) ||
                   (abs(kappa_0 - kappa) < tol * (kappa_0 + tol)))
                 break
               kappa_0 <- kappa
             }
             kappa
           })
}

solve_kappa_Newton_Fourier <-
function(Rbar, d, tol = 1e-6, maxiter = 100L)
{
    if(any(Rbar >= 1))
        stop("Cannot handle infinite concentration parameters")
    n <- max(length(Rbar), length(d))
    d <- rep_len(d, n)
    Rbar <- rep_len(Rbar, n)

    sapply(seq_along(Rbar),
           function(i) {
               r <- Rbar[i]
               D <- d[i]
               nu <- D / 2 - 1
               lower <- Rinv_lower_Amos_bound(r, nu)
               upper <- Rinv_upper_Amos_bound(r, nu)
               iter <- 1L
               while(iter <= maxiter) {
                   A <- A(lower, D)
                   Aprime <- Aprime(lower, D, A = A)
                   lower <- lower - (A - r) / Aprime
                   A <- A(upper, D)                   
                   upper <- upper - (A - r) / Aprime
                   if((upper - lower) < tol * (lower + upper)) {
                       if ((upper - lower) < - tol * (lower + upper))
                           stop("no convergence")
                       break
                   }
                   iter <- iter + 1L
               }
               ## <FIXME>
               ## What should we really return?
               (lower + upper) / 2
               ## </FIXME>
           })
}

## Combination of Newton with bisection
## See Press et al., Numerical Recipes
## Modified to
## allow different methods for determining steps (Newton, Halley)

step_Newton <-
function(x, d, Rbar, A = NULL)
{
    if (is.null(A)) 
        A <- A(x, d)
    fx <- A - Rbar
    dfx <- Aprime(x, d, A = A)
    fx / dfx
}

step_Halley <-
function(x, d, Rbar, A = NULL)
{
    if (is.null(A))
        A <- A(x, d)
    Adiff <- A - Rbar
    Adiff2 <- Adiff^2
    Adiff3 <- Adiff2 * Adiff
    dx <- (d - 1) / x
    Rbar2 <- Rbar^2
    Rbar3 <- Rbar2 * Rbar
    Rbar4 <- Rbar3 * Rbar
    
    numerator <- 2 * (Adiff * (1 - Rbar * dx - Rbar2) - Adiff2 * (2 * Rbar + dx) - Adiff3)
    denominator <-
        2 * (Rbar4 + Rbar3 * 2 * dx + Rbar2 * (-2 + dx^2) -  Rbar * 2 * dx + 1) + 
            Adiff * (6 * Rbar3 + 9 * Rbar2 * dx + Rbar * (-6 + dx * (3 * d - 4) / x) - 3 * dx) + 
                Adiff2 * (6 * Rbar2 + Rbar * 6 * dx + dx * (d - 2)/x - 2) + 
                    Adiff3 * (2 * Rbar + dx)
    numerator / denominator
}

solve_kappa_hybrid <-
function(Rbar, d, tol = 1e-6, maxiter = 100, step = step_Halley)
{
    if(any(Rbar >= 1))
        stop("Cannot handle infinite concentration parameters")
    
    n <- max(length(Rbar), length(d))
    d <- rep_len(d, n)
    Rbar <- rep_len(Rbar, n)
    
    sapply(seq_along(Rbar),
           function(i) {
             r <- Rbar[i]
             D <- d[i]
             nu <- D / 2 - 1
             rtsold <- xl <- rts <- Rinv_lower_Amos_bound(r, nu)
             xh <- Rinv_upper_Amos_bound(r, nu)
             dxold <- abs(xh - xl)
             Arts <- A(rts, D)
             for (j in seq_len(maxiter)) {
               dx <- step(rts, D, r, Arts)
               rts <- rts - dx
               if ((rts > xh) || (rts < xl) || (abs(2 * dx) > abs(dxold))) {
                 dx <- 0.5 * (xh - xl)
                 rts <- xl + dx
                 if(any(xl == rts))
                   break
               }
               else {
                 if (rtsold == rts)
                   break
               }
               if (abs(dx) < tol)
                 break
               dxold <- dx
               rtsold <- rts
               Arts <- A(rts, D)
               if ((Arts  - r) < 0.0)
                 xl <- rts
               else 
                 xh <- rts
             }
             rts
           })
  }

## Utility functions for computing A, H and its logarithm (see the
## implementation notes).
##
## With R_\nu = I_{\nu+1} / I_\nu the commonly considered modified
## Bessel function ratio and
##
##   log(H_\nu)(\kappa) = \int_0^\kappa R_\nu(t) dt
##
## (so that (log(H_\nu))' = R_\nu and H_\nu(0) = 1), we have
##
##   A_d = I_{d/2} / I_{d/2-1} = R_{d/2-1}
##
## and
##
##   log(f(x|theta)) = theta' x - log(H_{d/2-1})(\|theta\|).
##
## For R_\nu we have Amos-type bounds
##
##   G_{\alpha,\beta}(t) = t / (alpha + sqrt(t^2 + beta^2))
##
## where G_{\nu+1/2,\nu+3/2} is the uniformly best lower Amos-type
## bound, whereas there is no uniformly optimal upper bound, with
## G_{\nu,\nu+2} and G_{\nu+1/2,beta_SS(\nu)} being second-order exact
## and optimal at 0 and infinity, respectively, where
##
##   beta_SS(\nu) = sqrt((\nu+1/2)(\nu+3/2)).
##
## Inverses and anti-derivatives of G are given by
##
##   G_{\alpha,\beta}^{-1}(\rho)
##     = \frac{\rho}{1 - \rho^2}
##       (\alpha + \sqrt{\alpha^2 \rho^2 + \beta^2 (1 - \rho^2)})
##
## and
##
##   S_{\alpha,\beta}(\kappa)
##     = \sqrt{\kappa^2 + \beta^2} - \beta
##       - \alpha \log(\alpha + \sqrt{\kappa^2 + \beta^2})
##       + \alpha \log(\alpha + \beta).
##
## We use
##
##   \max(G_{\nu,\nu+2}^{-1}, G_{\nu+1/2,\beta_{SS}(\nu)}^{-1})
##
## as lower bound and approximation for R_\nu^{-1},
##
##   G_{\nu+1/2,\nu+3/2}^{-1}
##
## as upper bound for R_{\nu}^{-1}, and the antiderivative of
## \min(G_{\nu,\nu+2},G_{\nu+1/2,beta_SS(\nu)}) which is given by
##
##   S_{\nu+1/2,\beta_{SS}(\nu)}(kappa)
##     - S_{\nu+1/2,\beta_{SS}(\nu)}(\min(kappa, \kappa_\nu))
##     + S_{\nu,\nu+2}(\min(kappa, \kappa_\nu))
##
## where
##
##   \kappa_\nu = \sqrt{(3\nu + 11/2) (\nu + 3/2)}
##
## is the positive solution of
##
##   G_{\nu,\nu+2}(\kappa) = G_{\nu+1/2,beta_SS(\nu)}(\kappa)
##
## as approximation for \log(H_\nu).

G <-
function(kappa, alpha, beta)
{
    n <- max(length(kappa), length(alpha), length(beta))
    kappa <- rep_len(kappa, n)
    alpha <- rep_len(alpha, n)
    beta <- rep_len(beta, n)

    y <- numeric(n)

    ind <- ((alpha == 0) & (beta == 0))
    if(any(ind)) y[ind] <- 1
    ind <- !ind
    if(any(ind)) {
        denom <- alpha[ind] + sqrt(kappa[ind]^2 + beta[ind]^2)
        y[ind] <- ifelse(denom == 0, NaN, kappa[ind] / denom)
    }

    y
}

beta_SS <-
function(nu)
    sqrt((nu + 1/2) * (nu + 3/2))

Ginv <-
function(rho, alpha, beta)
{
    ## Perhaps recycle arguments eventually ...
    
    sigma <- rho^2
    rho * (alpha + sqrt(alpha^2 * sigma + beta^2 * (1 - sigma))) /
        (1 - sigma)
}

Rinv_lower_Amos_bound <-
function(rho, nu)
{
    ## Perhaps recycle arguments eventually ...
    pmax(Ginv(rho, nu, nu + 2),
         Ginv(rho, nu + 1/2, beta_SS(nu)))
}

Rinv_upper_Amos_bound <-
function(rho, nu)
{
    ## Perhaps recycle arguments eventually ...
    Ginv(rho, nu + 1/2, nu + 3/2)
}

S <-
function(kappa, alpha, beta)
{
    n <- max(length(kappa), length(alpha), length(beta))
    kappa <- rep_len(abs(kappa), n)
    alpha <- rep_len(alpha, n)
    beta <- rep_len(abs(beta), n)

    ab <- alpha + beta
    s <- double(n)
    ind <- (ab < 0) & (kappa^2 > alpha^2 - beta^2)
    if (any(ind))
        s[ind] <- NaN
    ind <- !ind
    if (any(ind)) {
        u <- sqrt(kappa[ind]^2 + beta[ind]^2)
        s[ind] <- u - beta[ind] - ifelse(alpha[ind] == 0, 0, alpha[ind] * log((alpha[ind] + u) / ab[ind]))
    }
    s
}

## Utility functions for computing H and log(H).

H <-
function(kappa, nu, v0 = 1)
{
    ## Compute v0 H_\nu(\kappa) by direct Taylor series summation.
    
    ## Add range checking eventually.
    n <- max(length(kappa), length(nu), length(v0))
    kappa <- rep_len(kappa, n)
    nu <- rep_len(nu, n)
    v0 <- rep_len(v0, n)
    
    .C(C_my0F1,
       as.integer(n),
       as.double(kappa ^ 2 / 4),
       as.double(nu + 1),
       as.double(v0),
       y = double(n))$y
}

lH_asymptotic <-
function(kappa, nu)
{
    ## Compute a suitable asymptotic approximation to
    ## \log(H_\nu(\kappa)).
    
    ## Add range checking eventually.    
    n <- max(length(kappa), length(nu))
    kappa <- rep_len(kappa, n)
    nu <- rep_len(nu, n)
    
    y <- double(n)
    ipk <- (kappa > 0)
    ipn <- (nu > 0)
    ind <- ipk & !ipn
    if(any(ind)) {
        ## For \log(H_0) = \log(I_0), use the asymptotic approximation
        ##   I_0(\kappa) \approx e^\kappa / \sqrt{2 \pi \kappa}
        ## (e.g., http://dlmf.nist.gov/10.40).
        y[ind] <- kappa[ind] - log(2 * pi * kappa[ind]) / 2
    }
    ind <- ipk & ipn
    if(any(ind)) {
        ## For \nu > 0, use the Amos-type based approximation discussed
        ## above.
        kappa <- kappa[ind]
        nu <- nu[ind]
        beta <- beta_SS(nu)
        kappa_L <- pmin(kappa, sqrt((3 * nu + 11 / 2) * (nu + 3 / 2)))
        y[ind] <- S(kappa, nu + 1/2, beta) +
            (S(kappa_L, nu, nu + 2) - S(kappa_L, nu + 1/2, beta))
    }
    y
}
        
lH <-
function(kappa, nu)
{
    ## Compute \log(H_\nu(\kappa)) (or an approximation to it).
    ## See the implementation notes for details.
    
    n <- max(length(kappa), length(nu))
    kappa <- rep_len(kappa, n)
    nu <- rep_len(nu, n)

    y <- lH_asymptotic(kappa, nu)
    ## If the value from the asymptotic approximation is small enough,
    ## we can use direct Taylor series summation.
    ind <- y <= 699.5
    if(any(ind))
        y[ind] <- log(H(kappa[ind], nu[ind]))
    ## For intermediate values of the asymptotic approximation, we use
    ## rescaling and direct Taylor series summation.
    ind <- !ind & (y <= 1399)
    if(any(ind)) {
        v <- y[ind] / 2
        y[ind] <- v + log(H(kappa[ind], nu[ind], exp(-v)))
    }
    ## (Otherwise, we use the asymptotic approximation.)
    y
}

## For testing this: delta should be close to 0.

delta <-
function(kappa, nu)
{
    (besselI(kappa, nu) -
     H(kappa, nu) * (kappa / 2)^nu / gamma(nu + 1))
}

## A, A prime and A double prime.

A <-
function(kappa, d, method = c("PCF", "GCF", "RH"), tol = 1e-6)
{
    method <- match.arg(method)
    n <- max(length(kappa), length(d))
    kappa <- rep_len(kappa, n)
    d <- rep_len(d, n)
    A <- vector("numeric", length = n)
    
    index <- kappa >= tol
    if (sum(index)) {
      method <- match.arg(method)
      if(method == "PCF") {
        ## Use the Perron continued fraction for R_{\nu-1}.
        ## Implementation based on Eqn 3.3' in Gautschi and Slavik
        ## (1978).
        ## Note that A_d = R_{d/2-1}.
        A[index] <- .C(C_mycfP,
                       as.integer(sum(index)),
                       as.double(kappa[index]),
                       as.double(d[index] / 2),
                       y = double(sum(index)))$y
      }
      else if(method == "GCF") {
        ## Use the Gauss continued fraction for R_{\nu-1}.
        ## Implementation based on Eqn 3.2' in Gautschi and Slavik
        ## (1978).
        ## Note that A_d = R_{d/2-1}.
        A[index] <- .C(C_mycfG,
                       as.integer(sum(index)),
                       as.double(kappa[index]),
                       as.double(d[index] / 2),
                       y = double(sum(index)))$y
      }
      else {
        ## Compute A_d via a ratio of H functions:
        ##   A_d(\kappa)
        ##     = (kappa/d) * H_{d/2}(\kappa) / H_{d/2-1}(\kappa).
        s <- d[index] / 2 - 1
        kappai <- kappa[index]
        y <- kappai / d[index]
        a <- lH_asymptotic(kappai, s + 1)
        ind <- (a <= 699.5)
        if(any(ind))
            y[ind] <- y[ind] *
                (H(kappai[ind], s + 1) /
                 H(kappai[ind], s))
        ind <- !ind & (a <= 1399)
        if(any(ind)) {
            v <- exp(- a[ind] / 2)
            y[ind] <- y[ind] *
                (H(kappai[ind], s + 1, v) /
                 H(kappai[ind], s, v))
        }
        ind <- (a > 1399)
        if(any(ind))
            y[ind] <- y[ind] *
                exp(a[ind] - lH_asymptotic(kappai[ind], s))
        if(any(y >= 1))
          stop("RH evaluation gave infeasible values which are not in the range [0, 1)")
        A[index] <- y
      }
    }
    if (sum(!index)) {
      di <- d[!index]
      kappai <- kappa[!index]
      A[!index] <- kappai / di - kappai^3 / (di^2 * (di + 2)) + 2 * kappai^5 / (di^3 * (di + 2) * (di + 4))
    }
    A 
}

Aprime <-
function(kappa, d, method = c("PCF", "GCF", "RH"), tol = 1e-6, A = NULL)
{
    n <- max(length(kappa), length(d), length(A))
    kappa <- rep_len(kappa, n)
    d <- rep_len(d, n)
    if(!is.null(A))
        A <- rep_len(A, n)
    aprime <- vector("numeric", length = n)
  
    index <- kappa >= tol
    if (sum(index)) {
        if(is.null(A)) 
            A <- A(kappa[index], d[index], method, tol)
        else
            A <- A[index]
        aprime[index] <- 1 - A ^ 2 - A * (d[index] - 1) / kappa[index]
    }
    if (sum(!index)) {
        di <- d[!index]
        aprime[!index] <- 1 / di - 3 / (di^2 * (di + 2)) * kappa[!index]^2 + 10 / (di^3 * (di + 2) * (di + 4)) * kappa[!index]^4
    }
    aprime
}

Adoubleprime <-
function(kappa, d, method = c("PCF", "GCF", "RH"), tol = 1e-6, A = NULL)
{
    n <- max(length(kappa), length(d), length(A))
    kappa <- rep_len(kappa, n)
    d <- rep_len(d, n)
    if(!is.null(A))
        A <- rep_len(A, n)
    adoubleprime <- vector("numeric", length = n)
  
    index <- kappa >= tol
    if (sum(index)) {
        di <- d[index]
        if (is.null(A)) 
            A <- A(kappa[index], di, method, tol)
        else
            A <- A[index]
        kappa2 <- kappa[index]^2
        adoubleprime[index] <- 2 * A^3 + 3 * (di - 1) / kappa[index] * A^2 + (di^2 - di - 2 * kappa2) / kappa2 * A - (di - 1) / kappa[index]
    }
    if(sum(!index)) {
        di <- d[!index]
        adoubleprime[!index] <- - 6 / (di^2 * (di + 2)) * kappa[!index] + 40 / (di^3 * (di + 2) * (di + 4)) * kappa[!index]^3
    }
    adoubleprime
}

## R, R prime and R double prime.

R <-
function(t, nu, ...)
    A(t, 2 * (nu + 1), ...)

Rprime <-
function(t, nu, ..., R = NULL)
    Aprime(t, 2 * (nu + 1), ..., A = R)

Rdoubleprime <-
function(t, nu, ..., R = NULL)
    Adoubleprime(t, 2 * (nu + 1), ..., A = R)

Rtripleprime <-
function(t, nu, ..., tol = 1e-6)
{
    n <- max(length(t), length(nu))
    t <- rep_len(t, n)
    nu <- rep_len(nu, n)
    rnu_tripleprime <- vector("numeric", length = n)
    ind <- t >= tol
    if(sum(ind)) {
        t_ind <- t[ind]
        nu_ind <- nu[ind]
        rnu <- R(t_ind, nu_ind, ..., tol = tol)
        ## No given R for now.
        rnu_prime <- Rprime(t_ind, nu_ind, tol = tol, R = rnu)
        rnu_doubleprime <- Rdoubleprime(t_ind, nu_ind, tol = tol, R = rnu)
        delta <- 2 * nu_ind + 1
        rnu_tripleprime[ind] <-
            (- 2 * delta * rnu / t_ind^3
             + 2 * delta * rnu_prime / t_ind^2
             -     delta * rnu_doubleprime / t_ind
             - 2 * (rnu_prime^2 + rnu * rnu_doubleprime))
    }
    ind <- !ind
    if(sum(ind)) {
        t_ind <- t[ind]
        nu_ind <- nu[ind]
        rnu_tripleprime[ind] <-
            (- 0.75 / ((nu_ind + 1)^2 * (nu_ind + 2))
             + 3.75 * t_ind^2 /
             ( (nu_ind + 1)^3 * (nu_ind + 2) * (nu_ind + 3) ))
    }

    rnu_tripleprime
}

## Log-likelihood

L <-
function(x, theta, alpha)
{
    sum(ldmovMF(x, theta, alpha))
}

## Utility for computing a binary posterior from the class ids.

.posterior_from_ids <-
function(ids, k = NULL)
{
    if(is.null(k))
        k <- max(ids)
    else if(max(ids) < k) {
        k <- max(ids)
        warning(gettextf("due to initialization number of components reduced to  %d",
                         k),
                domain = NA)
    }
    else if(max(ids) > k)
        stop("number of components k smaller than those provided for initialization")
    n <- length(ids)
    P <- matrix(0, n, k)
    P[cbind(seq_len(n), ids)] <- 1
    P
}
