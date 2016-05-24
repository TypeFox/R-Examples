optimMethods <- eval(formals(optim)$method)
## nlminb(): "direct", using "gradient", also using "Hessian"
nlminbMethods <- c("nlminb", "nlminb_1D","nlminb_2D")


### TODO:
###
### 0) Check nlminb() methods; add an environment() solution for gradient !
### 1) look at speed improvement of   nlm(..., check.analyticals = FALSE)
### *) Hessian both for  nlm() and nlminb()


L1median <- function(X, m.init = apply(X, 2, median), weights = NULL,
		     method = c("nlm", "HoCrJo", "VardiZhang",
				optimMethods , nlminbMethods),
		     pscale = apply(abs(centr(X, m.init)), 2, mean, trim = 0.40),
		     tol = 1e-08, maxit = 200, trace = FALSE,
		     zero.tol = 1e-15, ...)
{
    ## Purpose: Compute the multivariate L1-median,  by different algorithms,
    ## -------  notably using an R-builtin optimizer such as nlms(), optim(),
    ##  or  nlminb() [not yet]
    ## ----------------------------------------------------------------------
    ## Arguments: X  : [n x p] data matrix
    ##        tol    : the convergence criterium:
    ##                 the iterative process stops when ||m_k - m_{k+1}|| < tol
    ##        init.m :
    ##        pscale : coordinate-wise scale (typical size of delta{m[j]})
    ##        ...    : optional arguments to nlm() or the control list args.
    ##		       of optim(), nlminb().
    ## ----------------------------------------------------------------------
    ## Author: Martin Maechler, 5 Dec 2005
    ##         For method = "HoCrJo": originally Kristel Joossens, see below

    ## slightly faster version of 'sweep(x, 2, m)':
    centr <- function(X,m) X - rep(m, each = n)

    d <- dim(X); n <- d[1]; p <- d[2]
    if(!is.numeric(m.init) || length(m.init) != p)
        stop(gettextf("'m.init' must be numeric of length p = %d", p))

    method <- match.arg(method)

    use.wts <- !is.null(weights)
    if(use.wts) {
        stopifnot(is.numeric(weights), length(weights) == n,
                  all(weights > 0))
        ## TODO: weights == 0 --> leave away those weights and observations

        ## FIXME: *most* methods, in particular the minimizer ones should work with weights!
        if(!(method %in% c("VardiZhang")))
           stop(gettextf("'weights' are not yet implemented for method '%s'",
                         method))
    }
    wts <- if(is.null(weights)) 1 else weights

    if(method != "HoCrJo") {
        stopifnot(length(pscale) == p, all(pscale > 0))
        ## Computes objective function in m based on X
        ## {automatically taken from the environment!}
        if(trace) cat("pscale=", paste(formatC(pscale),collapse=","),"\n")
    }

    ## FIXME: add *weights*
    L1obj <- function(m) sum(sqrt(rowSums(centr(X,m)^2)))

    if(method == "nlm") { ## TODO: also allow Hessian ..
        L1f.g <- function(m) {
            X. <- centr(X,m)
            d.i <- sqrt(rowSums(X.^2))
            f <- sum(d.i)
            attr(f, "gradient") <- - colSums(X. / d.i)
            f
        }
        r <- nlm(p = m.init, f = L1f.g,
                 gradtol = tol, print.level = trace,
                 steptol = tol, # <- not sure about that (takes more iter)
                 typsize = pscale, iterlim = maxit, ...)
    }
    else if(method %in% optimMethods) {

        ## optim without gradient {well, now try *with* gradient;
        ## the nice feature of nlm() above is that function and gradient
        ## *share* most computations
        r <- optim(par = m.init, fn = L1obj, method = method,
                   control = list(parscale = pscale, maxit = maxit,
                   trace = trace, reltol = tol, ...))

    }
    else if(method %in% nlminbMethods) {

        Ctrl <- list(iter.max = maxit, rel.tol = tol,
                     ## abs.tol = ???
                     ## sing.tol = ???
                     trace = trace,# < not sure if PORT's notion match ours
                     ...)
        r <- if(method == "nlminb") { # no gradient, no Hessian
            nlminb(start = m.init, objective = L1obj,
                   scale = pscale, control = Ctrl)
                                        # lower = -Inf, upper = Inf
        } else { ## f + Gradient ( + Hessian , maybe )
            ## Idea: f(), its grad() and maybe Hess()
            ##       should  *share* an environment, and hence become
            ## similarly efficient as nlm() !

            stop("nlminb with Gradient -- not yet implemented")
### FIXME
            Fenv <- new.env()
            L1f <- function(m) {
                X. <- centr(X,m)
                d.i <- sqrt(rowSums(X.^2))
                f <- sum(d.i)
                ## attr(f, "gradient") <- - colSums(X. / d.i)
                ##  f
            }

            nlminb(..., ..., ...)
        }

    }
    else if(method == "HoCrJo") {
        ## original code by Kristel Joossens, still available at
        ## http://www.econ.kuleuven.be/public/NDBAE06/programs/pca/L1median.r.txt
        ## (MM: ~/R/MM/STATISTICS/robust/Croux+Kristel/Programs/pca/L1median.r.txt
        ##   -- 2005-12-07)
        ##  with improvwments by MM
        ##     - Fix typo replace 'mX' by 'm'
        ##     - s/maxstep/maxit/  max. #{iterations}

        ## MM: Note however: the 'tol' is not used as it should be __at all__
        ## --  (it's used only for the stephalfing, but not for the main
        ##      convergence check..)

        ##
        ## Ref: Hossjer and Croux (1995)
        ##  "Generalizing Univariate Signed Rank Statistics for Testing
        ##   and Estimating a Multivariate Location Parameter";
        ##   Non-parametric Statistics, 4, 293-308.

        m <- m.init
        k <- 1
        if(trace) nstps <- 0
        while (k <= maxit) {
            mold <- m
            obj.old <- if(k == 1) L1obj(mold) else obj
            X. <- centr(X, m)
            Xnorms <- sqrt(rowSums(X. ^ 2))
            inorms <- order(Xnorms)
            dx <- Xnorms[inorms] # smallest first, i.e., 0's if there are
            X  <- X [inorms,]
            X. <- X.[inorms,]
            ## using 1/D weights {gradient!} when D > 0
            w <- ## (0 norm -> 0 weight) :
                if (all(dn0 <- dx != 0))  1/dx
                else c(rep.int(0, length(dx)- sum(dn0)), 1/dx[dn0])
            delta <- colSums(X. * rep(w,p)) / sum(w)
            nd <- sqrt(sum(delta^2))

            maxhalf <- if (nd < tol) 0 else ceiling(log2(nd/tol))
            m <- mold + delta          # computation of a new estimate
            ## If step 'delta' is too far, we try halving the stepsize
            nstep <- 0
            while ((obj <- L1obj(m)) >= obj.old && nstep <= maxhalf) {
                nstep <- nstep+1
                m <- mold + delta/(2^nstep)
            }
            if(trace) {
                if(trace >= 2)
                    cat(sprintf("k=%3d obj=%19.12g m=(",k,obj),
                        paste(formatC(m),collapse=","),
                        ")", if(nstep) sprintf(" nstep=%2d halvings",nstep) else "",
                        "\n", sep="")
                nstps[k] <- nstep
            }
            if (nstep > maxhalf) { ## step halving failed; keep old
                m <- mold
                ## warning("step halving failed in ", maxhalf, " steps")
                break
            }
            k <- k+1
        }
        if (k > maxit) warning("iterations did not converge in ", maxit, " steps")
        if(trace == 1)
            cat("needed", k, "iterations with a total of",
                sum(nstps), "stepsize halvings\n")

        r <- m

    } else if(method == "VardiZhang") {

        ## *Weighted* L1-median  -- always converging --
        ##
        ## Vardi, Y. and Zhang, C.-H.  (2000).
        ## The multivariate $L_1$-median and associated data depth,
        ## \emph{Proc. National Academy of Science} {\bf 97}(4): 1423--1426.

        ## --------------------------------------------------------------------
        ## Author: Martin Maechler, Date: 16 Dec 2005, 18:21

        ## NB: originally ~/R/MM/STATISTICS/robust/L1median-Vardi.R
        ##  function(X, weights = NULL, m.init = apply(X, 2, median),
        ##           tol = 1e-15, maxit = 2000, trace = FALSE,
        ##           zero.tol = 1e-15, ...)
        ## i.e. *smaller* tol and larger maxit

        if(!is.null(weights))
            stopifnot(is.numeric(weights), length(weights) == n,
                      all(weights > 0))
        ## TODO: weights == 0 --> leave away those weights and observations
        wts <- if(is.null(weights)) 1 else weights

        T.t <- function(y) { ## T~() - the original iterator function  [2.3]
            Id <- wts/sqrt(rowSums(centr(X, y) ^ 2))
            Id[!is.finite(Id)] <- 0     # where denominator was 0
            colSums(X * Id) / sum(Id)
        }

        Tnew <- function(y) { ## T() - the new iterator function      [2.6]
            X. <- centr(X, y)
            di.y <- sqrt(rowSums(X. ^ 2))
            ## y == X[i,] -- numerically: smallest is "almost eq"
            i.min <- which.min(di.y)
            ## FIXME: the RHS 'zero.tol * median(di.y)' should be set initially!
            if(di.y[i.min] < zero.tol * median(di.y)) { ## y "==" x_[i.min]
                if(trace) {
                    if(trace >= 2) cat(" y == x_",i.min," ", sep="")
                    eqs <- eqs+1
                }
                eta.y <- if(is.null(weights)) 1 else weights[i.min]
                Id <- (wts/di.y)[- i.min]
                R.t <- colSums(X.[- i.min, ,drop=FALSE] * Id) ## R.tilde
                T.t <- colSums(X [- i.min, ,drop=FALSE] * Id) / sum(Id) ## T.tilde
                e.div.r <- eta.y / sqrt(sum(R.t ^ 2))
                max(0, 1 - e.div.r) * T.t + min(1, e.div.r) * y
            }
            else { ## y is different from all data points
                Id <- wts/di.y
                colSums(X * Id) / sum(Id)
            }
        }

        m <- m.init
        k <- 1
        if(trace) eqs <- 0
        while (k <= maxit) {
            mold <- m
            m <- Tnew(m)
            if((s1 <- sum(abs(mold - m))) < tol * (s2 <- sum(abs(m))))
                break
            if(trace) {
                if(trace >= 2)
                    cat(sprintf("k=%3d rel.chg=%12.7g, m=(",k, s1/s2),
                        paste(formatC(m, digits=5),collapse=", "), ")\n", sep="")
                else cat(".")
            }
            k <- k+1
        }
        if (k > maxit) warning("iterations did not converge in ", maxit, " steps")
        if(trace)
            cat(" needed ", k, " iterations (", eqs, " of which had y==x_k)\n",
                sep="")

        r <- m

    } ## end{ VardiZhang }

### FIXME: use common result structure !

    r
}
