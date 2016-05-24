### Tools for computing the worst VaR_alpha for given margins ##################

### 1) Crude VaR bounds (for both best and worst VaR) ##########################

##' @title Crude bounds for any VaR_alpha
##' @param alpha confidence level
##' @param qF (list of) marginal quantile functions
##' @param ... ellipsis argument passed to qF()
##' @return 2-vector containing crude VaR_alpha bounds
##' @author Marius Hofert
crude_VaR_bounds <- function(alpha, qF, ...)
{
    ## ... are passed to *all* qF()
    if(!is.list(qF))
        stop("qF has to be a list of (quantile) functions")
    d <- length(qF)
    qF.low <- sapply(qF, function(qF.) qF.(alpha/d, ...))
    qF.up  <- sapply(qF, function(qF.) qF.((d-1+alpha)/d, ...))
    d * c(min(qF.low), max(qF.up))
}


### 2) Explicit worst VaR in the homogeneous case ##############################

### Dual bound #################################################################

##' @title D(s,t) = d \int_{t}^{s-(d-1)t} \bar{F}(x) dx / (s-dt)
##' @param s real number
##' @param t real number < s/d
##' @param d dimension (integer > 2)
##' @param pF marginal distribution function (same for all d)
##' @param ... ellipsis argument passed to integrate()
##' @return D(s,t)
##' @author Marius Hofert
##' @note If t -> s/d-, l'Hospital's Rule shows that D(s, s/d) = d\bar{F}(s/d)
dual_bound_2 <- function(s, t, d, pF, ...)
{
    stopifnot(length(t) == 1)
    if(t > s/d) stop("t must be <= s/d")
    ## use D(s,t) = d( 1-\int_{t}^{s-(d-1)t} F(x) dx/(s-d*t) ) in this case
    if(t == s/d) d*(1-pF(s/d)) else
    d * (1 - (1/(s-d*t)) * integrate(pF, lower=t, upper=s-(d-1)*t, ...)$value)
}

##' @title Auxiliary function \bar{F}(t) + (d-1) * \bar{F}(s-(d-1)*t)
##' @param s real number
##' @param t real number < s/d
##' @param d dimension (integer > 2)
##' @param pF marginal distribution function (same for all d)
##' @return \bar{F}(t) + (d-1) * \bar{F}(s-(d-1)*t)
##' @author Marius Hofert
dual_bound_2_deriv_term <- function(s, t, d, pF)
    1-pF(t) + (d-1)*(1-pF(s-(d-1)*t))

##' @title Dual bound D(s)
##' @param s real number
##' @param d dimension (integer > 2)
##' @param pF marginal distribution function (same for all d)
##' @param ... ellipsis argument passed to dual_bound_2()'s integrate()
##' @return D(s)
##' @author Marius Hofert
##' @note The "first-order condition" (second equality in (14) in 2)) comes from the
##'       fact that
##'       (d/dt) D(s,t) = [ (-d)[\bar{F}(s-(d-1)t)(d-1)+\bar{F}(t)](s-dt) +
##'                         d^2 \int_{t}^{s-(d-1)t} \bar{F}(x) dx ] / (s-dt)^2 = 0
##'       if and only if
##'       d (\int_{t}^{s-(d-1)t} \bar{F}(x) dx) / (s-dt) = \bar{F}(s-(d-1)t)(d-1)-\bar{F}(t)
##'       => solving d (\int_{t}^{s-(d-1)t} \bar{F}(x) dx) / (s-dt) -
##'                  (\bar{F}(s-(d-1)t)(d-1)-\bar{F}(t)) = 0
##'          as a function in t for sufficiently large s leads to D(s)
dual_bound <- function(s, d, pF, tol=.Machine$double.eps^0.25, ...)
{
    stopifnot(length(s) == 1, s >= 0)
    if(s > 0) {
        ## h(s, t)
        h <- function(t) dual_bound_2(s, t=t, d=d, pF=pF, ...) -
            dual_bound_2_deriv_term(s, t=t, d=d, pF=pF)
        ## Note: h(t) -> 0 for h -> s/d- which is bad for uniroot() as the
        ##       latter will simply stop with the root t=s/d => we thus set f.upper > 0
        h.up <- -h(0) # guarantee that uniroot() doesn't fail due to root s/d
        t. <- uniroot(h, interval=c(0, s/d), f.upper=h.up, tol=tol)$root # optimal t in Equation (12) [= arginf]
        dual_bound_2_deriv_term(s, t=t., d=d, pF=pF) # dual bound D(s) in Equation (12) [= inf]
    } else {
        ## If s = 0, then t in [0, s/d] requires t to be 0 *and* f(0) = 0, so
        ## 0 is a root (as s/d). Furthermore, at t=0 (and with s=0),
        ## dual_bound_2_deriv_term(...) = d
        d
    }
}


### Wang's methods #############################################################

##' @title Scaled right-hand side term in the objective function for computing
##'        worst VaR as in McNeil, Frey, Embrechts (2015, Prop. 8.32)
##' @param c evaluation point
##' @param alpha confidence level alpha
##' @param d dimension d
##' @param method character string giving the method
##'        generic = numerical integration; Wang.Par = Pareto distibution
##' @param ... ellipsis argument containing theta (for method="Wang.Par")
##'        or qF (for method="generic")
##' @return Right-hand side term in Prop. 3.1
##' @author Marius Hofert
##' @note for the correct 'c', this is the conditional expectation
Wang_h_aux <- function(c, alpha, d, method=c("generic", "Wang.Par"), ...)
{
    ddd <- list(...)
    method <- match.arg(method)
    switch(method,
    "generic" = {
        qF <- ddd$qF # needs 'qF()'
        a <- alpha + (d-1)*c
        b <- 1-c
        qF(a)*(d-1)/d + qF(b)/d
    },
    "Wang.Par" = {
        ## We don't use qF(a)*(d-1)/d + qF(b)/d for qF(x) = qPar(x, theta=theta)
        ## here as qF(b) = qF(1-c) and 1-c==1 for small c > 0 => numerically,
        ## qF(b) = Inf then.
        th <- ddd$theta # needs 'theta'
        t1 <- (1-alpha)/c-(d-1)
        (c^(-1/th)/d) * ((d-1)*t1^(-1/th) + 1) - 1 # checked (= qF(a)*(d-1)/d + qF(b)/d)
    },
    stop("Wrong method"))
}

##' @title Objective function for computing the worst VaR as in
##'        McNeil, Frey, Embrechts (2015, Prop. 8.32)
##' @param c evaluation point
##' @param alpha confidence level alpha
##' @param d dimension d
##' @param method character string giving the method
##' @param ... ellipsis argument passed to Wang_h_aux() and integrate()
##' @return objective function for computing the worst VaR
##' @author Marius Hofert
Wang_h <- function(c, alpha, d, method=c("generic", "Wang.Par"), ...)
{
    stopifnot(0 <= c, c <= (1-alpha)/d) # sanity check (otherwise b > a)
    method <- match.arg(method)
    ddd <- list(...)

    ## Compute \bar{I}(a, b) = 1/(b-a)\int_a^b qF(y) dy =(subs) IE[L|L\in [qF(a), aF(b)]]
    Ibar <- switch(method,
    "generic" = {
        qF <- ddd$qF # needs 'qF()'
        if(c == (1-alpha)/d) { # Properly deal with limit c=(1-alpha)/d
            qF(1-(1-alpha)/d)
        } else {
            a <- alpha + (d-1)*c
            b <- 1-c
            ddd$qF <- NULL # remove from '...'
            int <- function(...)
                integrate(qF, lower=a, upper=b, ...)$value / (b-a)
            do.call(int, ddd)
        }
    },
    "Wang.Par" = {
        th <- ddd$theta # needs 'theta'
        if(c == (1-alpha)/d) { # Properly deal with limit c=(1-alpha)/d
            ((1-alpha)/d)^(-1/th) - 1
        } else {
            t1 <- (1-alpha)/c-(d-1)
            t2 <- 1-alpha-d*c
            if(th == 1) log(t1)/t2 - 1
            else (th/(1-th))*c^(1-1/th)*(1-t1^(1-1/th))/t2 - 1
        }
    },
    stop("Wrong method"))

    ## Return
    Ibar - Wang_h_aux(c, alpha=alpha, d=d, method=method, ...)
}


### Main wrapper function for computing the best/worst VaR in the homogeneous case

## Assumptions:
## - d=2: ultimately decreasing density (for x >= x0), alpha >= F(x0)
## - "Wang": F needs to live on [0, Inf), admitting a positive density which is
##           ultimately decreasing (for x >= x0), alpha >= F(x0)
## - "dual": F needs to be continuous with unbounded support and and ultimately
##           decreasing density, F(0) = 0 (otherwise, 0 as a lower bound for
##           uniroot() in dual_bound() is not valid)

##' @title Compute the best/worst VaR_\alpha in the homogeneous case with:
##'        1) d=2: Embrechts, Puccetti, Rueschendorf (2013, Proposition 2)
##'        2) d>=3:
##'           "Wang": McNeil, Frey, Embrechts (2015, Prop. 8.32)
##'                   Integral evaluated numerically; needs smaller default
##'                   tolerance for uniroot()!
##'           "Wang.Par": The same, just with explicit formula for the integral
##'                       in the Pareto case. Note that this requires to
##'                       extend the (theoretically correct) initial interval and
##'                       a smaller tolerance (see vignette).
##'           "dual": Embrechts, Puccetti, Rueschendorf (2013, Proposition 4)
##'                   Numerically less stable; no formula for best VaR known (=> NA)
##' @param alpha confidence level
##' @param d dimension
##' @param method character string giving the method
##' @param interval initial interval
##' @param tol uniroot() x-tolerance
##' @param ... ellipsis arguments passed to Wang_h()
##' @return (best VaR, worst VaR) in the homogeneous case
##' @author Marius Hofert
##' @note (*) Typos:
##'       - Wang, Peng, Yang (2013): best VaR wrong
##'       - Published version of Embrechts, Puccetti, Rueschendorf, Wang, Beleraj (2014, Prop. 3.1):
##'         Best VaR and worst VaR formulas wrong
##'       - Updated version of Embrechts, Puccetti, Rueschendorf, Wang, Beleraj (2014)
##'         on Ruodu's website: Correct worst VaR but still wrong best VaR (Eq. (3.4)).
##'       - Both best and worst VaR are correct in McNeil, Frey, Embrechts (2015, Prop. 8.32)
VaR_bounds_hom <- function(alpha, d, method=c("Wang", "Wang.Par", "dual"),
                           interval=NULL, tol=NULL, ...)
{
    stopifnot(0<alpha, alpha<1, d>=2)
    method <- match.arg(method)

    ## Deal with d==2 first ####################################################

    if(d==2) { # See Embrechts, Puccetti, Rueschendorf (2013, Prop. 2)
        if(method == "Wang.Par") {
            theta <- NULL # make CRAN check happy
            if(!hasArg(theta))
                stop("The Pareto case requires the parameter theta")
            th <- list(...)$theta
            stopifnot(length(th) == 1, th > 0) # check theta here
            qF <- function(p) qPar(p, theta=th)
            return( c((1-alpha)^(-1/th)-1, 2*(((1-alpha)/2)^(-1/th)-1)) )
        } else {
            qF <- NULL # make CRAN check happy
            if(!hasArg(qF))
                stop("The quantile function qF of F is required")
            qF <- list(...)$qF
            return(c(qF(alpha), 2*qF((1+alpha)/2)))
        }
    }

    ## Best VaR for d >= 3 #####################################################

    best <-
        switch(method,
               "Wang" = {
                   qF <- NULL # make CRAN check happy
                   if(!hasArg(qF))
                       stop("Method 'Wang' requires the quantile function qF of F")
                   ddd <- list(...)
                   qF <- ddd$qF # get qF()
                   ddd$qF <- NULL # remove from '...'
                   int <- function(...)
                       integrate(qF, lower=0, upper=alpha, ...)$value / alpha
                   max((d-1)*qF(0)+qF(alpha), # See (*) above
                       d * do.call(int, ddd))
               },
               "Wang.Par" = {
                   theta <- NULL # make CRAN check happy
                   if(!hasArg(theta))
                       stop("Method 'Wang.Par' requires the parameter theta")
                   th <- list(...)$theta
                   stopifnot(length(th) == 1, th > 0) # check theta here
                   Ibar <- if(th == 1) {
                       -log1p(-alpha) - alpha
                   } else {
                       ((1-alpha)^(1-1/th)-1)/(1-1/th) - alpha
                   }
                   max((d-1)*0 + (1-alpha)^(-1/th)-1, # See (*) above
                       d * Ibar)
               },
               "dual" = { # "dual" only provides worst VaR
                   NA
               },
               stop("Wrong method"))

    ## Worst VaR for d >= 3  ###################################################

    if(is.null(tol)) # use smaller tol
        tol <- if(method=="Wang" || method=="Wang.Par") 2.2204e-16 # MATLAB default
               else .Machine$double.eps^0.25 # uniroot() default
    worst <- switch(method,
           "Wang" = {

               ## Check qF()
               qF <- NULL # make CRAN check happy
               if(!hasArg(qF))
                   stop("Method 'Wang' requires the quantile function qF of F")
               ## Check 'interval'
               if(is.null(interval)) interval <- c(0, (1-alpha)/d)
               else {
                   if(interval[1] < 0) stop("interval[1] needs to be >= 0")
                   if(interval[1] > (1-alpha)/d) stop("interval[2] needs to be <= (1-alpha)/d")
                   if(interval[1] >= interval[2]) stop("interval[1] needs to be smaller than interval[2]")
               }

               ## Compute (and adjust) function values at endpoints
               h.low <- Wang_h(interval[1], alpha=alpha, d=d, ...)
               if(is.na(h.low))
                   stop("Objective function at interval[1] is NA or NaN. Provide a larger interval[1].")
               h.up <- -h.low # avoid that uniroot() fails due to 0 at upper interval endpoint

               ## Root-finding on 'interval'
               c. <- uniroot(function(c) Wang_h(c, alpha=alpha, d=d, ...),
                             interval=interval, f.lower=h.low, f.upper=h.up, tol=tol)$root
               d * Wang_h_aux(c., alpha=alpha, d=d, ...)

           },
           "Wang.Par" = {

               ## Critical here (see vignette): smaller tolerance and extending the initial interval

               ## Check 'theta'
               theta <- NULL # make CRAN check happy
               if(!hasArg(theta))
                   stop("Method 'Wang.Par' requires the parameter theta")
               th <- list(...)$theta
               stopifnot(length(th) == 1, th > 0) # check theta here

               ## Compute uniroot() initial interval
               if(is.null(interval)) {
                   low <- if(th > 1) {
                       r <- (1-alpha)/((d/(th-1)+1)^th + d-1)
                       r/2 # adjustment to guarantee numerically that h is of opposite sign (required for very large theta)
                   } else if(th == 1) {
                       e <- exp(1)
                       (1-alpha)/((d+1)^(e/(e-1))+d-1)
                   } else {
                       r <- (1-th)*(1-alpha)/d
                       r/2 # adjustment to guarantee numerically that h is of opposite sign (required for very small theta)
                   }
                   up <- if(th == 1) {
                       (1-alpha)/(3*d/2-1)
                   } else {
                       (1-alpha)*(d-1+th)/((d-1)*(2*th+d))
                   }
                   interval <- c(low, up)
               } else {
                   if(interval[1] < 0) stop("interval[1] needs to be >= 0")
                   if(interval[1] > (1-alpha)/d) stop("interval[2] needs to be <= (1-alpha)/d")
                   if(interval[1] >= interval[2]) stop("interval[1] needs to be smaller than interval[2]")
               }

               ## Root-finding on 'interval'
               h <- function(c) Wang_h(c, alpha=alpha, d=d, method="Wang.Par", ...)
               c <- uniroot(h, interval=interval, tol=tol)$root
               d * Wang_h_aux(c, alpha=alpha, d=d, method="Wang.Par", theta=th)

           },
           "dual" = {

               pF <- NULL # make CRAN check happy
               if(!hasArg(pF))
                   stop("Method 'dual' requires the distribution function pF")
               if(!hasArg(interval))
                   stop("Method 'dual' requires an initial interval c(s_l, s_u) to be given")
               uniroot(function(s) dual_bound(s, d=d, tol=tol, ...)-(1-alpha),
                       interval=interval, tol=tol)$root # s interval
               ## Note: We can't pass arguments to the inner root-finding

           },
           stop("Wrong method"))

           ## Return
           c(best, worst)
}


### 3) Worst VaR in the inhomogeneous case #####################################

##' @title Determine the indices which order any increasing (!) vector y
##'        oppositely to x
##' @param x A vector
##' @return order(order(x, decreasing=TRUE)) (= N+1-rank(x))
##' @author Marius Hofert
##' @note For convergence of rearrange() it is crucial to have a stable sorting
##'       procedure underlying (as then no swaps on ties back and forth until
##'       eternity take place which decreases the probability of non-convergence).
##'       The various methods like qsort() in C or rsort_with_index() are *not*
##'       stable. In the way we need it here, rank(, ties.method="last") would
##'       be as well, but internally uses order() and thus is not faster.
##'       However, we can make order() faster by using R_orderVector1() instead
##'       of R_orderVector() (available from 3.2.3 onwards)
##'       => For d=1000 and N=16384, this only brought an improvement of 1.3%, though
##'       If that turns out to be faster (in the future), use:
##'       indices_opp_ordered_to <- if(getRversion() >= "3.2.3")
##'       {
##'          function(x) .Call(C_indices_opp_ordered_to, x)
##'       } else {
##'          function(x) order(order(x, decreasing=TRUE))
##'       }
indices_opp_ordered_to <- function(x) order(order(x, decreasing=TRUE))

##' @title Compute the number of columns oppositely ordered to the sum of all others
##' @param x (N, d)-matrix
##' @return Number of columns oppositely ordered to the sum of all others
##' @author Marius Hofert
##' @note Same time-saving tricks as behind rearrange(), RA() and ARA() (work
##'       with list of columns of x)
num_of_opp_ordered_cols <- function(x) {
    x.rs <- .rowSums(x, nrow(x), ncol(x)) # faster than rowSums()
    x.lst <- .Call(C_col_split, x) # to avoid indexing the jth column, we work with a list!
    x.lst.sorted <- lapply(x.lst, sort.int) # sorting is only necessary once!
    sum(vapply(seq_len(ncol(x)),
               function(j) {
                   xj <- x.lst[[j]]
                   all(x.lst.sorted[[j]][indices_opp_ordered_to(x.rs - xj)]
                       == xj)
               }, NA))
}

##' @title Basic rearrangement function for (A)RA
##' @param X (N, d)-matrix \underline{X}^\alpha or \overline{X}^\alpha
##' @param tol Tolerance to determine (the individual) convergence;
##'        if NULL, column rearrangements are done until the matrix doesn't
##'        change anymore d consecutive times
##' @param tol.type Character string indicating the tolerance function used
##'        ("relative" or "absolute")
##' @param max.ra Maximal number of column rearrangements
##' @param method Character indicating which VaR is approximated (worst/best)
##'        determines optimizing function (min for worst VaR; max
##'        for best VaR)
##' @param sample A logical indicating whether each column of the working
##'        matrix is randomly permuted before the rearrangements begin
##' @param is.sorted A logical indicating whether X is columnwise sorted in
##'        increasing order
##' @param trace A logical indicating whether the underlying matrix is
##'        printed after each rearrangement step
##' @return List containing the
##'         1) Computed (lower or upper [depending on X]) bound for (worst or
##'            best [depending on method]) VaR
##'         2) (Individual) tolerance reached
##'         3) Logical indicating whether the algorithm has converged
##'         4) Vector of minimal [for worst VaR] or maximal [for best VaR]
##'            row sums after each considered column rearrangement
##'         5) The (optimally) rearranged (N, d)-matrix
##' @author Marius Hofert and Kurt Hornik
##' @note - We use "<= tol" to determine convergence instead of "< tol" as
##'         this then also nicely works with "= 0" (if tol=0) which stops in
##'         case the matrices are identical (no change at all).
##'       - We conduct checks of convergence after rearranging each column after the
##'         dth (not only after rearranging all d columns)
##'       - The columns of X have to be given in increasing order if is.sorted=TRUE
##'       - No checking here due to speed! Note that max.ra must be > ncol(X)
rearrange <- function(X, tol=0, tol.type=c("relative", "absolute"),
                      max.ra=Inf, method=c("worst", "best"),
                      sample=TRUE, is.sorted=FALSE, trace=FALSE)
{
    ## Setup
    N <- nrow(X)
    d <- ncol(X)
    tol.type <- match.arg(tol.type)
    method <- match.arg(method)

    ## Define helper functions
    optim.fun <- if(method=="worst") min else max
    tol.fun <- if(tol.type=="absolute") {
        function(x, y) abs(x-y)
    } else {
        function(x, y) abs((x-y)/y)
    }

    ## Tracing
    if(trace) {
        B <- X
        colnames(B) <- rep("", d)
        print(B)
    }

    ## Keep the sorted X
    X.lst.sorted <- if(is.sorted) {
        .Call(C_col_split, X)
    } else {
        .Call(C_col_split, apply(X, 2, sort)) # need to sort first
    }

    ## Sample the columns (if chosen), compute the initial row sum
    ## and the corresponding min/max row sum
    if(sample) {
        X.lst <- lapply(X.lst.sorted, sample) # list of (resampled) columns of X
        X.rs <- .rowSums(do.call(cbind, X.lst), N, d) # row sums of X
    } else {
        X.lst <- X.lst.sorted # list of columns of X
        X.rs <- .rowSums(X, m=N, n=d) # initial row sum
    }

    ## Go through the columns and rearrange one at a time
    iter <- 0 # current iteration number
    j <- 0 # current column number
    num.cols.no.change <- 0 # number of consecutively rearranged columns with no change
    m.row.sums <- c() # vector of minimal/maximal row sums after each rearranged column
    is.null.tol <- is.null(tol)
    while (TRUE) {

        ## Update the running indices
        iter <- iter+1 # current iteration number (in IN)
        j <- if(j >= d) 1 else j+1 # current column

        ## Update the working 'matrix'
        Y.lst <- X.lst # define 'matrix' Y (former 'matrix' X) to work with
        Y.rs <- X.rs # row sum of Y (= row sum of X)

        ## Oppositely order the jth column to the sum of all others
        yj <- Y.lst[[j]] # pick out jth column
        rs <- Y.rs - yj # sum over all other columns (but the jth)
        ## Note: The elements of X.lst.sorted are sorted in increasing order
        ##       which is required for oppositely reordering them
        yj. <- X.lst.sorted[[j]][indices_opp_ordered_to(rs)] # oppositely reorder Y_j

        ## Update the working 'matrix'
        Y.lst[[j]] <- yj. # update with rearranged jth column
        Y.rs <- rs + yj. # update row sum of Y

        ## Tracing
        if(trace) {
            B <- do.call(cbind, Y.lst)
            colnames(B) <- rep("", d)
            no.change <- identical(yj, yj.)
            colnames(B)[j] <- if(no.change) "=" else "|"
            B <- cbind(B, rs, sum=.rowSums(B, m=N, n=d))
            colnames(B)[d+1] <- paste0("-",j)
            print(B)
        }

        ## Update the vector of computed minimal/maximal row sums
        m.rs.cur.col <- optim.fun(Y.rs) # compute new minimal/maximal row sum
        m.row.sums <- c(m.row.sums, m.rs.cur.col) # append it

        ## Check convergence
        ## Idea: After a column has been rearranged, compute the tol (and thus
        ##       determine convergence) between the minimal/maximal row sum
        ##       after that rearrangement and from d steps before when that
        ##       column was rearranged the last time. The earliest we check for
        ##       convergence is when iter > d.
        ## Note: - This is a bit more elegant than the original RA which checked only
        ##         on j=d, not after rearranging *each* column.
        ##       - Checking only *two* consecutive columns led to a bad behavior for ARA()
        ##         in some cases (e.g., real OpRisk data): Both the individual and the joint
        ##         relative tolerances were satisfied but far off (with reltol[1]=0.001).
        ##         Of course one could check d consecutive columns for *all* of them to
        ##         fulfill the 'convergence' criterion, but then what's the reached
        ##         tolerance tol if more than two columns are involved? Maybe the maximum
        ##         tolerance computed over all previous d many rearranged columns?
        ##         There's probably no gain in doing that.
        if(is.null.tol) { # tol = NULL
            num.cols.no.change <- if(identical(yj, yj.)) num.cols.no.change + 1 else 0
            if(num.cols.no.change == d) { # => matrix has not changed in d consecutive col rearrangements
                tol. <- 0 # as there was no change
                tol.reached <- TRUE # as we reached 'no change' in d consecutive steps (we don't care whether max.ra has been reached)
                break
            } else { # check whether we have to stop due to max.ra
                if(iter == max.ra) { # need max.ra > d (as otherwise (*) is wrong)
                    ## Note: iter = number of columns we have already rearranged
                    m.rs.d.col.ago <- m.row.sums[iter-d] # (*)
                    tol. <- tol.fun(m.rs.cur.col, m.rs.d.col.ago) # compute the attained tolerance (in comparison to the last time the jth column was rearranged)
                    tol.reached <- FALSE # as num.cols.no.change < d
                    break
                }
            }
        } else { # tol >= 0
            if(iter > d) {
                m.rs.d.col.ago <- m.row.sums[iter-d]
                tol. <- tol.fun(m.rs.cur.col, m.rs.d.col.ago) # compute the attained tolerance (in comparison to the last time the jth column was rearranged)
                tol.reached <- tol. <= tol
                if(iter == max.ra || tol.reached) break # also here we need max.ra > d; see (*)
            }
        }

        ## Updates for the next column rearrangement
        X.lst <- Y.lst # update the working 'matrix'
        X.rs <- Y.rs # update the row sums

    }

    ## Return
    list(bound=m.rs.cur.col, # computed bound (\underline{s}_N or \overline{s}_N)
         tol=tol., # tolerance for the computed bound
         converged=tol.reached, # indicating whether converged
         m.row.sums=m.row.sums, # the computed minimal/maximal row sums after each column rearrangement
         X.rearranged=do.call(cbind, Y.lst)) # the rearranged matrix X
}

##' @title Computing lower/upper bounds for the worst VaR with the RA
##' @param alpha Confidence level
##' @param qF d-list of marginal quantile functions
##' @param N Number of discretization points
##' @param abstol Absolute convergence tolerance (to determine convergence)
##' @param max.ra Maximal number of column rearrangements
##' @param method Character indicating which VaR is approximated (worst/best)
##' @param sample Logical indicating whether each column of the two working
##'        matrices is randomly permuted before the rearrangements begin
##' @return List containing the
##'         1) Computed lower and upper bound for (worst or best) VaR
##'         2) The relative rearrangement gap
##'            "|(upper bound - lower bound) / upper bound|"
##'         3) Individual absolute tolerances reached (for each bound)
##'         4) 2-vector of logicals indicating whether the individual bounds reached
##'            the desired tolerances (=> convergence)
##'         5) Number of columns considered for rearrangement
##'         6) Vectors of minimal [for worst VaR] or maximal [for best VaR] row sums
##'            after each considered column rearrangement
##'         7) List of (N, d) input matrices X (for each bound)
##'         8) List of rearranged Xs (for each bound)
##' @author Marius Hofert
##' @note Notation is from p. 2757 in Embrechts, Puccetti, Rueschendorf (2013);
##'       variables are named according to the 'worst' VaR case.
RA <- function(alpha, qF, N, abstol=0, max.ra=Inf,
               method=c("worst", "best"), sample=TRUE)
{
    ## Checks and Step 1 (get N, abstol)
    stopifnot(0 < alpha, alpha < 1, is.null(abstol) || abstol >= 0,
              length(N) >= 1, N >= 2, is.logical(sample),
              is.list(qF), sapply(qF, is.function), (d <- length(qF)) >= 2, max.ra > d)
    method <- match.arg(method)

    ## Compute lower bound

    ## Step 2 (build \underline{X}^\alpha)
    p <- if(method=="worst") alpha + (1-alpha)*(0:(N-1))/N else alpha*(0:(N-1))/N # N-vector of prob. in *increasing* order
    X.low <- sapply(qF, function(qF) qF(p))
    ## adjust those that are -Inf (for method="best")
    ## use alpha*((0+1)/2 / N) = alpha/(2N) instead of 0 quantile
    if(method == "best")
        X.low[1,] <- sapply(1:d, function(j)
            if(is.infinite(X.low[1,j])) qF[[j]](alpha/(2*N)) else X.low[1,j])

    ## Steps 3--7 (determine \underline{X}^*)
    ## randomly permute each column of \underline{X}^\alpha and
    ## repeat oppositely ordering \underline{X}^\alpha until there is only an
    ## abstol change in the min (method="worst") or max (method="best") row sum
    ## or until we reached max.ra number of column rearrangements
    res.low <- rearrange(X.low, tol=abstol, tol.type="absolute",
                         max.ra=max.ra, method=method,
                         sample=sample, is.sorted=TRUE)

    ## Compute upper bound

    ## Step 2 (build \overline{X}^\alpha)
    p <- if(method=="worst") alpha + (1-alpha)*(1:N)/N else alpha*(1:N)/N # N-vector of prob. in *increasing* order
    X.up <- sapply(qF, function(qF) qF(p))
    ## adjust those that are Inf (for method="worst")
    ## use alpha+(1-alpha)*(N-1+N)/(2*N) = alpha+(1-alpha)*(1-1/(2*N)) instead of 1 quantile
    if(method == "worst")
        X.up[N,] <- sapply(1:d, function(j)
            if(is.infinite(X.up[N,j])) qF[[j]](alpha+(1-alpha)*(1-1/(2*N))) else X.up[N,j])

    ## Step 3--7 (determine \overline{X}^*)
    ## randomly permute each column of \overline{X}^\alpha and
    ## repeat oppositely ordering \overline{X}^\alpha until there is only an
    ## abstol change in the min (method="worst") or max (method="best") row sum
    ## or until we reached max.ra number of column rearrangements
    res.up <- rearrange(X.up, tol=abstol, tol.type="absolute",
                        max.ra=max.ra, method=method,
                        sample=sample, is.sorted=TRUE)

    ## Return
    optim.fun <- if(method=="worst") min else max
    list(bounds=c(low=res.low$bound, up=res.up$bound), # (\underline{s}_N, \overline{s}_N)
         rel.ra.gap=abs((res.up$bound-res.low$bound)/res.up$bound), # relative RA gap
         ind.abs.tol=c(low=res.low$tol, up=res.up$tol), # individual absolute tolerances
         converged=c(low=res.low$converged, up=res.up$converged), # converged?
         num.ra=c(low=length(res.low$m.row.sums), up=length(res.up$m.row.sums)), # number of considered column rearrangements (low, up)
         m.row.sums=list(low=res.low$m.row.sum, up=res.up$m.row.sums), # optimal row sums (low, up)
         X=list(low=X.low, up=X.up), # input matrices X (low, up)
         X.rearranged=list(low=res.low$X.rearranged, up=res.up$X.rearranged)) # rearranged Xs (low, up)
}

##' @title Computing lower/upper bounds for the worst VaR with the ARA
##' @param alpha Confidence level
##' @param qF d-list of marginal quantile functions
##' @param N.exp Vector of exponents of 2 used as discretization points
##' @param reltol 2-vector of relative convergence tolerances
##'        for determining the individual relative tolerance (i.e., the relative
##'        tolerance in the minimal/maximal row sum for each of the bounds) and
##'        the joint relative tolerance (i.e., the relative
##'        tolerance between the computed lower and upper bounds).
##' @param max.ra Maximal number of column rearrangements per N
##' @param method Character indicating which VaR is approximated (worst/best)
##' @param sample Logical indicating whether each column of the two working
##'        matrices is randomly permuted before the rearrangements begin
##' @return List containing the
##'          1) Computed lower and upper bound for (worst or best) VaR
##'          2) The relative rearrangement gap
##'             "|(upper bound - lower bound) / upper bound|"
##'          3) Relative tolerances reached (individually for each bound and jointly
##'             between the bounds)
##'          4) 3-vector of logicals indicating whether the individual bounds and
##'             the two bounds jointly reached the desired tolerances (=> convergence)
##'          5) The number of discretization points used
##'          6) Number of columns considered for rearrangement
##'          7) Vectors of minimal [for worst VaR] or maximal [for best VaR] row sums
##'             after each considered column rearrangement
##'          8) List of (N, d) input matrices X (for each bound)
##'          9) List of rearranged Xs (for each bound)
##' @author Marius Hofert
ARA <- function(alpha, qF, N.exp=seq(8, 19, by=1), reltol=c(0, 0.01),
                max.ra=10*length(qF), method=c("worst", "best"), sample=TRUE)
{
    ## Checks and Step 1 (get N, reltol)
    lreltol <- length(reltol)
    stopifnot(0 < alpha, alpha < 1, lreltol==1 || lreltol==2, reltol >= 0,
              length(N.exp) >= 1, N.exp >= 1, is.logical(sample),
              is.list(qF), sapply(qF, is.function), (d <- length(qF)) >= 2, max.ra > d)
    method <- match.arg(method)

    ## Determine tolerances
    itol <- if(lreltol == 2) reltol[1] else NULL # individual tolerance
    jtol <- if(lreltol == 2) reltol[2] else reltol[1] # joint tolerance

    ## Loop over N
    for(N in 2^N.exp) {

        ## Compute lower bound

        ## Step 2 (build \underline{X}^\alpha)
        p <- if(method=="worst") alpha + (1-alpha)*(0:(N-1))/N else alpha*(0:(N-1))/N # N-vector of prob. in *increasing* order
        X.low <- sapply(qF, function(qF) qF(p))
        ## adjust those that are -Inf (for method="best")
        ## use alpha*((0+1)/2 / N) = alpha/(2N) instead of 0 quantile
        if(method == "best")
            X.low[1,] <- sapply(1:d, function(j)
                if(is.infinite(X.low[1,j])) qF[[j]](alpha/(2*N)) else X.low[1,j])

        ## Steps 3--7 (determine \underline{X}^*)
        ## randomly permute each column of \underline{X}^\alpha and
        ## repeat oppositely ordering \underline{X}^\alpha until there is only an
        ## itol change in the min (method="worst") or max (method="best") row sum
        ## or until we reached max.ra number of column rearrangements
        res.low <- rearrange(X.low, tol=itol, tol.type="relative",
                             max.ra=max.ra, method=method,
                             sample=sample, is.sorted=TRUE)

        ## Compute upper bound

        ## Step 2 (build \overline{X}^\alpha)
        p <- if(method=="worst") alpha + (1-alpha)*(1:N)/N else alpha*(1:N)/N # N-vector of prob. in *increasing* order
        X.up <- sapply(qF, function(qF) qF(p))
        ## adjust those that are Inf (for method="worst")
        ## use alpha+(1-alpha)*(N-1+N)/(2*N) = alpha+(1-alpha)*(1-1/(2*N)) instead of 1 quantile
        if(method == "worst")
            X.up[N,] <- sapply(1:d, function(j)
                if(is.infinite(X.up[N,j])) qF[[j]](alpha+(1-alpha)*(1-1/(2*N))) else X.up[N,j])

        ## Step 3--7 (determine \overline{X}^*)
        ## randomly permute each column of \overline{X}^\alpha and
        ## repeat oppositely ordering \overline{X}^\alpha until there is only an
        ## itol change in the min (method="worst") or max (method="best") row sum
        ## or until we reached max.ra number of column rearrangements
        res.up <- rearrange(X.up, tol=itol, tol.type="relative",
                            max.ra=max.ra, method=method,
                            sample=sample, is.sorted=TRUE)

        ## Determine (individual and joint) convergence
        joint.tol <- abs((res.low$bound-res.up$bound)/res.up$bound)
        joint.tol.reached <- joint.tol <= jtol
        if(res.low$converged && res.up$converged && joint.tol.reached) break

    }

    ## Return
    optim.fun <- if(method=="worst") min else max
    list(bounds=c(low=res.low$bound, up=res.up$bound), # (\underline{s}_N, \overline{s}_N)
         rel.ra.gap=abs((res.up$bound-res.low$bound)/res.up$bound), # relative RA gap
         rel.tol=c(low=res.low$tol, up=res.up$tol, joint=joint.tol), # individual and joint relative tolerances
         converged=c(low=res.low$converged, up=res.up$converged, joint=joint.tol.reached), # converged?
         N.used=N, # number of discretization points used
         num.ra=c(low=length(res.low$m.row.sums), up=length(res.up$m.row.sums)), # number of considered column rearrangements (low, up)
         m.row.sums=list(low=res.low$m.row.sums,
                         up=res.up$m.row.sums), # optimal row sums (low, up) for the N used
         X=list(low=X.low, up=X.up), # input matrices X (low, up)
         X.rearranged=list(low=res.low$X.rearranged, up=res.up$X.rearranged)) # rearranged Xs (low, up)
}
