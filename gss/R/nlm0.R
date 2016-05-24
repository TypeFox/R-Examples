## minimization of univariate function on finite intervals
## using 3-point quadratic fit with golden-section safe-guard
nlm0 <- function(fun,range,prec=1e-7)
{
    ratio <- 2/(sqrt(5)+1)
    ll.x <- min(range)
    uu.x <- max(range)
    if (uu.x-ll.x<prec) {
        sol <- (ll.x+uu.x)/2
        fit <- fun(sol)
        return(list(estimate=sol,minimum=fit,evaluations=1))
    }
    ml.x <- uu.x - ratio*(uu.x-ll.x)
    mu.x <- ll.x + ratio*(uu.x-ll.x)
    ## Initialization
    uu.fit <- fun(uu.x)
    mu.fit <- fun(mu.x)
    ml.fit <- fun(ml.x)
    ll.fit <- fun(ll.x)
    neval <- 4
    ## Iteration
    repeat {
        ## Fit a parabola to the 3 best points and find its minimum
        if (ll.fit<uu.fit) {
            delta.l <- ml.x-ll.x
            sigma.l <- ml.x+ll.x
            d.l <- (ml.fit-ll.fit)/delta.l
            delta.u <- mu.x-ml.x
            d.u <- (mu.fit-ml.fit)/delta.u
        }
        else {
            delta.l <- mu.x-ml.x
            sigma.l <- mu.x+ml.x
            d.l <- (mu.fit-ml.fit)/delta.l
            delta.u <- uu.x-mu.x
            d.u <- (uu.fit-mu.fit)/delta.u
        }
        a <- (d.u-d.l)/(delta.l+delta.u)
        b <- d.l-a*sigma.l
        if (a<=0) nn.x <- max(range)
        else nn.x <- -b/2/a
        ## New bracket
        if (ml.fit<mu.fit) {
            uu.x <- mu.x
            uu.fit <- mu.fit
            mm.x <- ml.x
            mm.fit <- ml.fit
        }
        else {
            ll.x <- ml.x
            ll.fit <- ml.fit
            mm.x <- mu.x
            mm.fit <- mu.fit
        }
        range.l <- mm.x-ll.x
        range.u <- uu.x-mm.x
        delta <- min(abs(nn.x-c(ll.x,mm.x,uu.x)))
        ## Safeguard
        if ((nn.x<ll.x)|(nn.x>uu.x)|(delta<prec)) {
            if (range.u>range.l) nn.x <- uu.x - ratio*range.u
            else nn.x <- ll.x + ratio*range.l
        }
        ## Update middle points
        nn.fit <- fun(nn.x)
        neval <- neval + 1
        if (nn.x<mm.x) {
            ml.x <- nn.x
            ml.fit <- nn.fit
            mu.x <- mm.x
            mu.fit <- mm.fit
        }
        else {
            ml.x <- mm.x
            ml.fit <- mm.fit
            mu.x <- nn.x
            mu.fit <- nn.fit
        }
        ## Return results
        if ((range.l+range.u<.5)&(abs(mm.x-nn.x)<sqrt(prec))) {
            if (nn.fit<mm.fit) {
                solution <- nn.x
                fit <- nn.fit
            }
            else {
                solution <- mm.x
                fit <- mm.fit
            }
            break
        }
    }
    list(estimate=solution,minimum=fit,evaluations=neval)
}
