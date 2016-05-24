diagDA <- function(ls, cll, ts, pool= TRUE)
{
    ## Purpose: Diagonal (Linear or Quadratic) Discriminant Analysis
    ## ----------------------------------------------------------------------
    ## Arguments: --> ?diagDA  (i.e. ../man/diagDA.Rd )
    ## ----------------------------------------------------------------------
    ## Authors: Sandrine Dudoit, sandrine@stat.berkeley.edu
    ##	        Jane Fridlyand, janef@stat.berkeley.edu
    ## as function  stat.diag.da() in package "sma"
    ##
    ## Modification (API and speed): Martin Maechler, Date: 19 Nov 2003, 15:34

### ---------------------- Fit Model ------------------------------
    ls <- data.matrix(ls)
    n <- nrow(ls)
    p <- ncol(ls)

    cl0 <- as.integer(min(cll, na.rm=TRUE) - 1)
    cll <- as.integer(cll) - cl0 ## cll now in 1:K
    inaC <- is.na(cll)
    clL <- cll[!inaC]
    K <- max(clL)
    if(K != length(unique(clL)))
        stop(sQuote("cll")," did not contain *consecutive* integers")

    nk <- integer(K)
    m <- v <- matrix(0,p,K)

    colVars <- function(x, means = colMeans(x, na.rm = na.rm), na.rm=FALSE) {
        x <- sweep(x, 2, means)
        colSums(x*x, na.rm = na.rm) / (nrow(x) - 1)
    }
    sum.na <- function(x) sum(x, na.rm=TRUE)

    ## Class means and variances
    for(k in 1:K) {
        which <- (cll == k)
        nk[k] <- sum.na(which)
        lsk <- ls[which, , drop = FALSE]
        m[,k] <- colMeans(lsk, na.rm = TRUE)
        if(nk[k] > 1)
            v[,k] <- colVars (lsk, na.rm = TRUE, means = m[,k]) ## else 0
    }

### ---------------------- Predict from Model -----------------------------

    ts <- data.matrix(ts)
    if(p != ncol(ts))
        stop("test set matrix must have same columns as learning one")
    ## any NA's in test set currently must give NA predictions
    ts <- na.exclude(ts)
    nt <- nrow(ts)
    disc <- matrix(0, nt,K)

    if(pool) { ## LDA
        ## Pooled estimates of variances
        vp <- rowSums(rep(nk - 1, each=p) * v) / (n - K)
        ## == apply(v, 1, function(z) sum.na((nk-1)*z))/(n-K)
        if(any(i0 <- vp == 0)) vp[i0] <- 1e-7 * min(vp[!i0])

        ivp <- rep(1/vp, each = nt) # to use in loop

        for(k in 1:K) {
            y <- ts - rep(m[,k], each=nt)
            disc[,k] <- rowSums(y*y * ivp)
            ## == apply(ts, 1, function(z) sum.na((z-m[,k])^2/vp))
        }
    }
    else { ## QDA
if(FALSE) { ## not yet quite : fails ../tests/dDA.R  -- FIXME
        for(k in 1:K) {
            ts <- ts - rep(m[,k], each=nt)
            disc[,k] <- rowSums((ts*ts) / rep(v[,k], each=nt)) + sum(log(v[,k]))
        }
} else {
        for(k in 1:K) {
            disc[,k] <-
                apply(ts,1, function(z) sum((z-m[,k])^2/v[,k])) +
                    sum.na(log(v[,k]))
        }
}
    }

    ## predictions

    pred <- cl0 + apply(disc, 1, which.min)
    if(inherits(attr(ts,"na.action"), "exclude")) # had missings in `ts'
        pred <- napredict(omit = attr(ts,"na.action"), pred)
    pred
}

## Cleaner: One function to estimate; one to predict :
## ------- (my tests give a time-penalty 5% for doing things two steps)

dDA <- function(x, cll, pool= TRUE)
{
    ## Purpose: Diagonal (Linear or Quadratic) Discriminant Analysis

    x <- data.matrix(x)
    n <- nrow(x)
    p <- ncol(x)

    cl0 <- as.integer(min(cll, na.rm=TRUE) - 1)
    cll <- as.integer(cll) - cl0 ## cll now in 1:K
    inaC <- is.na(cll)
    clL <- cll[!inaC]
    K <- max(clL)
    if(K != length(unique(clL)))
        stop(sQuote("cll")," did not contain *consecutive* integers")

    nk <- integer(K)
    m <- v <- matrix(0,p,K)

    colVars <- function(x, means = colMeans(x, na.rm = na.rm), na.rm=FALSE) {
        x <- sweep(x, 2, means)
        colSums(x*x, na.rm = na.rm) / (nrow(x) - 1)
    }
    sum.na <- function(x) sum(x, na.rm=TRUE)

    ## Class means and variances
    for(k in 1:K) {
        which <- (cll == k)
        nk[k] <- sum.na(which)
        lsk <- x[which, , drop = FALSE]
        m[,k] <- colMeans(lsk, na.rm = TRUE)
        if(nk[k] > 1)
            v[,k] <- colVars (lsk, na.rm = TRUE, means = m[,k]) ## else 0
    }
    structure(list(call = match.call(), cl0 = cl0, n=n, p=p, K=K,
                   means=m, vars=v, nk=nk, pool=pool),
              class = "dDA")
}

print.dDA <- function(x, ...)
{
    cat(if(x$pool)"Linear (pooled var)" else "Quadratic (no pooling)",
        "Diagonal Discriminant Analysis,\n ", deparse(x$call),"\n")
    with(x,
         cat("  (n= ",n,") x (p= ",p,") data in K=",K," classes of [",
             paste(nk, collapse=", "),"] observations each\n", sep=""))
    cat("\n")
    invisible(x)
}

predict.dDA <- function(object, newdata, pool = object$pool, ...)
{
    newdata <- data.matrix(newdata)
    n <- object$n
    p <- object$p
    K <- object$K
    ## means and vars  are  (p x K) matrices:
    mu <- object$means
    Vr <- object$vars
    if(p != ncol(newdata))
        stop("test set matrix must have same columns as learning one")
    ## any NA's in test set currently must give NA predictions
    newdata <- na.exclude(newdata)
    nt <- nrow(newdata)
    disc <- matrix(0, nt,K)

    if(pool) { ## LDA
        ## Pooled estimates of variances
        vp <- rowSums(Vr * rep(object$nk - 1, each=p)) / (n - K)
        ## == apply(Vr, 1, function(z) sum.na((nk-1)*z))/(n-K)
        if(any(i0 <- vp == 0)) vp[i0] <- 1e-7 * min(vp[!i0])

        ivp <- rep(1/vp, each = nt) # to use in loop

        for(k in 1:K) {
            y <- newdata - rep(mu[,k], each=nt)
            disc[,k] <- rowSums(y*y * ivp)
            ## == apply(newdata, 1, function(z) sum.na((z-mu[,k])^2/vp))
        }
    }
    else { ## QDA
        sum.na <- function(x) sum(x, na.rm=TRUE)
        ## zero - variances are not acceptable later
	if(any(i0 <- Vr == 0)) {
	    if(all(i0))
		stop("all variances are 0 -- cannot predict")
	    Vr[i0] <- 1e-7 * min(Vr[!i0])
	}

if(FALSE) { ## not yet quite : fails ../tests/dDA.R  -- FIXME
        for(k in 1:K) {
            y <- newdata - rep(mu[,k], each=nt)
            disc[,k] <- rowSums((y*y) / rep(Vr[,k], each=nt)) + sum(log(Vr[,k]))
        }
} else {
        for(k in 1:K) {
            disc[,k] <-
                apply(newdata,1, function(z) sum((z-mu[,k])^2/Vr[,k])) +
                    sum.na(log(Vr[,k]))
        }
}
    }

    ## predictions

    pred <- object$cl0 + apply(disc, 1, which.min)
    if(inherits(attr(newdata,"na.action"), "exclude")) {
        ## had missings in `newdata'
        pred <- napredict(omit = attr(newdata,"na.action"), pred)
    } ##        ^^^^^^^^^ typically stats:::napredict.exclude()
    pred
}

