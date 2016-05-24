###
##
##  Qrot(), originally in  /u/maechler/R/MM/STATISTICS/robust/MV-WSt-ex.R
##  ----
##
## Idea: (alternative:  that rotates  e_i {i-th unit vector} on to
##	e_jk where (e_jk)_i = 1_[i \in {j,k}]


### Construct the Rotation matrix, that rotates (1,0,....0) onto ### (1,1,1,...1)/sqrt(p), or more generally to  u / ||u||  (u := unit.image)
Qrot <- function(p, transpose = FALSE, unit.image = rep(1,p))
{
    ## Purpose: Construct a  p x p  rotation matrix which rotates
    ##	(1,0,....0) on (1,1,,...1)/sqrt(p), i.e. the x1-axis onto the diagonal
    ## ----------------------------------------------------------------------
    ## Arguments: p : dimension
    ##          unit.image: vector onto which (1, 0,...,0) is rotated
    ## ----------------------------------------------------------------------
    ## Author: Martin Maechler, Date: 3 Dec 2005; 28 Nov 2008
    D <- diag(p)
    ## for some reason, the " - " below also keeps the sign correct:
    Q <- qr.qy(qr(cbind(unit.image, D[,-1])), - D)
    if(transpose) Q else t(Q)
}


rbwheel <- function(n,		# observations
		    p,		# variables
		    frac = 1/p, # proportion of outliers
		    sig1 = .05, # thickness of the 'wheel' = sigma(good[,1])
		    sig2 = 1/10,# diameter of the 'axle' (compared to 1)
                    rGood = rnorm,# generator for "good" observations
                    rOut = function(n)sqrt(rchisq(n,p-1))*sign(runif(n,-1,1)),
		    U1 = rep(1, p), ## Vector to which (1,0,...,0) is rotated
                    scaleAfter = TRUE, scaleBefore = FALSE, spherize = FALSE,
		    fullResult = FALSE)
{
    ## Purpose: simulate data according to the 'wheel' distribution
    ## ----------------------------------------------------------------------
    ## Arguments:
    ## ----------------------------------------------------------------------
    ## Author: Werner Stahel, Martin Maechler, Date: 28 Nov 2008, 15:27
    stopifnot(is.numeric(frac), 0 <= frac, frac < 1,
	      n >= 1, p >= 2)
    ## compatibility-warning -- at most once per session :
    if(missing(scaleAfter) &&
       (is.null(w <- getOption("rbwheel.warn.scaleA")) || isTRUE(w))) {
	if(is.null(w)) options( rbwheel.warn.scaleA = FALSE)
	warning("Note: rbwheel() now uses scaleAfter = TRUE  by default")
    }

    ## a simplified version of scale.default :
    scale.simply <- function(x) {
        x <- sweep(x, 2, colMeans(x, na.rm = TRUE), check.margin = FALSE)
        sdev <- apply(x, 2,
                       function(v) {
                           v <- v[!is.na(v)]
                           sqrt(sum(v^2)/max(1, length(v) - 1L))
                       })
        sweep(x, 2, sdev, "/", check.margin = FALSE)
    }

    n1 <- pmax(0, pmin(n, round((1-frac)*n)))
    n2 <- n-n1 ## ~= frac * n
    i <- if(n1 < n) (n1+1):n else integer(0) # index of n2 "outliers"

    d0 <- matrix(0, n,p)
    d0[ , -1] <- rGood(n*(p-1))
    d0[i, -1] <- sig2 * d0[i,-1]
    d0[-i, 1] <- sig1 * rGood(n1)
    d0[ i, 1] <- rOut(n2)

    d1 <- {
	if(spherize) { # use Chol(), such that X_1 remains unchanged:
            d. <- scale.simply(d0)
            t(backsolve(chol(cov(d.)), t(d.)))
        }
	else if(scaleBefore)
	    scale.simply(d0)
	else d0
    }

    maybeScale <- function(x) if(scaleAfter) scale.simply(x) else x
    if(fullResult) { ## for didactical reasons mainly, see example
	A <- Qrot(p, unit.image = U1)
	list(X = maybeScale(d1 %*% A), X0 = d0, A = A, n1 = n1, n2 = n2,
             scale = c(before=scaleBefore, after=scaleAfter, spherize=spherize))
    }
    else ## by default -- 'n1' as attribute :
	structure(maybeScale(d1 %*% Qrot(p, unit.image = U1)), n1 = n1)
}

