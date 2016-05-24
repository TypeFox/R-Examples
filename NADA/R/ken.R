#-->> BEGIN cenken code

## Generics

setGeneric("cenken", function(y, ycen, x, xcen) standardGeneric("cenken"))

## Classes 

setClass("cenken", "NADAList")

#setClass("cenken", 
#         representation(slope="numeric", intercept="numeric",
#                        tau="numeric", p="numeric"))


## Methods

setMethod("cenken",
          signature(y="numeric", ycen="logical", x="numeric", xcen="logical"),
          function(y, ycen, x, xcen)
{
    kendallATS(y, ycen, x, xcen) 
})

setMethod("cenken",
          signature(y="numeric", ycen="logical", x="numeric", xcen="missing"),
          function(y, ycen, x, xcen)
{
    xcen = rep(FALSE, length(x))
    kendallATS(y, ycen, x, xcen)
})

setMethod("cenken",
          signature(y="formula", ycen="missing", x="missing", xcen="missing"),
          function(y, ycen, x, xcen)
{
    f = y

    y    = eval.parent(f[[2]][[2]])
    ycen = eval.parent(f[[2]][[3]])

    x    = eval.parent(f[[3]])
    xcen = rep(F, length(y))

    kendallATS(y, ycen, x, xcen)
})

setMethod("show", signature(object="cenken"), function(object)
{
    show(as(object, "NADAList"))
})

setMethod("lines", signature(x="cenken"), function (x, ...)
{
     abline(a=x$intercept, b=x$slope, ...)
})

## Broken for the time being -- use lines
#setMethod("abline", signature(a="cenken"),
##function(a, b, h, v, reg, coef, untf, col, par, lty, ...)
#function (a = NULL, b = NULL, h = NULL, v = NULL, reg = NULL,
#          coef = NULL, untf = FALSE, col = par("col"), lty = par("lty"),
#          lwd = NULL, ...)
#{
#    #abline(a=a$intercept, b=a$slope, ...)
#    abline(a=a$intercept, b=a$slope, h=h, v=v, reg=reg, coef=coef, 
#           untf=untf, col=col, lty=lty, lwd=lwd, ...)
#})

## kendallATS function -- the heart of the cenken routines
#  Original S code written by D. Lorenz for S-Plus.
#  Port to R and R-native ktau function by L. Lee and D. Helsel.
#
# Kendall's tau used as an estimate of the relation between y and x
#    with the ATS slope estimator (y and x left-censored).
#
kendallATS =
function(y, ycen, x, xcen, tol=1e-7, iter=1e+3) 
{
    # Jitter y so that detects and nondetects don't tie.
    # This is only needed in ktau_s but is done here
    # because the original Fortran code (the "standard") does this.
    y = y - (min(y)/1000) * ycen 

    iter_s = 
    function(lb, ub, x, xcen, y, ycen, tol=tol, iter=iter) 
    {
        step_s = 
        function(lb, ub, x, xcen, y, ycen)
        {
          b = c(lb, (lb+ub)/2, ub)
          res = sapply(b, function(i) y - i * x)
          s = apply(res, 2, ktau_s, x=x, xcen=xcen, ycen=ycen)
          list(b=b, s=s)
        }

        bs = step_s(lb, ub, x, xcen, y, ycen)

        for (i in 1:iter)
          {
            b = bs$b
            s = bs$s

            if ((s[1] * s[2]) <= 0) { s[3] = s[2]; b[3] = b[2] }
            else                    { s[1] = s[2]; b[1] = b[2] }

            if ( (s[2] == 0) || (abs(b[3] - b[1]) <= tol) ) break()

            bs = step_s(b[1], b[3], x, xcen, y, ycen)
          }
        return(bs)
    }

    k = ktau_b(x, xcen, y, ycen)
    bs = iter_s(k[1], k[2], x, xcen, y, ycen, tol=tol, iter=iter)

    ubs = iter_s(bs$b[2], bs$b[3], x, xcen, y, ycen, tol=tol, iter=iter)
    lbs = iter_s(bs$b[1], bs$b[2], x, xcen, y, ycen, tol=tol, iter=iter)

    slope = 0.5 * (ubs$b[2] + lbs$b[2])
    
    int = turnbull(y, ycen, x, xcen, slope)

    p_tau = ktau_p(x, xcen, y, ycen)

    ret = list(slope=slope, intercept=int, tau=p_tau$tau, p=p_tau$p)
    class(ret) = "cenken"

    return(ret)
}

turnbull =
function(y, ycen, x, xcen, slope, tol=.Machine$double.eps)
{
    # TODO: make sure no negative residuals by bumping up resids
    resid = y[!xcen] - slope * x[!xcen]

    A=cbind(L = resid - (y[!xcen] * ycen[!xcen]), R = resid)

    em = EM(A, tol=tol)
    surv = rev(cumsum(rev(em$pf)))[-1]
    int = em$intmap[2,][min(which(surv <= 0.5))]

    return(int)
}

ktau_b =
function (x, xcen, y, ycen) 
{
    cx = !xcen
    cy = !ycen

    slopes = unlist(lapply(seq(along = x), function(i, y, x) ((y[i] - 
        y[1:i])/(x[i] - x[1:i])), y * cy, x * cx))

    return(range(slopes[is.finite(slopes)]))
}

ktau_s =
function(x, xcen, y, ycen)
{
    # Jitter y so that detects and nondetects don't tie
    #y = y - (min(y)/1000) * ycen

    ## Find signs of all x<->y slopes (xy2)
    xy1 = sign(outer(y, y, "-"))       
    xy2 = xy1 * sign(outer(x, x, "-"))

    itot = 0.5 * sum( (outer(ycen, -ycen, "-") == 0) * xy2)

    xy3 = outer(ycen, ycen, "-")
    xy4 = outer(y, y, "<=")

    itot = itot + sum(((xy3 * xy4) == 1) * xy2)

    return(itot)
}

#  Original S code written by D. Lorenz for S-Plus.
#  To do: simplify.
ktau_p =
function (x, xcen, y, ycen) 
{
    xx <- x
    cx <- xcen
    yy <- y
    cy <- ycen
    n  <- length(xx)

    delx <- min(diff(sort(unique(xx))))/1000
    dely <- min(diff(sort(unique(yy))))/1000
    dupx <- xx - delx * cx
    diffx <- outer(dupx, dupx, "-")
    diffcx <- outer(cx, cx, "-")
    xplus <- outer(cx, -cx, "-")
    dupy <- yy - dely * cy
    diffy <- outer(dupy, dupy, "-")
    diffcy <- outer(cy, cy, "-")
    yplus <- outer(cy, -cy, "-")
    signyx <- sign(diffy * diffx)
    tt <- (sum(1 - abs(sign(diffx))) - n)/2
    uu <- (sum(1 - abs(sign(diffy))) - n)/2
    cix <- sign(diffcx) * sign(diffx)
    cix <- ifelse(cix <= 0, 0, 1)
    tt <- tt + sum(cix)/2
    signyx <- signyx * (1 - cix)
    ciy <- sign(diffcy) * sign(diffy)
    ciy <- ifelse(ciy <= 0, 0, 1)
    uu <- uu + sum(ciy)/2
    signyx <- signyx * (1 - ciy)
    xplus <- ifelse(xplus <= 1, 0, 1)
    yplus <- ifelse(yplus <= 1, 0, 1)
    diffx <- abs(sign(diffx))
    diffy <- abs(sign(diffy))
    tplus <- xplus * diffx
    uplus <- yplus * diffy
    tt <- tt + sum(tplus)/2
    uu <- uu + sum(uplus)/2

    itot <- sum(signyx * (1 - xplus) * (1 - yplus))
    kenS <- itot/2
    tau <- (itot)/(n * (n - 1))

    J <- n * (n - 1)/2
    taub <- kenS/(sqrt(J - tt) * sqrt(J - uu))
    
    varS <- n * (n - 1) * (2 * n + 5)/18

    intg <- 1:n
    dupx <- xx - delx * cx
    dupy <- yy - dely * cy
    dorder <- order(dupx)
    dxx <- dupx[dorder]
    dcx <- cx[dorder]
    dorder <- order(dupy)
    dyy <- dupy[dorder]
    dcy <- cy[dorder]
    tmpx <- dxx - intg * (1 - dcx) * delx
    tmpy <- dyy - intg * (1 - dcy) * dely
    rxlng <- rle(rank(tmpx))$lengths
    nrxlng <- table(rxlng)
    rxlng <- as.integer(names(nrxlng))
    x1 <- nrxlng * rxlng * (rxlng - 1) * (2 * rxlng + 5)
    x2 <- nrxlng * rxlng * (rxlng - 1) * (rxlng - 2)
    x3 <- nrxlng * rxlng * (rxlng - 1)
    rylng <- rle(rank(tmpy))$lengths
    nrylng <- table(rylng)
    rylng <- as.integer(names(nrylng))
    y1 <- nrylng * rylng * (rylng - 1) * (2 * rylng + 5)
    y2 <- nrylng * rylng * (rylng - 1) * (rylng - 2)
    y3 <- nrylng * rylng * (rylng - 1)
    delc <- (sum(x1) + sum(y1))/18 - sum(x2) * sum(y2)/(9 * n * 
        (n - 1) * (n - 2)) - sum(x3) * sum(y3)/(2 * n * (n - 
        1))
    x4 <- nrxlng * (rxlng - 1)
    y4 <- nrylng * (rylng - 1)
    tmpx <- intg * dcx - 1
    tmpx <- ifelse(tmpx < 0, 0, tmpx)
    nrxlng <- sum(tmpx)
    rxlng <- 2
    x1 <- nrxlng * rxlng * (rxlng - 1) * (2 * rxlng + 5)
    x2 <- nrxlng * rxlng * (rxlng - 1) * (rxlng - 2)
    x3 <- nrxlng * rxlng * (rxlng - 1)
    tmpy <- intg * dcy - 1
    tmpy <- ifelse(tmpy < 0, 0, tmpy)
    nrylng <- sum(tmpy)
    rylng <- 2
    y1 <- nrylng * rylng * (rylng - 1) * (2 * rylng + 5)
    y2 <- nrylng * rylng * (rylng - 1) * (rylng - 2)
    y3 <- nrylng * rylng * (rylng - 1)
    deluc <- (sum(x1) + sum(y1))/18 - sum(x2) * sum(y2)/(9 * 
        n * (n - 1) * (n - 2)) - sum(x3) * sum(y3)/(2 * n * (n - 
        1)) - (sum(x4) + sum(y4))
    dxx <- dxx - intg * dcx * delx
    dyy <- dyy - intg * dcy * dely
    rxlng <- rle(rank(dxx))$lengths
    nrxlng <- table(rxlng)
    rxlng <- as.integer(names(nrxlng))
    x1 <- nrxlng * rxlng * (rxlng - 1) * (2 * rxlng + 5)
    x2 <- nrxlng * rxlng * (rxlng - 1) * (rxlng - 2)
    x3 <- nrxlng * rxlng * (rxlng - 1)
    rylng <- rle(rank(dyy))$lengths
    nrylng <- table(rylng)
    rylng <- as.integer(names(nrylng))
    y1 <- nrylng * rylng * (rylng - 1) * (2 * rylng + 5)
    y2 <- nrylng * rylng * (rylng - 1) * (rylng - 2)
    y3 <- nrylng * rylng * (rylng - 1)
    delu <- (sum(x1) + sum(y1))/18 - sum(x2) * sum(y2)/(9 * n * 
        (n - 1) * (n - 2)) - sum(x3) * sum(y3)/(2 * n * (n - 
        1))

    varS <- varS - delc - deluc - delu

    p.val <- 2 * (1 - pnorm((abs(kenS - sign(kenS)))/sqrt(varS)))

    return(list(tau=tau, p=p.val))
}

#-->> END cenken code
