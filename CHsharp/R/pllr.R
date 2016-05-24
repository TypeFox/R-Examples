pllr <- function(x, y) {
    # penalized local linear regression using a roughness penalty
    M<- 401
    yy <- y[order(x)];  xx <- sort(x); n <- length(x)
    xgrid <- seq(min(xx), max(xx),length=M)
# Choosing lambda and h based on MISE
    CVscore <- function(par, x, y, xgrid) {
        htry <- par[1]
        lambdatry <- par[2]
        n <- length(x)
        indices <- seq(2,n-1,2)
        xs <- x[indices]
        ys <- y[indices]
        ys.lp <- penlocreg(xs, ys, xgrid, degree = 1, h = htry, lambda = lambdatry,
           L = SecondDerivativePenalty)
        ys.lpr <- locpoly(xs,ys.lp$y,degree=1,bandwidth=htry)
        ghat <- function(x) approx(ys.lpr$x, ys.lpr$y, xout=x)$y
        xtesting <- x[-indices]
        ytesting <- y[-indices]
        mean(na.omit((ghat(xtesting)-ytesting)^2))
    }
    h.dpill <- dpill(xx,yy)
    h.min <- max(h.dpill/10, max(diff(sort(xx)))/2)
    h.dpill <- max(h.dpill, h.min)
    out <- optim(c(h.dpill, 1e-7), CVscore, lower=c(h.min, 1e-10), method="L-BFGS-B",
           upper=c(h.dpill*10, 100), x=xx, y=yy, xgrid=xgrid)$par
    h.out <- out[1]
    lambda.out <- out[2]
    yy.lp <- penlocreg(xx, yy, xgrid, degree = 1, h = h.out/2^.2, lambda = lambda.out,
        L=SecondDerivativePenalty)
    lambda.out <- lambda(xx, yy, h.out/2^.2, 1, xgrid, yy.lp$A, yy.lp$B, 3)[3]
    yy.lp <- penlocreg(xx, yy, xgrid, degree = 1, h = h.out/2^.2, lambda = lambda.out,
        L=SecondDerivativePenalty)
    yy.lpr <- locpoly(xx,yy.lp$y,degree=1,bandwidth=h.out/2^.2)
    list(x = yy.lpr$x, y = yy.lpr$y, h = h.out, lambda = lambda.out)
}

