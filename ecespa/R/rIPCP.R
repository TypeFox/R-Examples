rIPCP <-
function (x, lambda = NULL, type = 1, lmax = NULL, win = owin(c(0,
    1), c(0, 1)), ...)
{
    x.ppp <- x$data
    xw <- x.ppp$w
    rho <- x$rho
    sigma <- sqrt(x$sigma2)
    mu <- x.ppp$n/(rho * area.owin(x.ppp$w))
    if (is.null(lambda))
        lambda <- x$lambda
    win <- if (is.im(lambda))
        rescue.rectangle(as.owin(lambda))
    else as.owin(win)
    x.ppp$window <- win
    if (is.null(lmax)) {
        imag <- as.im(lambda, win, ...)
        summ <- summary(imag)
        lmax <- summ$max
    }
    if (is.im(lambda)) {
        probm <- mean(eval.im(lambda/lmax))
        n.preth <- x.ppp$n/probm
        mu.preth <- n.preth * mu/x.ppp$n
        rho.preth <- rho/probm
        if (type == 1)
            X <- rThomas(win = x.ppp$w, kappa = rho, sigma = sigma,
                mu=mu.preth)
        if (type == 2)
            X <- rThomas(win = x.ppp$w, kappa = rho.preth, sigma = sigma,
                mu = mu)
        if (X$n == 0)
            return(X)
        prob <- lambda[X]/lmax
        u <- runif(X$n)
        retain <- (u <= prob)
        X <- X[retain, ]
        return(X)
    }
    stop(paste(sQuote("lambda"), "must be a constant, a function or an image"))
} 
