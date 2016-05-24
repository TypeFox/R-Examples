linesHyperb.lm <-
    function(object, c.prob = .95, confidence = FALSE,
             k = if(confidence) Inf else 1,
             col = 2, lty = 2, do.abline = TRUE)
{
    n <- length(Res <- residuals(object))
    df <- object $ df.resid
    s2 <- sum(Res^2)/df
    s <- sqrt(s2)
    if(is.null(R <- object $ R)) ## in R
        R <- qr.R(object $ qr)
    Xm <- R[1,2]/R[1,1] # = mean(x_i) : (X'X)[1,] = (R'R)[1,] = [n  sum_{x_i}]
    ##-- S_{xx} = sum_i{(x_i - mean(x_i))^2} : you can prove this: (R'R) = ...
    S.xx <- R[2,2]^2

    ux <- par("usr")[1:2]
    d.xs <- data.frame(x = xs <- seq(ux[1],ux[2], length = 100))
    names(d.xs) <-  attr(object$terms,"term.labels") #-- proper x-variable name
    ys <- predict(object, new = d.xs)
    pred.err <- qt(1-(1-c.prob)/2, df) * s * sqrt(1/k + 1/n + (xs-Xm)^2/S.xx)
    o.p <- par(err=-1)
    on.exit(par(o.p))
    if(do.abline)
        abline(object)
    lines(xs, ys - pred.err, col=col, lty=lty)
    lines(xs, ys + pred.err, col=col, lty=lty)
}

