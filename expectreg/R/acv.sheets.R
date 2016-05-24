acv.sheets <-
function (penalty, yy, B, pp, DD, nb, center) 
{
    penalty = matrix(abs(penalty), ncol = 2)
    print(penalty)
    np = length(pp)
    m = length(yy)
    ynp <- rep(yy, np)
    ps <- rep(pp, each = m)
    vv <- log(ps/(1 - ps))
    w <- runif(m * np)
    p2f.new <- pspfit2d.new(B, DD, ps, ynp, w, penalty[, 1], 
        penalty[, 2], center)
    score = w * (ynp - p2f.new$fit)^2/(1 - hat(p2f.new$hatma)[1:length(ynp)])^2
    mean(score[which(is.finite(score))], na.rm = TRUE)
}
