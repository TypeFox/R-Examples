"glm.diag" <-
function (glmfit) 
{
    if (is.null(glmfit$prior.weights)) 
        w <- rep(1, length(glmfit$residuals))
    else w <- glmfit$prior.weights
    sd <- sqrt(summary(glmfit)$dispersion)
    dev <- residuals(glmfit, type = "deviance")/sd
    pear <- residuals(glmfit, type = "pearson")/sd
    h <- rep(0, length(w))
    h[w != 0] <- lm.influence(glmfit)$hat
    p <- glmfit$rank
    rp <- pear/sqrt(1 - h)
    rd <- dev/sqrt(1 - h)
    cook <- (h * rp^2)/((1 - h) * p)
    res <- sign(dev) * sqrt(dev^2 + h * rp^2)
    list(res = res, rd = rd, rp = rp, cook = cook, h = h, sd = sd)
}

