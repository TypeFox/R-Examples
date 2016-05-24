# last modified 2014-08-09 by J. Fox

summarySandwich <- function(model, ...){
    UseMethod("summarySandwich")
}

summarySandwich.lm <- function(model, type=c("hc3", "hc0", "hc1", "hc2", "hc4", "hac"), ...){
    s <- summary(model)
    c <- coef(s)
    type <- match.arg(type)
    v <- if (type != "hac") hccm(model, type=type, ...)
    else vcovHAC(model, ...)
    c[, 2] <- sqrt(diag(v))
    c[, 3] <- c[,1]/c[,2]
    c[, 4] <- 2*pt(abs(c[,3]), df=s$df[2], lower.tail=FALSE)
    colnames(c)[2] <- paste("Std.Err(", type, ")", sep="")
    s$coefficients <- c
    coefs <- names(coef(model))
    coefs <- coefs[coefs != "(Intercept)"]
    h <- linearHypothesis(model, coefs, vcov.=v)
    s$fstatistic <- c(value=h$F[2], numdf=length(coefs), dendf=s$df[2])
    s
}