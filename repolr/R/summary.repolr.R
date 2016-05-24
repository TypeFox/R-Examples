summary.repolr <-
function (object, digits = 4, robust.var = TRUE, ...)  {

    # summary data
    coef <- object$coefficients
    cnames <- names(coef)
    coef <- matrix(rep(round(coef, digits = digits), 4), ncol = 4)
    if(robust.var == TRUE){
      se <- sqrt(diag(object$robust.var))
      dimnames(coef) <- list(cnames, c("coeff", "se.robust", "z.robust", "p.value"))
    } else {
      se <- sqrt(diag(object$naive.var))
      dimnames(coef) <- list(cnames, c("coeff", "se.naive", "z.naive", "p.value"))
    }
    coef[,2] <- round(se, digits = digits)
    coef[,3] <- round(coef[,1] / coef[,2], digits = digits)
    coef[,4] <- round(2 * pnorm(abs(coef[, 3]), lower.tail = F), digits = digits)

    # output
    summary <- list()
    summary$call <- object$call
    summary$title <- object$title
    summary$coefficients <- coef
    summary$alpha <- round(object$alpha, digits = digits)
    summary$times <- object$times
    summary$fixed <- object$fixed
    summary$corr.mod <- object$corr.mod
    summary$po.test <- object$po.test
    summary$po.test$po.stat <- round(object$po.test$po.stat, digits = digits)
    summary$po.test$po.chi <- round(object$po.test$po.chi, digits = digits)
    attr(summary,"class") <- "summary.repolr"
    summary
}
