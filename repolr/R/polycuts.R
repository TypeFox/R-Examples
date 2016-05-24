polycuts <-
function(object, digits = 3, robust.var = TRUE){

   # summary data
   if(is.null(object$poly$poly) == FALSE){
    coef <- object$coefficients[1:(object$poly$poly + 1)]
    cnames <- names(coef)
    coef <- matrix(rep(round(coef, digits = digits), 2), ncol = 2)
    pcnames <- paste("cuts",paste(as.character(1:object$categories)[-object$categories],
                             as.character(1:object$categories)[-1], sep = "|"), sep = "")
    pcoef <- object$poly$polycuts$coeff
    pcoef <- matrix(rep(round(pcoef, digits = digits), 2), ncol = 2)
    dimnames(pcoef) <- list(pcnames, c("coeff", "se.robust"))
    if(robust.var == TRUE){
      se <- sqrt(diag(object$robust.var))
      dimnames(coef) <- list(cnames, c("coeff", "se.robust"))
      coef[,2] <- round(se[1:(object$poly$poly + 1)], digits = digits)
      dimnames(pcoef) <- list(pcnames, c("coeff", "se.robust"))
      pcoef[,2] <- round(object$poly$polycuts$robust, digits = digits)
    } else {
      se <- sqrt(diag(object$naive.var))
      dimnames(coef) <- list(cnames, c("coeff", "se.naive"))
      coef[,2] <- round(se[1:(object$poly$poly + 1)], digits = digits)
      dimnames(pcoef) <- list(pcnames, c("coeff", "se.naive"))
      pcoef[,2] <- round(object$poly$polycuts$naive, digits = digits)
    }
   } else {
    coef <- object$coefficients[1:(object$categories - 1)]
    cnames <- names(coef)
    coef <- matrix(rep(round(coef, digits = digits), 2), ncol = 2)
    if(robust.var == TRUE){
      se <- sqrt(diag(object$robust.var))
      dimnames(coef) <- list(cnames, c("coeff", "se.robust"))
    } else {
      se <- sqrt(diag(object$naive.var))
      dimnames(coef) <- list(cnames, c("coeff", "se.naive"))
    }
    coef[,2] <- round(se[1:(object$categories - 1)], digits = digits)
   }

   # output coefficients and standard errors
   cat("\n")
   cat(object$title,"\n")
   cat("\nCall:\n", paste(deparse(object$call), sep = "\n", collapse = "\n"), 
        "\n\n", sep = "")
    if (length(coef)) {
        cat("Coefficients:", "\n")
        print.default(format(coef), print.gap = 2L, quote = FALSE)
    }
    else {cat("No Coefficients", "\n")}
    cat("\n")
   if(is.null(object$poly$poly) == FALSE){
    if (length(pcoef)) {
        cat("Cut-point Estimates:", "\n")
        print.default(format(pcoef), print.gap = 2L, quote = FALSE)
        cat("\n")
        cat("Polynomial Order: ", object$poly$poly, "\n", sep = "")
    }
    else {cat("No Cut-point Estimates", "\n")}
    cat("\n")
    invisible(list(coef = coef, poly = pcoef, order = object$poly$poly))
   } else {
    invisible(list(coef = coef))
   }
}
