print.summary.repolr <-
function(x, digits = 4, robust.var = TRUE, ...){

   # output coefficients, standard errors, po-test and correlation
   cat("\n")
   cat(x$title,"\n")
   cat("\nCall:\n", paste(deparse(x$call), sep = "\n", collapse = "\n"), 
        "\n\n", sep = "")
    if (length(x$coef)) {
        cat("Coefficients:", "\n")
        print.default(format(x$coef), print.gap = 2L, 
            quote = FALSE)
    }
    else {cat("No coefficients")}
    cat("\n")
    if(length(x$times) > 1) {
     cat("Correlation Structure: ",x$corr.mod, "\n")
     if(x$fixed == FALSE){
      cat("Estimated Correlation: ", x$alpha, "\n")
     } else {
      cat("Fixed Correlation: ", x$alpha, "\n")
     }
    }
    if(is.na(x$po.test$po.stat) == FALSE){
     cat("PO Score Test: ", x$po.test$po.stat, " (d.f. = ",
            as.integer(x$po.test$po.df), " and p.value = ",
            x$po.test$po.chi, ")", "\n", sep = "")
    }
    cat("\n")
    invisible(x)
}
