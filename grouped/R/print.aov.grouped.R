"print.aov.grouped" <-
function(x, ...){
    if(!inherits(x, "aov.grouped"))
        stop("Use only with 'aov.grouped' objects.\n")
    p.val <- x$p.value
    p.val <- if(p.val < 0.0001) "<0.0001" else if(p.val > 0.9999) ">0.9999" else as.character(round(p.val, 4))
    dat <- data.frame("logLik" = c(x$L0, x$L1), "AIC" = c(x$AIC0, x$AIC1), "BIC" = c(x$BIC0, x$BIC1), 
                        "L.ratio" = c("", as.character(round(x$L01, 3))), "df" = c("", x$df1 - x$df0), 
				"p.value" = c("", p.val), row.names = c(x$name0, x$name1))
    cat("\nLikelihood Ratio Test\n")
    print(dat)
    cat("\n")
    invisible(x)
}

