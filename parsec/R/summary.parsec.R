summary.parsec <-
function(object, ...) {
    
    object$inequality <- NA
    
    if (!is.null(object$profiles)) {
        res <- object$profiles$profiles    
    } else {
        res <- data.frame(tmp = 1:object$number_of_profiles)
    }
    
    res$weights <- object$prof_w
    res$"threshold" <- object$threshold
    res$"id. function" <- object$idn_f
    res$"average rank" <- apply(
        object$rank_dist,
        2,
        function(x) sum((1:object$number_of_profiles)*x)
    )
    res$"abs. severity" <- object$svr_abs
    res$"rel. severity" <- object$svr_rel
    res$"abs. wealth gap" <- object$wea_abs
    res$"rel. wealth gap" <- object$wea_rel
    
    if ("tmp" %in% names(res)) res$tmp <- NULL
    
    rownames(res) <- rownames(object$incidence)
    
    ord <- order(-res$"id. function", res$"average rank")
    res <- res[ord,]
    
    if (object$number_of_profiles < 11)
        print(res)
    else {
        print(res[1:5,])
        cat("\n...\n\n")
        print(res[object$number_of_profiles - 4:0,])
    }
    try(
#         cat("\nhead count ratio  = ", object$head_count_ratio, "\n",
        cat("poverty gap  = ", object$poverty_gap, "\n",
            "wealth gap   = ", object$wealth_gap, "\n", sep = "")
    , silent = TRUE)
    try(
    {
    if (!is.na(object$inequality)) {
        try(cat("inequality = ", object$inequality, sep = ""), silent = TRUE)
    } else {
        cat("inequality has not been evaluated")
    }
    cat("\n\nthis function returns a data.frame that summarize each profile\nyou can also summarize\n\n")
    }, silent = TRUE)
    invisible(res)
}
