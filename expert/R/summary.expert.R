### ===== expert =====
###
### Summary method for object of class "expert"
###
### AUTHORS: Vincent Goulet <vincent.goulet@act.ulaval.ca>,
###          Mathieu Pigeon <mathieu.pigeon.3@ulaval.ca>

summary.expert <- function(object, ...)
{
    class(object) <- c("summary.expert", class(object))
    object
}

print.summary.expert <- function(x, ...)
{
    cat("Call:\n")
    print(attr(x, "call"))
    NextMethod()                       # print.expert()
    if (is.null(x$alpha)) cat("\n")    # alpha not displayed
    cat(" Number of experts: ", x$nexp,
        ",\tNumber of seed variables: ", x$nseed, "\n", sep = "")
    q <- x$quantiles
    qa <- unlist(attributes(q))
    cat(" Quantiles: ",
        if (qa["set.lower"]) "0*, ",
        paste(x$quantiles, collapse = ", "),
        if (qa["set.upper"]) ", 1*",
        if (any(qa)) "\t(* were set)",
        "\n", sep = "")
    invisible(x)
}
