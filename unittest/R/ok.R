ok <- function(
    test,
    description
) {
    if(missing(description)) description <- strtrim(paste0(deparse(substitute(test)), collapse = " "), 60)
    if(! is.character(description) || length(description) > 1) stop('\'description\' must be of type \'chr\' and not a vector.')

    result <- tryCatch(test, error = function(e) e)

    outcome <- data.frame()
    if(identical(result, TRUE) ) {
        outcome <- data.frame(
            status = TRUE,
            output = paste('ok -', description, collapse = " "),
            stringsAsFactors = FALSE
        )
    }
    else if(inherits(result,'error')) {
        outcome <- data.frame(
            status = FALSE,
            output = paste(
                paste('not ok -', description, collapse = " "),
                paste("# Test resulted in error:", result$message, collapse = " "),
                paste("#  ->", result$call, collapse = "\n"),
                sep = "\n", collapse = "\n"
            ),
            stringsAsFactors = FALSE
        )
    }
    else if(is.character(result)) {
        outcome <- data.frame(
            status = FALSE,
            output = paste(
                paste('not ok -', description, collapse = " "),
                "# Test returned non-TRUE value:",
                paste("#", unlist(strsplit(result, split = "\n")), collapse = "\n"),
                sep = "\n", collapse = "\n"
            ),
            stringsAsFactors = FALSE
        )
    }
    else {
        outcome <- data.frame(
            status = FALSE,
            output = paste(
                paste('not ok -', description, collapse = " "),
                "# Test returned non-TRUE value:",
                paste("#", capture.output( print(result) ), collapse = "\n"),
                sep = "\n", collapse = "\n"
            ),
            stringsAsFactors = FALSE
        )
    }
    assign_outcome(outcome)
    rv <- paste0(outcome['output'], "\n")
    cat(rv)
    invisible(result)
}
