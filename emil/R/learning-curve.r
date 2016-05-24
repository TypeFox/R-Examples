#' Learning curve analysis
#' 
#' This function studies the change in permformance as the sizes of the training
#' set is varied. In case the studied modeling procedures cannot
#' produce models on the smallest training sets, please use
#' \code{.return_error=TRUE} (see \code{\link{evaluate}}.
#' 
#' @param procedure \code{\link{modeling_procedure}}.
#' @param x Dataset descriptors.
#' @param y Response.
#' @param test_fraction Fraction of dataset to hold out, i.e. use as test set. Defaults
#'   20 logarithmically distributed values ranging from all but 5 observations
#'   per class in the largest test set to only 5 observations per class in
#'   the smallest test set.
#' @param nfold How many holdout folds that should be calculated.
#' @param ... Sent to \code{\link{evaluate}}.
#' @param .verbose Whether to print an activity log. Set to \code{-1} to
#'   also suppress output generated from the procedure's functions.
#' @examples
#' options(emil_max_indent=3)
#' lc <- learning_curve(c(Linear="lda", Quadratic="qda"),
#'                      iris[-5], iris$Species, test_fraction=7:3/10)
#' plot(lc)
#' @references Richard O Duda, Peter E Hart, and David G Stork. Pattern
#'   Classification. Wiley, 2nd edition, 2000. ISBN 978-0-471-05669-0.
#' @author Christofer \enc{Bäcklin}{Backlin}
#' @export
learning_curve <- function(procedure, x, y, test_fraction, nfold=100, ..., .verbose=TRUE){
    if(missing(test_fraction)){
        min.fraction <- switch(class(y),
            factor = min(5/table(y)),
            5/length(y)
        )
        if(min.fraction >= .5) stop("Dataset is too small to do a meaningful learning curve analysis.")
        test_fraction <- 1-exp(seq(log(min.fraction), log(1-min.fraction), length.out=20))
        test_fraction <- test_fraction[!duplicated(round(test_fraction*length(y)))]
    }

    resamples <- lapply(test_fraction, function(f) resample("holdout", y, test_fraction=f, nfold=nfold))
    counter <- 0
    log_message(.verbose, "Learning curve analysis")
    result <- lapply(resamples, function(r){
        counter <<- counter + 1
        log_message(indent(.verbose, 1), "Test set test_fraction %i of %i (%.4g)",
                  counter, length(test_fraction), test_fraction[counter])
        evaluate(procedure, x, y, ..., resample=r, .verbose=indent(.verbose, 2))
    })
    structure(list(n = length(y), test_fraction = test_fraction, result=result),
              class="learning_curve")
}

#' Plot results from learning curve analysis
#'
#' @method plot learning_curve
#' @param x Results from \code{\link{learning_curve}}.
#' @param ... Ignored, kept for S3 consistency.
#' @param summaries Named list of summary functions that can reduce a vector of
#'   performance estimates to a single quantity.
#' @return A \code{\link[ggplot2]{ggplot}} object.
#' @author Christofer \enc{Bäcklin}{Backlin}
#' @importFrom dplyr summarize_
#' @importFrom ggplot2 ggplot aes_string geom_point geom_line xlab facet_wrap
#' @export
plot.learning_curve <- function(x, ..., summaries=list(mean = mean, `95-percentile`=function(x) quantile(x, .95))){
    if(is_multi_procedure(x$result[[1]])){
        plot.data <- x$result %>%
            select(test_fraction = TRUE, fold = TRUE, method = TRUE, performance = "error")
    } else {
        plot.data <- x$result %>%
            select(test_fraction = TRUE, fold = TRUE, performance = "error") %>%
            mutate_(method = NA)
    }
    plot.data$test_fraction <- x$test_fraction[plot.data$test_fraction]

    data.summary <- do.call(rbind, Map(function(method, fun){
        data.frame(summaries=method,
            summarize_(group_by_(plot.data, "test_fraction", "method"),
                      performance = ~fun(performance)))
    }, names(summaries), summaries))
    p <- ggplot(plot.data, aes_string(x="1-test_fraction", y="performance")) + 
        geom_point(colour="grey") +
        geom_line(data=data.summary, aes_string(colour="summaries")) +
        xlab("Relative training set size")
    if(is_multi_procedure(x$result[[1]]))
        p <- p + facet_wrap(~method)
    p
}
