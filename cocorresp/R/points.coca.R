`points.coca` <- function(x, display = c("sites","species"),
                          which = c("response","predictor"),
                          choices = c(1,2), scaling = FALSE,
                          select, ...) {
    if(length(display) > 1) {
        warning("Only one set of scores can be plotted at a time.")
    }
    display <- match.arg(display)
    ## what are we plotting, response or predictor?
    which <- match.arg(which)
    ## and map to X and Y for extraction
    WHICH <- ifelse(which == "response", "Y", "X")
    ## should the scores be rescaled - only for species though
    if(is.logical(scaling))
        scaling <- ifelse(scaling, 2, 1)
    ## need two and only two axes to plot
    if(length(choices) != 2)
        stop("Exactly two axes should be specified in `choices`")
    ## get the scores for plotting
    scrs <- scores(x, choices = choices, display = display,
                   scaling = scaling)
    ## then extract the response or predictor scores
    scrs <- lapply(scrs, `[[`, WHICH)
    ## alter the names of scrs so xy.coords knows what to do
    scrs <- scrs[[1]]
    colnames(scrs) <- c("x","y")
    ## draw the scores with points()
    points(scrs, ...)
    invisible(scrs)
}
