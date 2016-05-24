kendallSeasonalTrendTest.formula <-
function (y, data = NULL, subset, na.action = na.pass, ...) 
{
    if (missing(y) || (length(y) != 3L)) 
        stop("formula missing or incorrect")
    m <- match.call(expand.dots = FALSE)
    if (is.matrix(eval(m$data, parent.frame()))) 
        m$data <- as.data.frame(data)
    m$... <- NULL
    m$formula <- m$y
    m$y <- NULL
    m$na.action <- na.action
    requireNamespace("stats", quietly = TRUE)
    m[[1L]] <- as.name("model.frame")
    mf <- eval(m, parent.frame())
    ncol.mf <- ncol(mf)
    if (ncol.mf != 3) 
        stop("Incorrect formula; it must be of the form y ~ season + year")
    response <- attr(attr(mf, "terms"), "response")
    season <- (1:ncol.mf)[-response][1]
    year <- (1:ncol.mf)[-c(response, season)]
    arg.list <- list(y = mf[, response], season = mf[, season], 
        year = mf[, year])
    names.mf <- names(mf)
    names.list <- list(data.name = names.mf[response], season.name = names.mf[season], 
        year.name = names.mf[year])
    dot.list <- list(...)
    match.vec <- pmatch(names(dot.list), c("data.name", "season.name", 
        "year.name"), nomatch = 0)
    if (length(match.vec) == 0 || all(match.vec == 0)) 
        arg.list <- c(arg.list, names.list, dot.list)
    else arg.list <- c(arg.list, names.list[-match.vec], dot.list)
    if (!missing(data)) 
        arg.list$parent.of.data <- deparse(substitute(data))
    if (!missing(subset)) 
        arg.list$subset.expression <- deparse(substitute(subset))
    do.call(kendallSeasonalTrendTest.default, arg.list)
}
