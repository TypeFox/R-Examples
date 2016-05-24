kendallSeasonalTrendTest.data.frame <-
function (y, ...) 
{
    if (ncol(y) == 1 || !all(sapply(y, is.numeric))) {
        stop(paste("When 'y' is a data frame,", "it must have 2 or more columns and they must", 
            "all be numeric."))
    }
    season <- as.vector(col(y))
    year <- as.vector(row(y))
    arg.list <- list(y = as.vector(unlist(y)), season = season, 
        year = year)
    dots.list <- list(...)
    data.name <- deparse(substitute(y))
    names.list <- list(data.name = data.name, season.name = paste("Columns of", 
        data.name), year.name = paste("Rows of", data.name))
    match.vec <- pmatch(names(dots.list), c("data.name", "season.name", 
        "year.name"), nomatch = 0)
    if (length(match.vec) == 0 || all(match.vec == 0)) 
        arg.list <- c(arg.list, names.list, dots.list)
    else arg.list <- c(arg.list, names.list[-match.vec], dots.list)
    do.call("kendallSeasonalTrendTest.default", arg.list)
}
