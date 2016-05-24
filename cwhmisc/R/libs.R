libs <- function (Lib)
{
    if (missing(Lib))
        print(.packages(all.available = TRUE), q = FALSE)
    else eval(parse(text = paste("library(help=",
as.character(substitute(Lib)),")")))
}
