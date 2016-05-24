
library(qtutils)

## test utility

showAndClose <- function(x, wait = 3, title = deparse(substitute(x)))
{
    x$show()
    x$windowTitle <- as.character(title)
    Sys.sleep(wait)
    x$close()
    invisible()
}    

showAndClose(data.browse())
showAndClose(qdataview(iris))

## showAndClose(qpager())



