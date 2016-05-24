`orgltext` <-
    function (object, text, display = "sites", choices = 1:3,
              adj = 0.5, col = "black", ...)
{
    x <- scores(object, display = display, choices = choices, 
                ...)
    if (missing(text)) 
        text <- rownames(x)
    ## colors
    if (is.factor(col))
        col <- as.numeric(col)
    col <- rep(col, length = nrow(x))
    rgl.texts(x[, 1], x[, 2], x[, 3], text, adj = adj,  col = col, ...)
    invisible()
}
