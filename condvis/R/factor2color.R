factor2color <-
function (x, colors = NULL)
{
    x <- if (!is.factor(x)){
        warning("'x' has been coerced to a factor.")
        as.factor(x)
    } else x
    n <- nlevels(x)
    colors <- if (is.null(colors))
        if (requireNamespace("RColorBrewer", quietly = TRUE))
		    RColorBrewer::brewer.pal(n = max(n, 3L, na.rm = TRUE), name = 
                "Set3")[1L:n]
		else rainbow(n)
    else rep(colors, length.out = n)
    vapply(x, function(y) colors[levels(x) == as.character(y)], character(1L))
}
