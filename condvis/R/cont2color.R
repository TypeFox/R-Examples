cont2color <-
function (x, xrange = NULL, breaks = NULL, colors = NULL)
{
    x <- if (!is.numeric(x)){
        as.numeric(x)
        warning("'x' has been coerced to numeric.")
    } else x
    xrange <- if (is.null(xrange))
        range(x)
    else xrange    
    breaks <- if(is.null(breaks))
	    11
    else breaks
	br <- c(min(x, min(xrange)) - 1, seq(min(xrange), max(xrange), length.out = 
        breaks - 1), max(x, max(xrange)) + 1)
    colors <- if (is.null(colors)){
        if (requireNamespace("RColorBrewer", quietly = TRUE))
		    RColorBrewer::brewer.pal(n = max(breaks, 3L, na.rm = TRUE), 
                name = "PiYG")
		else cm.colors(breaks)
    } else rep(colors, length.out = breaks)
    as.character(cut(x, br, labels = colors, include.lowest = TRUE))
}
