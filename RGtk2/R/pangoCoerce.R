#as.PangoMatrix <-
#function(x)
#{
#    x <- as.struct(x, "PangoMatrix", c("xx", "xy", "yx", "yy", "x0", "y0"))
#    x[[1]] <- as.numeric(x[[1]])
#    x[[2]] <- as.numeric(x[[2]])
#    x[[3]] <- as.numeric(x[[3]])
#    x[[4]] <- as.numeric(x[[4]])
#    x[[5]] <- as.numeric(x[[5]])
#    x[[6]] <- as.numeric(x[[6]])
#
#    return(x)
#}

as.PangoRectangle <- function(x)
{
	x <- as.GdkRectangle(x)
	class(x) <- "PangoRectangle"
	x
}
