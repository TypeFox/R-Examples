fp <- function(x, df = 4, select = NA, alpha = NA, scale=TRUE)
{
# Version 1.0.1	09 dec 03
#
    name <- deparse(substitute(x))
    attr(x, "df") <- df
    attr(x, "alpha") <- alpha
    attr(x, "select") <- select
    attr(x, "scale") <- scale
    attr(x, "name") <- name
    x
}
