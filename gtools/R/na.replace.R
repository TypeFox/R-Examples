na.replace <- function(x, replace)
{
  x[is.na(x)] <- replace
  x
}
