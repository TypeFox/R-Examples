num.ident <- function(x,y) {
  x==y | is.nan(x) & is.nan(y) | is.na(x) & !is.nan(x) & is.na(y) & !is.nan(y)
}
