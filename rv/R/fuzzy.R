
is.fuzzy <- function (x) {
  UseMethod("is.fuzzy")
}

is.fuzzy.rv <- function (x) {
  # NAME
  #  is.fuzzy - Is a Vector Component Logical But Random
  # 
  component.is.logical <- rvsimapply(x, is.logical)
  component.prop <- rvmean(x)
  (component.is.logical & component.prop>0 & component.prop<1)
}

is.fuzzy.default <- function (x)
{
  return(FALSE)
}
