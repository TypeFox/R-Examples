.pastenames <- function (...){
  if (nargs () == 1L & is.list (..1))
    dots <- ..1
  else
    dots <- list (...)

  names <- names (dots)
  names <- sapply (names,
                   function (x){
                     if (nchar (x) > 0L)
                       sprintf ("%s = ", x)
                     else
                       ""
                     })

  paste (names, dots, collapse = ", ", sep = "")
}
