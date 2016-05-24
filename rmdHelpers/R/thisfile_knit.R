thisfile_knit <-
function() {
  # borrowed from https://github.com/krlmlr/kimisc/blob/master/R/thisfile.R
  if (requireNamespace("knitr"))
    return (knitr::current_input())
  
  NULL
}
