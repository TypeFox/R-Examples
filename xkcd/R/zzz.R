## Emilio Torres Manzanera
## University of Oviedo
## Time-stamp: <2014-04-20 dom 17:01 emilio on emilio-Satellite-P100>
## ============================================================
    
## Borrowed from dplyr
## Otherwise, the print options of trunc_mat fails
.onAttach <- function(libname, pkgname) {
  op <- options()
  op.dplyr <- list(
    dplyr.show_sql = FALSE,
    dplyr.explain_sql = FALSE,
    dplyr.strict_sql = FALSE,
    dplyr.print_min = 10L,
    dplyr.print_max = 100L
  )
  toset <- !(names(op.dplyr) %in% names(op))
  if(any(toset)) options(op.dplyr[toset])

  invisible()
}
