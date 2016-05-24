## Ensure consistent "diss.." class --- make "namespace-private-global !
dissiCl <- c("dissimilarity", "dist")

if((Rv <- getRversion()) < "3.1.0") {
  anyNA <- function(x) any(is.na(x))
  ## if(Rv < "3.0.0") {
  ##     rep_len <- function(x, length.out) rep(x, length.out=length.out)
  ##     ## if(Rv < "2.15")
  ##     ##     paste0 <- function(...) paste(..., sep = '')
  ## }
}; rm(Rv)

## Not exported, and only used because CRAN checks must be faster
doExtras <- function() {
    interactive() || nzchar(Sys.getenv("R_CLUSTER_CHECK_EXTRA")) ||
        identical("true", unname(Sys.getenv("R_PKG_CHECKING_doExtras")))
}

