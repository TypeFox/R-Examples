
## no longer needed since we sent the PR to cpptoml
## .sort <- function(s) {
##     if (length(s) > 0) {
##         for (elem in seq_along(s)) {
##             if (is.list(s[[elem]])) {
##                 s[[elem]] <- .sort(s[[elem]])
##             }
##         }
##         if (length(names(s)) > 0) {
##             s <- s[order(names(s))]
##         }
##     }
##     s
## }
