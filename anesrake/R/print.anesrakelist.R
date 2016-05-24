print.anesrakelist <-
function(x, ...) {
    out <- cbind(x$caseid, x$weightvec)
    colnames(out) <- c("caseid", "weightvec")
    out
}

