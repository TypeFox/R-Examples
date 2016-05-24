print.anesrake <-
function(x, ...) {
    out <- cbind(x$caseid, x$weightvec)
    colnames(out) <- c("caseid", "weightvec")
    as.data.frame(out)
    print(out)
}

