gamboostLSS_intern <- function(..., fun = c("check", "do_trace")) {

    fun <- match.arg(fun)
    do.call(fun, list(...))
}
