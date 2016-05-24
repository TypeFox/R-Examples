geom_mean <- function (x, na.rm = TRUE) {
    ans <- exp(mean(log(x), na.rm = TRUE))
    ans
}
