password <- function(n = 8, numbers = TRUE, case = TRUE,
                     special = c("?", "!", "&", "%", "$")) {
    resample <- function(x, ...)
        x[sample.int(length(x), ... , replace = TRUE)]
    from <- letters
    if (numbers)
        from <- c(from, 0:9)
    if (!identical(FALSE, special))
        from <- c(from, special)
    if (case) from <- c(from, LETTERS)
    res <- resample(from, n)
    paste(res, collapse = "")
}
