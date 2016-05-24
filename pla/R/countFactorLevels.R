
.is.digits <- function(x) all(nchar(gsub("\\d",   "", x, perl = TRUE)) == 0)
.is.llower <- function(x) all(nchar(gsub("[a-l]", "", x, perl = TRUE)) == 0)
.is.lower  <- function(x) all(nchar(gsub("[a-z]", "", x, perl = TRUE)) == 0)
.is.upper  <- function(x) all(nchar(gsub("[A-Z]", "", x, perl = TRUE)) == 0)
.is.uupper <- function(x) all(nchar(gsub("[M-Z]", "", x, perl = TRUE)) == 0)

.classFactor <- function(factor) {
    class <- 99
    if (all(nchar(factor) == 1)) {
       if (.is.digits(factor)) class <- 0
       if (.is.llower(factor)) class <- 1
       if (.is.lower (factor)) class <- 2
       if (.is.upper (factor)) class <- 21
       if (.is.uupper(factor)) class <- 22
    } else {
       if (.is.digits(factor)) class <- 10
    }
    return(class)
}

.countBlocks <- function(m) {
    m1 <- m[1]
    x <- 0
    count <- 1
    for (i in 2:length(m))
        if (m[i] == m1) {
            count <- count + x
            x <- 0
        } else
            x <- 1
    return(count)
}

.countRun <- function(m, warnings = TRUE) {
    r <- NULL
    j <- m[1]
    count <- 1
    for (i in 2:length(m))
        if (m[i] == j)
            count <- count + 1
        else {
            r <- c(r, count)
            j <- m[i]
            count <- 1
        }
    r <- c(r, count)
    if (min(r) != max(r))
        if (FALSE & warnings)
            warning(paste0("Unbalanced factor, length of 'run': ",
                           paste0(m, collapse = "/"), " - ",
                           paste0(r, collapse = "/")),
                    call. = FALSE)
    return(r)
}

.countFactorLevels <- function(factor, warnings = TRUE) {
    name <- paste0(unique(factor), collapse = "/")
    r <- c(before = .countBlocks(factor),
           inter = length(unique(factor)),
           after = max(unique(.countRun(factor, warnings = warnings))))
    tab <- table(factor)
    if (min(tab) != max(tab))
        if (warnings)
            warning(paste0("Unbalanced factor (", name, "): ",
                           paste0(factor, collapse = "/"), " - ",
                           paste0(tab, collapse = "/")),
                    call. = FALSE)
    p <- prod(r)
    if (p != length(factor))
        if (FALSE & warnings)
            warning(paste0("Unbalanced factor (", name, "): ",
                           paste0(factor, collapse = "/")),
                    call. = FALSE)
        else
            ## Patch for differente, but some identical,
            ## doses between samples:
            r["inter"] <- length(factor) / (r["before"] * r["after"])
    r <- c(r, prod = p, class = .classFactor(factor))
    return(list(values = r, name = name))
}
