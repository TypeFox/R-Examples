.clean_up_dependencies2 <- function (x, installed, available)
{

    .split_op_version <- function (x)
        {
            pat <- "^([^\\([:space:]]+)[[:space:]]*\\(([^\\)]+)\\).*"
            x1 <- sub(pat, "\\1", x)
            x2 <- sub(pat, "\\2", x)
            if (x2 != x1) {
                pat <- "[[:space:]]*([[<>=!]+)[[:space:]]+(.*)"
                version <- sub(pat, "\\2", x2)
                if (!grepl("^r", version))
                    version <- package_version(version)
                list(name = x1, op = sub(pat, "\\1", x2), version = version)
            }
            else list(name = x1)
        }

    .split_dependencies <- function(x) {
        .split2 <- function(x) {
            x <- sub("[[:space:]]+$", "", x)
            x <- unique(sub("^[[:space:]]*(.*)", "\\1", x))
            names(x) <- sub("^([[:alnum:].]+).*$", "\\1", x)
            x <- x[names(x) != "R"]
            x <- x[nzchar(x)]
            x <- x[!duplicated(names(x))]
            lapply(x, .split_op_version)
        }
        if (!any(nzchar(x)))
            return(list())
        unlist(lapply(strsplit(x, ","), .split2), FALSE, FALSE)
    }
    x <- x[!is.na(x)]
    if (!length(x))
        return(list(character(), character()))
    xx <- .split_dependencies(x)
    if (!length(xx))
        return(list(character(), character()))
    pkgs <- installed[, "Package"]
    have <- sapply(xx, function(x) {
        if (length(x) == 3L) {
            if (!x[[1L]] %in% pkgs)
                return(FALSE)
            if (x[[2L]] != ">=")
                return(TRUE)
            current <- as.package_version(installed[pkgs == x[[1L]],
                "Version"])
            target <- as.package_version(x[[3L]])
            eval(parse(text = paste("any(current", x$op, "target)")))
        }
        else x[[1L]] %in% pkgs
    })
    xx <- xx[!have]
    if (!length(xx))
        return(list(character(), character()))
    pkgs <- row.names(available)
    canget <- miss <- character()
    for (i in seq_along(xx)) {
        x <- xx[[i]]
        if (length(x) == 3L) {
            if (!x[[1L]] %in% pkgs) {
                miss <- c(miss, x[[1L]])
                next
            }
            if (x[[2L]] != ">=") {
                canget <- c(canget, x[[1L]])
                next
            }
            current <- as.package_version(available[pkgs == x[[1L]],
                "Version"])
            target <- as.package_version(x[[3L]])
            res <- eval(parse(text = paste("any(current", x$op,
                "target)")))
            if (res)
                canget <- c(canget, x[[1L]])
            else miss <- c(miss, paste0(x[[1L]], " (>= ", x[[3L]],
                ")"))
        }
        else if (x[[1L]] %in% pkgs)
            canget <- c(canget, x[[1L]])
        else miss <- c(miss, x[[1L]])
    }
    list(canget, miss)
}
