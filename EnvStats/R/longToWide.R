longToWide <-
function (x, data.var, row.var, col.var, row.labels = levels(factor(x[, 
    row.var])), col.labels = levels(factor(x[, col.var])), paste.row.name = FALSE, 
    paste.col.name = FALSE, sep = ".", check.names = FALSE, ...) 
{
    if (!((is.df <- is.data.frame(x)) || is.matrix(x)) || ncol(x) < 
        3) 
        stop("'x' must be a data frame or matrix with at least 3 columns")
    if (!is.df) 
        x <- data.frame(x)
    names.x <- names(x)
    all.char <- is.character(data.var) && is.character(row.var) && 
        is.character(col.var)
    all.num <- is.numeric(data.var) && is.numeric(row.var) && 
        is.numeric(col.var)
    if (length(data.var) != 1 || length(row.var) != 1 || length(col.var) != 
        1 || !(all.char || all.num) || length(unique(c(row.var, 
        col.var, data.var))) != 3) 
        stop(paste("'data.var', 'row.var', and 'col.var' must be", 
            "all character scalars or all numeric scalars, and must be unique"))
    if (all.char) {
        if (!all(c(row.var, col.var, data.var) %in% names.x)) 
            stop(paste("When 'row.var', 'col.var', and 'data.var' are character scalars,", 
                "they must match column names in 'x'"))
    }
    else {
        if (!all(c(row.var, col.var, data.var) %in% 1:ncol(x))) 
            stop(paste("When 'row.var', 'col.var', and 'data.var' are numeric scalars,", 
                "they must match column numbers in 'x'"))
    }
    row.vals <- x[, row.var]
    col.vals <- x[, col.var]
    if (any(is.na(row.vals)) || any(is.na(col.vals))) 
        stop("Missing values not allowed in the variables used to define rows and columns")
    if (length(unique(paste(row.vals, col.vals, sep = "."))) != 
        nrow(x)) 
        stop("The value combinations using 'row.var' and 'col.var' must be unique")
    if (!missing(row.labels) && !is.character(row.labels) || 
        ((lrl <- length(row.labels)) != length(unique(row.labels)) || 
            lrl != length(levels(factor(row.vals))))) 
        stop(paste("'row.labels' must be a character vector and the number of", 
            "unique values in 'row.labels' must match the number of unique values", 
            "in the column specified by 'row.var'"))
    if (!missing(col.labels) && !is.character(col.labels) || 
        ((lcl <- length(col.labels)) != length(unique(col.labels)) || 
            lcl != length(levels(factor(col.vals))))) 
        stop(paste("'col.labels' must be a character vector and the number of", 
            "unique values in 'col.labels' must match the number of unique values", 
            "in the column specified by 'col.var'"))
    if (paste.row.name || paste.col.name) {
        if (all.num) {
            row.var <- names.x[row.var]
            col.var <- names.x[col.var]
        }
        if (paste.row.name) 
            row.labels <- paste(row.var, row.labels, sep = sep)
        if (paste.col.name) 
            col.labels <- paste(col.var, col.labels, sep = sep)
    }
    ret.val <- tapply(x[, data.var], list(row.vals, col.vals), 
        FUN = function(x) x)
    dimnames(ret.val) <- list(as.character(row.labels), as.character(col.labels))
    if (is.df) 
        ret.val <- data.frame(ret.val, check.names = check.names, 
            ...)
    ret.val
}
