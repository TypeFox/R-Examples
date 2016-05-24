write.dat <-
function (x, path) 
{
    if (is.na(dim(x)[3]) == FALSE) {
        ifelse(is.null(dimnames(x)[[3]]) == TRUE, lb <- 1:dim(x)[3], 
            lb <- dimnames(x)[[3]])
    }
    dir.create(path, showWarnings = FALSE)
    if (is.na(dim(x)[3]) == FALSE) {
        for (i in 1:dim(x)[3]) {
            filename <- sprintf(paste(substitute(x), "_", lb[i], 
                ".dat", sep = ""), i)
            pathname <- file.path(path, filename)
            obj <- x[, , i]
            utils::write.table(obj, file = pathname, row.names = FALSE, 
                col.names = FALSE)
        }
        rm(i)
    }
    else if (is.na(dim(x)[3]) == TRUE) {
        pathname <- file.path(path, paste(deparse(substitute(x)), 
            "dat", sep = "."))
        utils::write.table(x, file = pathname, row.names = FALSE, col.names = FALSE)
    }
}
