env2stack <-
function (foldenv, vars = "", fext) 
{
    extens <- paste("*.", fext, sep = "")
    files <- list.files(path = foldenv, pattern = extens, full.names = FALSE)
    n <- nchar(vars)
    if (any(n) == 0) {
        fn1 <- file.path(dirname(foldenv), basename(foldenv), 
            files[1], fsep = .Platform$file.sep)
        fe <- file.exists(fn1)
        if (fe == FALSE) 
            stop("file does not exist")
        s <- raster(fn1)
        for (i in 2:length(files)) {
            fn1 <- file.path(dirname(foldenv), basename(foldenv), 
                files[i], fsep = .Platform$file.sep)
            fe <- file.exists(fn1)
            if (fe == FALSE) 
                stop("file does not exist")
            bb <- raster(fn1)
            s <- stack(s, bb)
        }
    }
    else {
        v1 <- paste(vars[1], fext, sep = ".")
        fn1 <- file.path(dirname(foldenv), basename(foldenv), 
            v1, fsep = .Platform$file.sep)
        fe <- file.exists(fn1)
        if (fe == FALSE) 
            stop("file does not exist")
        s <- raster(fn1)
        if (length(vars) > 1) {
            for (j in 2:length(vars)) {
                fn <- paste(vars[j], fext, sep = ".")
                fn1 <- file.path(dirname(foldenv), basename(foldenv), 
                  fn, fsep = .Platform$file.sep)
                fe <- file.exists(fn1)
                if (fe == FALSE) 
                  stop("file does not exist")
                bb <- raster(fn1)
                s <- stack(s, bb)
            }
        }
    }
    return(s)
}
