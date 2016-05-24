## --- Auxiliary routines ---------------------------------------------------

## check for warnings but suppress warning message
expect_warning_suppress <- function (x, y) {
        opt.warn.save <- getOption("warn")
        options(warn=-1)
        expect_warning(x,y)
        options(warn=opt.warn.save)
}

## --------------------------------------------------------------------------
