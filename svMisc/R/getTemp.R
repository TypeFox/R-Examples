getTemp <- function (x, default = NULL, mode = "any", item = NULL)
{
    if (is.null(item)) Mode <- mode else Mode <- "any"
    if  (exists(x, envir = TempEnv(), mode = Mode, inherits = FALSE)) {
        dat <- get(x, envir = TempEnv(), mode = Mode, inherits = FALSE)
        if (is.null(item)) return(dat) else {
            item <- as.character(item)[1]
            if (inherits(dat, "list") && item %in% names(dat)) {
                dat <- dat[[item]]
                if (mode != "any" && mode(dat) != mode) dat <- default
                return(dat)
            } else {
                return(default)
            }
        }
    } else { # Variable not found, return the default value
        return(default)
    }
}
