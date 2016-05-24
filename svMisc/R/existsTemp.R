existsTemp <- function (x, mode = "any")
    exists(x, envir = TempEnv(), mode = mode, inherits = FALSE)
