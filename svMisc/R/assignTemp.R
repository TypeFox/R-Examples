assignTemp <- function (x, value, replace.existing = TRUE)
    if (replace.existing || !exists(x, envir = TempEnv(), mode = "any",
		inherits = FALSE))
        assign(x, value, envir = TempEnv())
