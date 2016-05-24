Log <- function (description = NULL)
{
	if (!exists(".Log", envir = .GlobalEnv, inherits = FALSE))
		createLog(description = description)
	return(get(".Log", envir = .GlobalEnv, inherits = FALSE))
}

createLog <- function (description = NULL, deleteExisting = FALSE)
{
    ## Create a log consisting in an environment with class svSuiteData
    if (isTRUE(deleteExisting) && exists(".Log", envir = .GlobalEnv,
        inherits = FALSE)) rm(.Log, envir = .GlobalEnv)
    if (!exists(".Log", envir = .GlobalEnv, inherits = FALSE)) {
        .Log <<- structure(new.env(parent = .GlobalEnv),
            class = c("svSuiteData", "environment"))
        ## Add general informations and variables to it
        .Log$.time <- Sys.time()            # Creation time of the log
        .Log$.R.version <- R.version        # R version
        .Log$.sessionInfo <- sessionInfo()  # Information about current session
        .Log$.description <- description    # Optional description of this log
        ## Create ..xxx variables used for test context
		## Note: never delete or put NULL in these variables, use "" instead
		.Log$..Unit <- ""
		.Log$..Msg <- ""
		.Log$..Obj <- ""
		.Log$..File <- ""
		.Log$..Tag <- ""
		## Create .lastTest that contains details from last check...()
        naChr <- as.character(NA)
        .Log$.lastTest <- structure(
            data.frame(msg = naChr, call = naChr,
                timing = as.numeric(NA), kind = .kind(NA), res = naChr,
                obj = naChr, file = naChr, tag = naChr,
                stringsAsFactors = FALSE),
            class = c("svTestData", "data.frame"))
		## Create .lastSuite with an empty list of test units to run
		.Log$.lastSuite <- list()
    }
}

clearLog <- function ()
{
	if (exists(".Log", envir = .GlobalEnv, inherits = FALSE)) {
		rm(list = ".Log", envir = .GlobalEnv)
		return(invisible(TRUE))
	} else return(invisible(FALSE))
}

errorLog <- function (stopit = TRUE, summarize = TRUE)
{
	.Log <- Log()
	Res <- table(stats(.Log)$kind)
	if (isTRUE(stopit) && any(Res[2:3] > 0)) {
		if (isTRUE(summarize)) summary(.Log)
		msg <- paste(Res[2], "failure(s) and", Res[3], "error(s)")
		stop(msg)
	} else if (interactive()) {
		cat("Summary statistics on all tests run:\n")
		print(Res)
	}
	return(invisible(Res))
}

lastTest <- function ()
{
    ## Return a svTestData object with data from last recorded test
	Log()$.lastTest
}

lastSuite <- function ()
{
    ## Return data about last suite run
	Log()$.lastSuite
}
