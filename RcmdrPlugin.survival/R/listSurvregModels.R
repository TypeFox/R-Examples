# last modified 7 December 2008 by J. Fox

listSurvregModels <-
function(envir=.GlobalEnv, ...) {
	objects <- ls(envir=envir, ...)
	if (length(objects) == 0) NULL
	else objects[sapply(objects,
						function(.x) "survreg" == (class(get(.x, envir=envir))[1]))]
}

survregP <- function() activeModelP() && class(get(ActiveModel()))[1] == 'survreg'