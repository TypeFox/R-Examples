progress <- function (value, max.value = NULL, progress.bar = FALSE, char = "|",
init = (value == 0), console = TRUE, gui = TRUE)
{
	if (!is.numeric(value))
		stop("'value' must be numeric!")
	if (is.null(max.value)) {
		max.value <- 100
		percent <- TRUE
	} else percent <- FALSE
	if (!is.numeric(max.value))
		stop("'max.value' must be numeric or NULL!")
	## If value is higher than max.value, we erase the message
	erase.only <- (value > max.value)
	## Get the saved data associated with this function
	CmdProgress <- getTemp(".progress", default = list(), mode = "list")

	if (console & progress.bar) {
		## The progress bar consists in two lines:
		## first is a "scale" (only drawn when init == TRUE),
		## second is filled with char according to the actual progression
		if (erase.only) {
			cat ("\n")
			CmdProgress$pos <- NULL
			CmdProgress$scale <- NULL
			assignTemp(".progress", CmdProgress)
		} else {
			if (init || is.null(CmdProgress$pos)) {
				msg1 <- gettext("Progress:")
				l1 <- nchar(msg1)
				if (percent) {
					scalebar <- " 0%---------25%---------50%---------75%--------100%\n"
					CmdProgress$Scale <- 2
				} else {
					## Calculate best scale
					w <- getOption("width") - l1
					CmdProgress$Scale <- (max.value %/% w) + 1
					sl <- round(max.value / CmdProgress$Scale)
					maxV <- as.character(round(max.value))
					ticks <- sl - 1 - nchar(maxV)
					if (ticks > 0) {
						scaleticks <- paste(rep("-", ticks), collapse = "")
					} else scaleticks <- "-"
					scalebar <- paste(" 0", scaleticks, maxV, "\n", sep = "")
				}
				cat(rep(" ", l1), scalebar, msg1, " ", sep = "")
				pos1 <- 0
				CmdProgress$pos <- 0
				assignTemp(".progress", CmdProgress)
			} else pos1 <- CmdProgress$pos
			pos2 <- round(value / CmdProgress$Scale)
			if (pos2 > pos1) {
				CmdProgress$pos <- pos2
				assignTemp(".progress", CmdProgress)
				cat(rep(as.character(char[1]), pos2 - pos1))
			}
		}
		## Under Windows or MacOS, make sure the message is actualized
		flush.console()
	} else if (console & !progress.bar) {
		## A progress indicator in the R console
		## We work only with integer part of the values
		## and transform them into strings of same length
		max.value <- as.character(round(max.value))
		l <- nchar(max.value)
		value <- formatC(round(value), width = l)
		msg1 <- gettext("Progress:")
		l1 <- nchar(msg1)
		msg2 <- gettext("on")
		l2 <- nchar(msg2)
		l3 <- def(CmdProgress$msglength, 0, mode = "numeric", length.out = 1)
		if (l3 < 0) l3 <- 0
		CmdProgress$msglength <- NULL  # Avoid using twice same data
		backspaces <- paste(rep("\b", l3), collapse = "")
		if (erase.only) {
			message <- ""
			cat(backspaces, rep(" ", l3), sep = "")
		} else {
			## Treatment is different if it is 'x%' or 'x on y' display type
			if (percent) {
				message <- paste(msg1, " ", value, "%  ", sep = "",
					collapse = "")
			} else {
				message <- paste(msg1, " ", value, " ", msg2, " ",
					max.value, "  ", sep = "", collapse = "")
			}
		}
		cat(backspaces, message, sep = "")
		CmdProgress$msglength <- nchar(message)
		assignTemp(".progress", CmdProgress)
		## Under Windows or MacOS, make sure the message is actualized
		flush.console()
	}

	## An additional, graphical display of progression may be implemented too
	## using custom functions as items in .progress in SciViews:TempEnv...
	## Here we look for and trigger them...
	if (gui && length(CmdProgress) > 1) {
		## Execute each item of the list that is a function
		for (i in 1:length(CmdProgress))
			if (mode(CmdProgress[[i]]) == "function")
				CmdProgress[[i]](value, max.value)
	}
	invisible(NULL)
}
