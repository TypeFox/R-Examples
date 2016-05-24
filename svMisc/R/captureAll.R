## Inspired by 'capture.output' and the old .try_silent in utils package
## Requires: R >= 2.13.0 [??]
captureAll <- function (expr, split = TRUE, echo = TRUE, file = NULL,
markStdErr = FALSE)
{
	if (!is.expression(expr))
		if (is.na(expr)) return(NA) else
		stop("expr must be an expression or NA")
	
	## TODO: support for 'file'
	## markStdErr: if TRUE, stderr is separated from sddout by STX/ETX character

	last.warning <- list()
	Traceback <- list()
	NframeOffset <- sys.nframe() + 23L	# frame of reference (used in traceback)
										# + length of the call stack when a
										# condition occurs
	## Note: if 'expr' is a call not expression, 'NframeOffset' is lower by 2
	## (i.e. 24): -1 for lapply, -1 for unwrapping 'expression()'

	getWarnLev <- function() options('warn')[[1L]]	# This may change in course
													# of evaluation, so must be
													# retrieved dynamically
	rval <- NULL
	tconn <- textConnection("rval", "w", local = TRUE)
	split <- isTRUE(split)
	if (split) {
		## This is required to print error messages when we are, say, in a
		## browser() environment
		sink(stdout(), type = "message")
	} else {
		## This is the conventional way to do it
		sink(tconn, type = "message")
	}
	sink(tconn, type = "output", split = split)
	#sink(tconn, type = "message")
	on.exit({
		sink(type = "message")
		sink(type = "output")
		close(tconn)
	})

	inStdOut <- TRUE

	if (isTRUE(markStdErr)) {
		putMark <- function (toStdout, id) {
			if (inStdOut) {
				if (!toStdout) {
					cat("\x03")
					inStdOut <<- FALSE
				}
			} else { # in StdErr stream
				if (toStdout) {
					cat("\x02")
					inStdOut <<- TRUE
				}
			}
			#cat("<", id, inStdOut, ">")
		}
	} else {
		putMark <- function (toStdout, id) {}
	}

	evalVis <- function (x) {
		## Do we print the command? (note that it is reformatted here)
		if (isTRUE(echo)) {
			## Reformat long commands... and possibly abbreviate them
			cmd <- deparse(x)
			l <- length(cmd)
			if (l > 7) cmd <- c(cmd[1:3], "[...]", cmd[(l-2):l])
			cat(":> ", paste(cmd, collapse = "\n:+ "), "\n", sep = "")
		}
		res <- withVisible(eval(x, .GlobalEnv))
		## Do we have result to print?
		if (inherits(res, "list") && res$visible) {
			## Printing is veeery slow on Windows when split = TRUE
			## => unsplit temporarily, and print twice instead!
			#print(res$value)
			
			if (split) {
				sink(type = "message")
				sink(type = "output")
				## Print first to the console
				try(print(res$value), silent = TRUE)
				sink(tconn, type = "message")
				sink(tconn, type = "output", split = FALSE)
				## Print a second time to the connection
				try(print(res$value), silent = TRUE)
				## Resink with split = TRUE
				sink(type = "message")
				sink(type = "output")
				sink(stdout(), type = "message")
				sink(tconn, type = "output", split = TRUE)
			} else {
				## This is the conventional way to do it
				print(res$value)	
			}	
		}
		
		return(res)
	}

	formatMsg <- function (msg) {
		## For some reasons, 'Error: ' and 'Error in ' are not translated,
		## although the rest of the message is correctly translated
		## This is a workaround for this little problem
		res <- .makeMessage(msg)
		res <- sub("^Error: ", ngettext(1, "Error: ", "Error: ", domain = "R"),
			res)
		res <- sub("^Error in ", ngettext(1, "Error in ", "Error in ",
			domain = "R"), res)
		return(res)
	}

	restartError <- function (e, calls) {
		## Remove call (eval(expr, envir, enclos)) from the message
		ncls <- length(calls)

		##DEBUG: cat("n calls: ", ncls, "NframeOffset: ", NframeOffset, "\n")
		if (isTRUE(all.equal(calls[[NframeOffset]], e$call,
			check.attributes = FALSE)))
			e$call <- NULL

		Traceback <<- rev(calls[-c(seq.int(NframeOffset), (ncls - 1L):ncls)])

#> cat(captureAll(expression(1:10, log(-1),log(""),1:10)), sep="\n")
#Error in calls[[NframeOffset]]: subscript out of bounds
#Warning message:
#In log(-1) : NaNs produced

		putMark(FALSE, 1)
		cat(formatMsg(e))
		if (getWarnLev() == 0L && length(last.warning) > 0L)
			cat(ngettext(1, "In addition: ", "In addition: ", domain = "R"))
	}

	res <- tryCatch(withRestarts(withCallingHandlers({
			## TODO: allow for multiple expressions and calls (like in
			## 'capture.output'). The problem here is how to tell 'expression'
			## from 'call' without evaluating it?
			##list(evalVis(expr))
			lapply(expr, evalVis)
		},

		error = function (e) invokeRestart("grmbl", e, sys.calls()),
		warning = function (e) {
			## Remove call (eval(expr, envir, enclos)) from the message
			if (isTRUE(all.equal(sys.call(NframeOffset), e$call,
				check.attributes = FALSE)))
				e$call <- NULL

			last.warning <<- c(last.warning, structure(list(e$call),
				names = e$message))

			if (getWarnLev() != 0L) {
				putMark(FALSE, 2)
				## Was:
				#.Internal(.signalCondition(e, conditionMessage(e), conditionCall(e)))
				#.Internal(.dfltWarn(conditionMessage(e), conditionCall(e)))
				.signalSimpleWarning(conditionMessage(e), conditionCall(e))
				putMark(TRUE, 3)
			}
			invokeRestart("muffleWarning")
		}),
		## Restarts:

		## Handling user interrupts. Currently it works only from within R.
		##TODO: how to trigger interrupt via socket connection?
		abort = function (...) {
			putMark(FALSE, 4)
			cat("<aborted!>\n") #DEBUG
		},

		interrupt = function (...) cat("<interrupted!>\n"), #DEBUG: this does not seem to be ever called.

		muffleWarning = function () NULL,
		grmbl = restartError),
		error = function (e) { ##XXX: this is called if warnLevel == 2
			putMark(FALSE, 5)
			cat(formatMsg(e))
			e #identity
		},
		finally = {}
	)

	if (getWarnLev() == 0L) {
		nwarn <- length(last.warning)
		assign("last.warning", last.warning, envir = baseenv())

		if (nwarn > 0L) putMark(FALSE, 6)
		if (nwarn <= 10L) {
			print.warnings(last.warning)
		} else if (nwarn < 50L) {
			## This is buggy and does not retrieve a translation of the message!
			#cat(gettextf("There were %d warnings (use warnings() to see them)\n",
			#	nwarn, domain = "R"))
			msg <- ngettext(1,
				"There were %d warnings (use warnings() to see them)\n",
				"There were %d warnings (use warnings() to see them)\n",
				domain = "R")
			cat(sprintf(msg, nwarn))		      
		} else {
			cat(ngettext(1,
				"There were 50 or more warnings (use warnings() to see the first 50)\n",
				"There were 50 or more warnings (use warnings() to see the first 50)\n",
				domain = "R"))
		}
	}

	putMark(TRUE, 7)

	sink(type = "message")
	sink(type = "output")
	close(tconn)
	on.exit()

	## Allow for tracebacks of this call stack:
	assign(".Traceback", lapply(Traceback, deparse), envir = baseenv())

	## Make sure last line ends up with \n
	l <- length(rval)
	if (l) rval[l] <- paste(rval[l], "\n", sep = "")
	rval
}
