### tclTask.R - Functions to schedule task to be executed later using Tcl after
### Copyright (c), Philippe Grosjean (phgrosjean@sciviews.org)
### Licensed under LGPL 3 or above
###
### Changes:
### - 2009-07-02: fisrt version (for tcltk2_1.1-0)

tclAfter <- function (wait, fun)
{
	## This is the basic Tcl command, do prefer tclTaskSchedule()!
	wait <- as.integer(wait)[1]
	if (wait <= 0) wait <- "idle"  # Schedule task on next event loop
	## Check fun
	if (!is.function(fun))
		stop("'fun' must be a function")
	## Install a new Tcl timer
	tcl("after", wait, fun)
}

tclAfterCancel <- function (task)
{
	## Cancel a Tcl timer (no effect if the timer does not exist)
	tcl("after", "cancel", as.character(task)[1])
}

tclAfterInfo <- function (task = NULL)
{
	## Get info about a Tcl timer, or list all current ones (using task = NULL)
	if (is.null(task)) {
		return(tcl("after", "info"))
	} else {
		## First check that task exists
		task <- as.character(task)[1]
		ok <- tclvalue(.Tcl(paste("catch {after info ", task, "}", sep = "")))
		if (ok == 0) {
			return(tcl("after", "info", task))
		} else return(NULL)
	}
}

print.tclTask <- function (x, ...)
{
	## Look when the task is run
	if (x$wait == "idle") {
		cat("tclTask '", x$id, "' scheduled on next event loop\n", sep = "")
	} else {
		cat("tclTask '", x$id, "' scheduled after ", x$wait, " ms ", sep = "")
		## Determine how much is remaining
		rem <- x$started + x$wait - proc.time()["elapsed"] * 1000
		if (rem <= 0) {
			cat("(elapsed)\n")
		} else {
			cat("(", as.integer(rem), " remaining)\n", sep = "")
		}
	}
	## Look if it is rescheduled
	if (isTRUE(x$redo)) {
		cat("Rescheduled forever\n")
	} else if (x$redo == FALSE || x$redo <= 0) {
		cat("Not rescheduled\n")
	} else if (x$redo == 1) {
		cat("Rescheduled once\n")
	} else {
		cat("Rescheduled", x$redo, "times\n")
	}
	## Print the command to be executed
	cat("runs:\n")
	print(x$expr)
	return(invisible(x))
}

tclTaskSchedule <- function (wait, expr, id = "task#", redo = FALSE)
{
	## Schedule a task to be executed after 'wait' ms
	## If 'wait' is <= 0, schedule for execution on the next event loop
	## Id is the task name to use (if the task already exists, it is deleted
	## and replaced by the new definition)

	wait <- as.integer(wait)[1]
	if (wait <= 0) wait <- "idle"  # Schedule task on next event loop

	Tasks <- .getTclTasks()
	TNames <- ls(Tasks, all.names = TRUE)

	id <- as.character(id)[1]
	## If 'id' contains '#', replace it by a number (first one available)
	## but don't allow more than 1000 tasks with same name (to avoid bad
	## situations with buggy code like infinite loops or so)
	if (grepl("#", id)) {
		for (i in 1:1000) {
			Id <- sub("#", i, id)
			if (!Id %in% TNames) break
		}
		if (Id %in% TNames)
			stop("Too many tclTasks!")
	} else {
		## Delete the task if it already exists
		if (id %in% TNames) tclTaskDelete(id)
		Id <- id
	}

	if (!isTRUE(redo)) {
		redo <- as.integer(redo)[1]
		if (redo <= 0) redo <- FALSE
	}

	## Schedule the task, but don't run expr directly, but through tclTaskRun()
	## Note: if I use tcl("after", wait, tclTaskRun(Id), R is blocked until the
	## task is done. Here, I must provide the name of a function without args)
	task <- .makeTclTask(id = Id, wait = wait)

	## Create a tclTask object containing all info about this task
	res <- list(task = task, id = Id, expr = substitute(expr),
		started = proc.time()["elapsed"] * 1000, wait = wait,
		redo = redo)
	class(res) <- c("tclTask", class(res))

	## Add this task to the list
	Tasks[[Id]] <- res

	return(invisible(res))
}

tclTaskRun <- function(id)
{
	## Execute the code associated with a given task and detemine if the task
	## should be rescheduled again (repeat argument)
	id <- as.character(id)[1]

	Tasks <- .getTclTasks()
	Task <- Tasks[[id]]
	if (is.null(Task)) {
		warning("tclTask '", id, "' is not found")
		return(invisible(FALSE))
	}
	## Make sure to indicate that we run it once
	if (!is.logical(Task$redo)) {
		Task$redo <- Task$redo - 1
		if (Task$redo < 1) Task$redo <- FALSE
	}
	## Update the original object too
	Tasks[[id]] <- Task

	## Run the code associate with this task
	eval(Task$expr, envir = .GlobalEnv)

	## Should we delete this task (if repeat is FALSE), or reschedule it?
	## Note, we read Task again, in case fun() would have changed something there!
	Task <- Tasks[[id]]
	## Make sure the tcl timer is destroyed (in case tclTaskRun() is
	## triggered otherwise)
	tclTaskDelete(id)
	if (Task$redo) {
		## Reschedule the task
		Task$task <- .makeTclTask(id = id, wait = Task$wait)
		## and update information in .tclTasks
		Tasks[[id]] <- Task
	}
	return(invisible(TRUE))
}

tclTaskGet <- function(id = NULL, all = FALSE)
{
	## If id is NULL, list all scheduled tasks, otherwise, give info about a
	## particular scheduled task
	if (is.null(id)) {
		return(ls(.getTclTasks(), all.names = all))
	} else {
		## Get the data associated with a scheduled task
		return(.getTclTasks()[[id]])
	}
}

tclTaskChange <- function (id, expr, wait, redo)
{
	## Change a characteristic of a scheduled task
	## Is there something to change?
	if (missing(expr) && missing(wait) && missing(redo))
		return(invisible(FALSE))
	## Get task and change it
	Tasks <- .getTclTasks()
	Task <- Tasks[[id]]
	if (is.null(Task)) {
		warning("tclTask '", id, "' is not found")
		return(invisible(FALSE))
	}
	if (!missing(expr)) Task$expr <- substitute(expr)
	if (!missing(wait )) {
		wait <- as.integer(wait)[1]
		if (wait <= 0) wait <- "idle"  # Schedule task on next event loop
		Task$wait <- wait
	}
	if (!missing(redo)) {
		if (!isTRUE(redo)) {
			redo <- as.integer(redo)[1]
			if (redo <= 0) redo <- FALSE
		}
		Task$redo <- redo
	}
	## Delete the task and recreate it with the new parameters
	tclTaskDelete(id)
	Task$task <- .makeTclTask(id = id, wait = Task$wait)

	## Update Tasks
	Tasks[[id]] <- Task

	return(invisible(TRUE))
}

tclTaskDelete <- function (id)
{
	Tasks <- .getTclTasks()
	## Remove a previously scheduled task (if id s NULL, then, remove all tasks)
	if (is.null(id)) {
		## Delete all current tasks
		for (Task in ls(Tasks, all.names = TRUE))
			tclAfterCancel(Tasks[[Task]]$task)
		## Eliminate .tclTasks environment from Sciviews:TempEnv
		rm(list = ".tclTasks", envir = .TempEnv())
	} else {
		## Delete only one task
		Task <- Tasks[[id]]
		if (!is.null(Task)) {  # The task exists
			tclAfterCancel(Task$task)
			rm(list = id, envir = Tasks)
		}
	}
}

.getTclTasks <- function ()
{
	## Retrieve references to all scheduled tasks
	res <- .getTemp(".tclTasks", default = NULL)
	if (is.null(res)) {
		res <- new.env(parent = .TempEnv())
		.assignTemp(".tclTasks", res)
	}
	return(res)
}

.makeTclTask <- function (id, wait)
{
	run <- function ()
		eval(parse(text = paste('tclTaskRun("', id, '")', sep = "")))
	task <- tclAfter(wait, run)
	return(task)
}
