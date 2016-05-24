### Copyright (C) 2012 Sylvain Mareschal <maressyl@gmail.com>
### 
### This program is free software: you can redistribute it and/or modify
### it under the terms of the GNU General Public License as published by
### the Free Software Foundation, either version 3 of the License, or
### (at your option) any later version.
### 
### This program is distributed in the hope that it will be useful,
### but WITHOUT ANY WARRANTY; without even the implied warranty of
### MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
### GNU General Public License for more details.
### 
### You should have received a copy of the GNU General Public License
### along with this program.  If not, see <http://www.gnu.org/licenses/>.

# Progression class

setClass(
	Class = "progress",
	representation = representation(
		main = "character",
		iMax = "integer",
		iCurrent = "integer"
	)
)

setClass(
	Class = "progress.console",
	representation = representation(
		"progress",
		pTimes = "numeric",
		eraseLength = "integer"
	)
)

setClass(
	Class = "progress.file",
	representation = representation(
		"progress",
		steps = "integer"
	)
)

setMethod(
	f = "initialize",
	signature = signature(.Object="progress.file"),
	definition = function(.Object, main, iMax, iCurrent=0, nSteps=10) {
		# Arguments
		if (missing(main)) {
			stop(call.=FALSE, "'iCurrent' argument needed")
		}
		if (missing(iMax)) {
			stop(call.=FALSE, "'iMax' argument needed")
		}
		
		# Steps
		steps = pretty(
			x = c(0, iMax),
			n = nSteps
		)
		steps = unique(c(steps, iMax))
		steps = as.integer(steps)
		
		# Filling
		.Object@main = main
		.Object@iMax = as.integer(iMax)
		.Object@iCurrent = as.integer(iCurrent)
		.Object@steps = as.integer(steps)
		
		# Title
		message(
			sprintf("%s : ", .Object@main),
			appendLF = TRUE
		)
		
		return(.Object)
	}
)

setMethod(
	f = "initialize",
	signature = signature(.Object="progress.console"),
	definition = function(.Object, main, iMax, iCurrent=0) {
		# Arguments
		if (missing(main)) {
			stop(call.=FALSE, "'iCurrent' argument needed")
		}
		if (missing(iMax)) {
			stop(call.=FALSE, "'iMax' argument needed")
		}
		
		# Filling
		.Object@main = main
		.Object@iMax = as.integer(iMax)
		.Object@iCurrent = as.integer(iCurrent)
		.Object@pTimes = as.double(proc.time()['elapsed'])
		.Object@eraseLength = as.integer(0)
		
		# Title
		message(
			sprintf("%s : ", .Object@main),
			appendLF = FALSE
		)
		
		return(.Object)
	}
)

setGeneric(
	name = "set",
	def = function(progress, iCurrent) {
		standardGeneric("set")
	}
)

setMethod(
	f = "set",
	signature = signature(progress="progress.file"),
	definition = function(progress, iCurrent) {
		# New i
		progress@iCurrent = iCurrent
		
		if (iCurrent %in% progress@steps) {
			# Message
			message(
				sprintf(
					paste("%0", nchar(progress@iMax), "i/%i", sep=""),
					progress@iCurrent,
					progress@iMax
				),
				appendLF = TRUE
			)
			
			# Updates console
			flush.console()
		}
		
		return(progress)
	}
)

setMethod(
	f = "set",
	signature = signature(progress="progress.console"),
	definition = function(progress, iCurrent) {
		# New i
		progress@iCurrent = iCurrent
		
		# New step time
		progress@pTimes = tail(
			x = c(progress@pTimes, as.double(proc.time()['elapsed'])),
			n = 20
		)
		
		# Erase previous output
		message(
			rep("\b", progress@eraseLength),
			appendLF = FALSE
		)
		
		# New output content
		mess = sprintf(
			"%i/%i (ETA: %s)",
			progress@iCurrent,
			progress@iMax,
			difftimeFmt(
				mean(diff(progress@pTimes)) * (progress@iMax - progress@iCurrent)
			)
		)
		progress@eraseLength = nchar(mess)
		
		# Output
		message(
			mess,
			appendLF = (progress@iCurrent == progress@iMax)
		)
		
		# Updates console
		flush.console()
		
		return(progress)
	}
)
