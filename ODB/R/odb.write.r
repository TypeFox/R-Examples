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

# Executes a series of SQL queries through an ODB connection (Write Only)
odb.write = function(
		odb,
		sqlQueries,
		onError = c("warning", "stop"),
		progress = c("console", "file", "none")
		)
	{
	# Class checks
	if (!is(odb, "ODB")) {
		stop("'odb' must be an 'ODB' object")
	}
	validObject(odb)
	
	# Checks
	if (!is.character(sqlQueries) || length(sqlQueries) == 0 || any(is.na(sqlQueries))) {
		stop(call.=FALSE, "'sqlQueries' must be a non NA character vector")
	}	
	
	# Args matching
	onError = match.arg(onError)
	
	# Block on error or continue
	if (onError == "stop") {
		errorFun = function(e) {
			stop(call.=FALSE, "Query #", i, " : \"", conditionMessage(e), "\"")
		}
	} else {
		errorFun = function(e) {
			warning(call.=FALSE, immediate.=TRUE, "Query #", i, " : \"", conditionMessage(e), "\"")
		}
	}
	
	# Initiating progression
	progress = match.arg(progress)
	if (progress != "none") {
		progressObj = new(paste("progress", progress, sep="."), main="Execution", iMax=length(sqlQueries))
	}
	
	# Queries execution
	for(i in 1:length(sqlQueries)) {
		tryCatch(
			dbSendUpdate(odb, sqlQueries[i]),
			error = errorFun
		)
		
		# Progression
		if (progress != "none") {
			progressObj = set(progressObj, i)
		}
	}
}
