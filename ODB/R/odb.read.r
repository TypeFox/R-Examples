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

# Executes a series of SQL queries through an ODB connection (Read Only)
odb.read = function(
		odb,
		sqlQuery,
		stringsAsFactors = FALSE,
		check.names = FALSE,
		encode = TRUE,
		autoLogical = TRUE
		)
	{
	# Class checks
	if (!is(odb, "ODB")) {
		stop("'odb' must be an 'ODB' object")
	}
	validObject(odb)
	
	# Query check
	if (!is.character(sqlQuery) || length(sqlQuery) != 1 || is.na(sqlQuery)) {
		stop(call.=FALSE, "'sqlQuery' must be a unique non NA character vector")
	}	
	
	# Execution
	tryCatch(
		query <- dbSendQuery(odb, sqlQuery),
		error = function(e) {
			stop(call.=FALSE, "Error while executing SQL query  : \"", conditionMessage(e), "\"")
		},
		warning = function(w) {
			stop(call.=FALSE, "Warning while executing SQL query  : \"", conditionMessage(w), "\"")
		}
	)
	
	# Results retrieving
	results = fetch(res=query, n=-1)
	
	## Clearing results (not currently supported by JDBC)
	# dbClearResult(res=query)
	
	# Reverts check.names
	if (!check.names) {
		columns = dbColumnInfo(res=query)
		names(results) = as.character(columns$field.name)
	}
	
	# Factors to Characters
	if (!stringsAsFactors) {
		for(k in 1:ncol(results)) {
			if (is.factor(results[,k])) {
				results[,k] = as.character(results[,k])
			}
		}
	}
	
	# Re-encoding, keeping NAs
	if (encode) {
		for(k in 1:ncol(results)) {
			if (is.character(results[,k]) | is.factor(results[,k])) {
				naVals = is.na(results[,k])
				results[,k] = iconv(results[,k], from="UTF-8", to="")
				if (length(naVals) > 0) {
					results[naVals,k] = NA
				}
			}
		}
	}
	
	# Logical
	if (autoLogical) {
		for(k in 1:ncol(results)) {
			if (all(is.na(results[,k]) | results[,k] %in% c("true", "false"))) {
				results[,k] = as.logical(results[,k])
			}
		}
	}
	
	return(results)
}
