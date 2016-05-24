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

# Generates INSERT INTO queries from a data.frame fitting an ODB base
odb.insert = function(
		odb,
		tableName,
		data,
		execute = TRUE,
		dateFormat = "%Y-%m-%d",
		...
		)
	{
	# Class checks
	if (!is(odb, "ODB")) {
		stop("'odb' must be an 'ODB' object")
	}
	validObject(odb)
	
	# Args checks
	if (!is.character(tableName) || length(tableName) != 1 || is.na(tableName)) {
		stop(call.=FALSE, "'tableName' must be a unique non NA character vector")
	}	
	data = as.data.frame(data)
	
	# Gets table definition
	tryCatch(
		query <- dbSendQuery(
			conn = odb,
			statement = paste("SELECT * FROM", tableName, "WHERE FALSE")
		),
		error = function(e) {
			stop(call.=FALSE, "Error while querying table '", tableName, "' : \"", conditionMessage(e), "\"")
		},
		warning = function(w) {
			stop(call.=FALSE, "Warning while querying table '", tableName, "' : \"", conditionMessage(w), "\"")
		}
	)
	overview = dbColumnInfo(res=query)
	
	# Size check
	if (nrow(overview) != ncol(data)) {
		stop("'tableName' table (", nrow(overview), ") and 'data' (", ncol(data), ") column counts don't match")
	}
	
	# Conversion from data.frame to varchar matrix
	mtx = as.matrix(
		data.frame(
			lapply(data, as.character),
			stringsAsFactors = FALSE
		)
	)
	
	# Cell processing
	isna = is.na(mtx)
	mtx[, overview$field.type == "DATE" ] = as.character(strptime(mtx[, overview$field.type == "DATE" ], dateFormat))
	mtx[, overview$data.type == "character" ] = paste(sep="", "'", gsub("'", "''", mtx[, overview$data.type == "character" ]), "'")
	mtx[, overview$data.type == "numeric" ] = gsub("[^0-9\\.-]", "", mtx[, overview$data.type == "numeric" ])
	mtx[ isna ] = "NULL"
	
	# Aggregation
	SQL = paste(sep="",
		"INSERT INTO ", tableName, " VALUES (",
		apply(mtx, 1, paste, collapse=", "),
		");"
	)
	
	# Execution
	if (execute) {
		odb.write(odb, SQL, ...)
		invisible(SQL)
	} else {
		return(SQL)
	}
}
