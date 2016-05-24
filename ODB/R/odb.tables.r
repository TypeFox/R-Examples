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

# Returns the list of tables, fields and comments from an odb connection
odb.tables = function(
		odb
		)
	{
	# Class checks
	if (!is(odb, "ODB")) {
		stop("'odb' must be an 'ODB' object")
	}
	validObject(odb)
	
	# Lists tables in DB
	tryCatch(
		tables <- dbListTables(conn=odb),
		error = function(e) {
			stop(call.=FALSE, "Error while listing tables : \"", conditionMessage(e), "\"")
		},
		warning = function(w) {
			stop(call.=FALSE, "Warning while listing tables : \"", conditionMessage(w), "\"")
		}
	)
	
	# Excludes SYSTEM tables
	tables = tables[ !grepl("^SYSTEM_[A-Z_]+$", tables) ]
	
	# Lists fields of each table
	overview = list()
	for(name in tables) {
		# Execution
		tryCatch(
			query <- dbSendQuery(conn=odb, statement=paste("SELECT * FROM \"", name, "\" WHERE FALSE", sep="")),
			error = function(e) {
				stop(call.=FALSE, "Error while querying table '", name, "' : \"", conditionMessage(e), "\"")
			},
			warning = function(w) {
				stop(call.=FALSE, "Warning while querying table '", name, "' : \"", conditionMessage(w), "\"")
			}
		)
		columns = dbColumnInfo(res=query)
		
		# Factors to Characters
		for(k in 1:ncol(columns)) {
			if (is.factor(columns[,k])) {
				columns[,k] = as.character(columns[,k])
			}
		}
		
		# Adds to output
		overview[[name]] = columns
	}
	
	# Comments
	comments = odb.comments(odb)
	for(tableName in names(overview)) {
		overview[[tableName]][,"comment"] = NA
		for(columnName in names(comments[[ tableName ]])) {
			overview[[ tableName ]][ overview[[ tableName ]]$field.name == columnName , "comment" ] = comments[[ tableName ]][ columnName ]
		}
	}
	
	return(overview)
}
