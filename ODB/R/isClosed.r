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

# Check if a DBI or odb connection is closed
isClosed = function(
		connection
		)
	{
	# Connection check
	if (!inherits(connection, "DBIConnection")) {
		stop("'connection' must inherit from 'DBIConnection', or be an 'ODB' connection")
	}
	
	# Test
	closed = tryCatch(
		dbSendUpdate(connection, ""),
		error = function(e) {
			return(grepl("Connection is closed", conditionMessage(e)))
		}
	)
	
	return(!is.null(closed))
}
