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

# Close a connection to an Open Document Database file
odb.close = function(
		odb,
		write = TRUE
		)
	{
	# Class checks
	if(!is(odb, "ODB")) {
		stop("'odb' must be an 'ODB' object")
	}
	validObject(odb)
	
	# Updates the .odb file properly
	if(isTRUE(write)) {
		# Zip check
		cat("test", file="test.txt")
		check <- utils::zip(
			zipfile = "test.zip",
			files = "test.txt",
			flags = "-r9Xq"
		)
		unlink("test.txt")
		if(check == 0) unlink("test.zip")
		else           stop("utils::zip does not seem to work on your platform, ODB file update can not be performed.")
		
		# HSQLDB disconnection with file compaction
		tryCatch(
			dbSendUpdate(odb, "SHUTDOWN COMPACT"),
			error = function(e) {
				stop("HSQL disconnection failed : \"", conditionMessage(e), "\"")
			}
		)
		
		# Renames HSQLDB files back
		dbFiles = c("backup", "data", "properties", "script")
		dbFiles.odb = paste(odb@directory, "/database/", dbFiles, sep="")
		dbFiles.jdbc = paste(odb@directory, "/database/ODB.", dbFiles, sep="")
		for(i in 1:length(dbFiles)) {
			if (file.exists(dbFiles.jdbc[i])) {
				file.rename(dbFiles.jdbc[i], dbFiles.odb[i])
			}
		}
		
		# Updates the ODB archive
		wd <- getwd()
		setwd(odb@directory)
		status <- utils::zip(
			zipfile = odb@odbFile,
			files = c("content.xml", "database"),
			flags = "-r9Xq"
		)
		if(status != 0) stop(sprintf("ODB file building failed (code %s)", status))
		setwd(wd)
	} else {
		# HSQLDB disconnection attempt
		tryCatch(
			dbSendUpdate(odb, "SHUTDOWN"),
			error = function(e) {
				warning("HSQL disconnection failed : \"", conditionMessage(e), "\"")
			}
		)
	}
	
	# DBI disconnection attempt
	error = tryCatch(
		dbDisconnect(odb),
		error = function(e) {
			warning("DBI disconnection failed : \"", conditionMessage(error), "\"")
		}
	)
	
	# Removes temporary files
	unlink(odb@directory, recursive=TRUE)
	
	invisible(TRUE)
}
