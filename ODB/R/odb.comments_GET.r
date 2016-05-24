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

# Returns a list (by table) of character vector (by column) of all column comments in an Open Document Database file
odb.comments = function(
		odb,
		tableNames = NULL,
		columnNames = NULL,
		simplify = TRUE
		)
	{
	# Class checks
	if(!is(odb, "ODB")) {
		stop("'odb' must be an 'ODB' object")
	}
	validObject(odb)
	
	# Text reading
	content <- scan(file=sprintf("%s/content.xml", odb@directory), what="", sep="\n", quiet=TRUE, fileEncoding="UTF-8")
	content <- paste(content, collapse="\n")
	
	# Target positions
	tableTags <- gregexpr("<db:table-representation[^s]", content)[[1]]
	columnTags <- gregexpr("<db:column[^s]", content)[[1]]
	nameParams <- gregexpr("name=\"(.*?)\"", content)[[1]]
	helpParams <- gregexpr("help-message=\"(.*?)\"", content)[[1]]
	
	# Table / Column linkage
	columnTable <- rep(as.integer(NA), length(columnTags))
	for(i in 1:length(tableTags)) {
		columnTable[ columnTags > tableTags[i] ] <- i
	}
	
	# Extraction
	results <- list()
	for(i in unique(columnTable)) {
		# Table name
		tName <- extractParam(content, nameParams, tableTags[i], "name=\"(.*?)\"")
		if(is.null(tableNames) || tName %in% tableNames) {
			# For each column in the table
			for(j in which(columnTable == i)) {
				# Column name
				cName <- extractParam(content, nameParams, columnTags[j], "name=\"(.*?)\"")
				if(is.null(columnNames) || cName %in% columnNames) {
					# Comment extraction
					comm <- extractParam(content, helpParams, columnTags[j], "help-message=\"(.*?)\"")
					results[[ tName ]][ cName ] <- comm
				}
			}
		}
	}
	
	# Simplification
	if(simplify) {
		if(length(tableNames) == 1) {
			# Single table queried
			if(length(results) == 1) {
				# Exists
				results <- results[[1]]
				if(length(columnNames) == 1 && length(results) == 1) {
					names(results) <- NULL
				}
			} else {
				# Does not exist
				results <- list()
				results[[ tableNames ]] <- character(0)
			}			
		}
	}
	
	return(results)
}

extractParam <- function(xmlText, paramMatches, position, pattern) {
	m <- which(paramMatches > position)[1]
	m <- substr(xmlText, paramMatches[m], paramMatches[m] + attr(paramMatches, "match.length")[m] - 1L)
	m <- sub(pattern, "\\1", m)
	m <- entityDecode(m)
	return(m)		
}
