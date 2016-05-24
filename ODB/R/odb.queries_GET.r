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

# Returns a character vector of SQL queries stored in an Open Document Database file
odb.queries = function(
		odb,
		queryNames = NULL
		)
	{
	# Class checks
	if (!is(odb, "ODB")) {
		stop("'odb' must be an 'ODB' object")
	}
	validObject(odb)
	
	# Text reading
	content <- scan(file=sprintf("%s/content.xml", odb@directory), what="", sep="\n", quiet=TRUE, fileEncoding="UTF-8")
	content <- paste(content, collapse="\n")
	
	# Target positions
	queryTags <- gregexpr("<db:query.*?/>", content)[[1]]
	
	# Extraction
	results <- character(0)
	for(i in 1:length(queryTags)) {
		if(queryTags[i] != -1) {
			tag <- substr(content, queryTags[i], queryTags[i] + attr(queryTags, "match.length")[i] - 1L)
			name <- entityDecode(sub("^.*db:name=\"(.*?)\".*$", "\\1", tag))
			if(is.null(queryNames) || name %in% queryNames) {
				command <- entityDecode(sub("^.*db:command=\"(.*?)\".*$", "\\1", tag))
				results[ name ] <- command
			}
		}
	}
	
	return(results)	
}
