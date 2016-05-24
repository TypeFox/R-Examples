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

# Inserts or updates an SQL query stored in an Open Document Database file
"odb.queries<-" = function(
		odb,
		queryNames,
		value
		)
	{
	# Class checks
	if(!is(odb, "ODB")) {
		stop("'odb' must be an 'ODB' object")
	}
	validObject(odb)
	
	# Args check
	if(length(queryNames) != length(value)) {
		stop("'queryNames' and value must have same lengths")
	}
	if(any(grepl("/", queryNames))) {
		stop("'queryNames' can not contain slashes")
	}
	
	# Encoding
	queryNames <- entityEncode(queryNames)
	value <- entityEncode(value)
	
	# Text reading
	content <- scan(file=sprintf("%s/content.xml", odb@directory), what="", sep="\n", quiet=TRUE, fileEncoding="UTF-8")
	content <- paste(content, collapse="\n")
	
	for(i in 1:length(value)) {
		# Look for the query
		queryTag <- gregexpr(sprintf("<db:query[^>]+db:name=\"%s\"[^>]*/>", queryNames[i]), content)[[1]]
		
		if(queryTag != -1) {
			# Extract
			tagStart <- queryTag[1]
			tagLength <- attr(queryTag, "match.length")[1]
			tag <- substr(content, tagStart, tagStart + tagLength - 1L)
			
			# Replace
			tag <- sub("(?<=db:command=\").*?(?=\")", value[i], tag, perl=TRUE)
			
			# Write
			content <- paste(
				substr(content, 1L, tagStart - 1L),
				tag,
				substr(content, tagStart + tagLength, nchar(content)),
				sep = ""
			)
		} else {
			# Create query section
			if(!grepl("<db:queries>", content)) {
				if(!grepl("</db:data-source>", content)) stop("\"</db:data-source>\" is not supposed to be missing in content.xml")
				content <- sub("(</db:data-source>)", "\\1<db:queries></db:queries>", content)
			}
			
			# Add a query
			content <- sub(
				pattern = "(</db:queries>)",
				replacement = sprintf("<db:query db:name=\"%s\" db:command=\"%s\" db:escape-processing=\"false\"/>\\1", queryNames[i], value[i]),
				x = content
			)			
		}		
	}
	
	# Save XML to file
	conn <- file(sprintf("%s/content.xml", odb@directory), open="wt", blocking=FALSE, encoding="UTF-8")
	cat(content, file=conn)
	close(conn)
	
	return(odb)
}
