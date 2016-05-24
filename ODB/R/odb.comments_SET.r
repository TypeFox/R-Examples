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

# Set column comments from specified tables in an ODB connection
"odb.comments<-" = function(
		odb,
		tableNames,
		columnNames,
		value
		)
	{
	# Class checks
	if (!is(odb, "ODB")) {
		stop("'odb' must be an 'ODB' object")
	}
	validObject(odb)
	
	# Encoding
	tableNames <- entityEncode(tableNames)
	columnNames <- entityEncode(columnNames)
	value <- entityEncode(value)
	
	# Recycling
	maxLength = max(length(tableNames), length(columnNames), length(value))
	columnNames = rep(columnNames, length.out=maxLength)
	tableNames = rep(tableNames, length.out=maxLength)
	value = rep(value, length.out=maxLength)
	
	# Text reading
	content <- scan(file=sprintf("%s/content.xml", odb@directory), what="", sep="\n", quiet=TRUE, fileEncoding="UTF-8")
	content <- paste(content, collapse="\n")
	
	# Create tables section
	if(!grepl("<db:table-representations>", content)) {
		if(!grepl("</office:database>", content)) stop("\"</office:database>\" is not supposed to be missing in content.xml")
		content <- sub("(</office:database>)", "<db:table-representations></db:table-representations>\\1", content)
	}
	
	for(i in 1:maxLength) {
		# Create target table section
		tableTag <- sprintf("<db:table-representation db:name=\"%s\"", tableNames[i])
		if(!grepl(tableTag, content, fixed=TRUE)) {
			content <- sub("(</db:table-representations>)", sprintf("%s/>\\1", tableTag), content)
		}
		
		# Expand target table section
		tableTag <- sprintf("(<db:table-representation db:name=\"%s\")/>", tableNames[i])
		if(grepl(tableTag, content)) {
			content <- sub(tableTag, sprintf("\\1><db:columns></db:columns></db:table-representation>", tableTag), content)
		}
		
		# Table node
		tablePos <- gregexpr(sprintf("<db:table-representation db:name=\"%s\">.*?</db:table-representation>", tableNames[i]), content)[[1]]
		tableTag <- substr(content, tablePos[1], tablePos[1] + attr(tablePos, "match.length")[1] - 1L)
		columnPos <- gregexpr(sprintf("<db:column[^>]+name=\"%s\"[^>]+/>", columnNames[i]), tableTag)[[1]]
		if(columnPos[1] != -1) {
			# Extract column
			tagStart <- tablePos[1] + columnPos[1] - 1L
			tagLength <- attr(columnPos, "match.length")[1]
			tag <- substr(content, tagStart, tagStart + tagLength - 1L)
			
			# Replace
			tag <- sub("(?<=db:help-message=\").*?(?=\")", value[i], tag, perl=TRUE)
			
			# Write
			content <- paste(
				substr(content, 1L, tagStart - 1L),
				tag,
				substr(content, tagStart + tagLength, nchar(content)),
				sep = ""
			)
		} else {
			# Extract whole table
			tagStart <- tablePos[1]
			tagLength <- attr(tablePos, "match.length")[1]
			tag <- substr(content, tagStart, tagStart + tagLength - 1L)
			
			# Replace
			tag <- sub(
				pattern = "(</db:columns>)",
				replacement = sprintf("<db:column db:name=\"%s\" db:help-message=\"%s\"/>\\1", columnNames[i], value[i]),
				x = tag
			)
			
			# Write
			content <- paste(
				substr(content, 1L, tagStart - 1L),
				tag,
				substr(content, tagStart + tagLength, nchar(content)),
				sep = ""
			)
		}
		
	}
	
	# Save XML to file
	conn <- file(sprintf("%s/content.xml", odb@directory), open="wt", blocking=FALSE, encoding="UTF-8")
	cat(content, file=conn)
	close(conn)
	
	return(odb)
}
