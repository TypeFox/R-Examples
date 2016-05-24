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

# Converts XML entities
entityEncode <- function(from) {
	from <- gsub("&", "&amp;", from)
	from <- gsub("'", "&apos;", from)
	from <- gsub(">", "&gt;", from)
	from <- gsub("<", "&lt;", from)
	from <- gsub("\"", "&quot;", from)
	from <- gsub("\n", "&#x0a;", from)
	from <- gsub("\t", "&#x09;", from)
	return(from)
}
entityDecode <- function(from) {
	from <- gsub("&amp;", "&", from)
	from <- gsub("&apos;", "'", from)
	from <- gsub("&gt;", ">", from)
	from <- gsub("&lt;", "<", from)
	from <- gsub("&quot;", "\"", from)
	from <- gsub("&#x0a;", "\n", from)
	from <- gsub("&#x09;", "\t", from)
	return(from)
}
