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

# Create a new Open Document Database from a template
odb.create = function(
		odbFile,
		template = NULL,
		overwrite = c("warning", "do", "skip", "stop")
		)
	{
	# Args matching
	overwrite = match.arg(overwrite)
	
	# Checks
	if (is.null(template)) {
		template = system.file("tools/template.odb", package="ODB")
	}
	if (!file.exists(template)) {
		stop(call.=FALSE, "'template' file doesn't exist, it must be an .odb file")
	}
	
	# Overwriting
	if (file.exists(odbFile)) {
		if (overwrite == "warning")   {
			over = TRUE
			warning(call.=FALSE, "'odbFile' has been overwritten")
		} else if (overwrite == "stop") {
			over = FALSE
			stop(call.=FALSE, "'odbFile' already exists")
		} else if (overwrite == "skip") {
			over = FALSE
		} else if (overwrite == "do") {
			over = TRUE
		}
	} else {
		over = FALSE
	}
	
	# Copy
	invisible(
		file.copy(from=template, to=odbFile, overwrite=over)
	)
}
