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

# Converts a time difference from seconds to multiple units
difftimeFmt = function(
		x
		)
	{
	# difftime objects
	if (class(x) == "difftime") {
		units(x) = "secs"
	}
	
	# Coercion
	x = as.double(x)
	
	# Check
	if (length(x) != 1) {
		stop("'x' must be single value")
	}
	
	# Always positive
	x = abs(x)
	
	if (!is.na(x)) {
		# Units
		mseconds = round((x %% 1) * 1000)
		x = floor(x)
		seconds = x %% 60
		x = x %/% 60
		minutes = x %% 60
		x = x %/% 60
		hours = x %% 24
		days = x %/% 24
		
		# Base format
		fmt = "%02i:%02i:%02i"
		args = list(
			hours,
			minutes,
			seconds
		)
		
		# Days
		if (days > 0) {
			fmt = paste("%id ", fmt, sep="")
			args = c(days, args)
		}
		
		# Milliseconds
		if (mseconds > 0) {
			fmt = paste(fmt, ".%03i", sep="")
			args = c(args, mseconds)
		}
		
		# Format
		args$fmt = fmt
		return(
			do.call(
				what = "sprintf",
				args = args
			)
		)
	} else {
		# NA
		return(NA_character_)
	}
}
