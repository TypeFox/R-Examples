## Copyright (c) 2004-2012, Ph. Grosjean <phgrosjean@sciviews.org>
##
## This file is part of ZooImage
## 
## ZooImage is free software: you can redistribute it and/or modify
## it under the terms of the GNU General Public License as published by
## the Free Software Foundation, either version 2 of the License, or
## (at your option) any later version.
## 
## ZooImage is distributed in the hope that it will be useful,
## but WITHOUT ANY WARRANTY; without even the implied warranty of
## MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
## GNU General Public License for more details.
## 
## You should have received a copy of the GNU General Public License
## along with ZooImage.  If not, see <http://www.gnu.org/licenses/>.

## Check that the file is a zic file
zicCheck <- function (zicfile)
{	
	zicfile <- as.character(zicfile)
	if (!length(zicfile)) {
		warning("No zicfile provided")
		return(FALSE)
	}
	if (length(zicfile) > 1) {
		warning("testing only first zicfile")
		zicfile <- zicfile[1]
	}
	
	## This should be a .zic file directly
	if (!checkFileExists(zicfile)) return(FALSE)
	
	## First line of the file must be "ZI1", "ZI2", or "ZI3"
	if (!checkFirstLine(zicfile)) return(FALSE) 
	
	## Second line must be [path]
	Line2 <- scan(zicfile, character(), skip = 1, nmax = 2, quiet = TRUE)
	if (tolower(Line2[1]) != "[path]") {
		warning("not a ZIC file, or corrupted at line #2!")
		return(FALSE)
	}
	if (length(Line2) < 2) {
		warning("empty ZIC file is not correct")
		return(FALSE)
	} else return(TRUE)
}
