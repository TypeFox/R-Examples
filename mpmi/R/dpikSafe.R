#  
#     Copyright 2013 Chris Pardy <cpardy@unsw.edu.au>
# 
#     This file is part of the mpmi R package.
# 
#     This program is free software: you can redistribute it and/or modify
#     it under the terms of the GNU General Public License as published by
#     the Free Software Foundation, version 3.
# 
#     This program is distributed in the hope that it will be useful,
#     but WITHOUT ANY WARRANTY; without even the implied warranty of
#     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#     GNU General Public License for more details.
# 
#     You should have received a copy of the GNU General Public License
#     along with this program.  If not, see <http://www.gnu.org/licenses/>.
#  


# Function to handle cases where dpik() is unable to estimate
# a bandwidth, usually because a data vector has an interquartile range of 0.
dpikSafe <- function(x, ...)
{
	result <- try(dpik(x, ...), silent = TRUE)

	if (class(result) == "try-error")
	{
		msg <- geterrmessage()
		if (grepl("scale estimate is zero for input data", msg))
		{
			warning("Using standard deviation as scale estimate, probably because IQR == 0")
			result <- dpik(x, scalest = "stdev", ...)	
		} else 
		{
			stop(msg)
		}		
		
	}
	return(result)
}
