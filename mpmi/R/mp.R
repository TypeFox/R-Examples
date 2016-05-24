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


mp <- function(mat, ...)
{
    image(max(mat, na.rm = TRUE) - t(mat[dim(mat)[1]:1, ]), ...)
}

