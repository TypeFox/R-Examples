#  Copyright (C) 2005 David Firth and Heather Turner
#
#  This program is free software; you can redistribute it and/or modify
#  it under the terms of the GNU General Public License as published by
#  the Free Software Foundation; either version 2 or 3 of the License
#  (at your option).
#
#  This program is distributed in the hope that it will be useful,
#  but WITHOUT ANY WARRANTY; without even the implied warranty of
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#  GNU General Public License for more details.
#
#  A copy of the GNU General Public License is available at
#  http://www.r-project.org/Licenses/

unlistOneLevel <- function(theList){
     result <- vector(length = sum(sapply(theList,
                      function(x) if(is.list(x)) length(x) else 1)),
                      mode = "list")
     count <- 0
     for (i in seq(theList)){
         theItem <- theList[[i]]
         if (is.list(theItem)){
             for (j in seq(theItem)){
                 count <- count + 1
                 result[[count]] <- theItem[[j]]
             }
         }
         else {
             count <- count + 1
             result[[count]] <- theItem
         }
     }
     return(result[1:count])
}
