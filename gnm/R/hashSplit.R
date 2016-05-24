#  Copyright (C) 2006 Heather Turner
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

hashSplit <- function(string){
    ## An adaptation of some Python code by 'tim'
    ## http://forum.textdrive.com/viewtopic.php?id=3095
    if (!length(string) || !nchar(string))
        return(string)
    s <- strsplit(string, "")[[1]]
    a <- 0
    ans <- vector("list", length(s))
    iq  <- FALSE
    for (z in seq(s)) {
        if (s[z] ==  "#" & !iq) {
            ans[z] <- paste(s[a:(z - 1)], collapse = "")
            a <- z + 1
        }
        else if (s[z] == "\""){
            iq <- !iq
        }
    }
    ans[z] <- paste(s[a:z], collapse = "")
    unlist(ans)
}
