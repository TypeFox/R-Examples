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

instances <- function(term, instances = 1){
    term <- match.call()$term
    if (!"inst" %in% names(formals(match.fun(term[[1]]))))
        stop(term[[1]], " has no inst argumnt")
    termList <- vector(mode = "list", length = instances)
    for (i in seq(instances)) {
        termList[[i]] <- term
        termList[[i]]$inst <- i
    }
    paste(unlist(termList), collapse = " + ")
}
