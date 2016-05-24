#  Copyright (C) 2005, 2006, 2008 David Firth and Heather Turner
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

Symm <- function(...){
    dots <- list(...)
    if (any(diff(sapply(dots, length)) != 0)) stop(
                "arguments to symm() must all have same length")
    dots <- lapply(dots, as.factor)
    facMatrix <- sapply(dots, as.character)
    f <- function(row){
        string <- paste(sort(row), collapse = "")
        if (any(is.na(row))) is.na(string) <- TRUE
        string
    }
    if (inherits(facMatrix, "matrix"))
        result <- factor(apply(facMatrix, 1, f))
    else
        result <- factor(f(facMatrix))
    result
}
