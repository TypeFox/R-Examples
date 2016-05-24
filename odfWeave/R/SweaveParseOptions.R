#  From src/library/utils/R/Sweave.R
#  Part of the R package, http://www.R-project.org
#
#  Copyright (C) 1995-2012 The R Core Team
#
#  This program is free software; you can redistribute it and/or modify
#  it under the terms of the GNU General Public License as published by
#  the Free Software Foundation; either version 2 of the License, or
#  (at your option) any later version.
#
#  This program is distributed in the hope that it will be useful,
#  but WITHOUT ANY WARRANTY; without even the implied warranty of
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#  GNU General Public License for more details.
#
#  A copy of the GNU General Public License is available at
#  http://www.r-project.org/Licenses/


SweaveParseOptions <- function(text, defaults = list(), check = NULL)
{
  x <- sub("^[[:space:]]*(.*)", "\\1", text)
  x <- sub("(.*[^[:space:]])[[:space:]]*$", "\\1", x)
  x <- unlist(strsplit(x, "[[:space:]]*,[[:space:]]*"))
  x <- strsplit(x, "[[:space:]]*=[[:space:]]*")
  
  ## only the first option may have no name: the chunk label
  if (length(x)) {
    if (length(x[[1L]]) == 1L) x[[1L]] <- c("label", x[[1L]])
  } else return(defaults)
  
  if (any(sapply(x, length) != 2L))
    stop(gettextf("parse error or empty option in\n%s", text), domain = NA)
  
  options <- defaults
  for (k in seq_along(x)) options[[ x[[k]][1L] ]] <- x[[k]][2L]
  
  ## This is undocumented
  if (!is.null(options[["label"]]) && !is.null(options[["engine"]]))
    options[["label"]] <-
    sub(paste0("\\.", options[["engine"]], "$"),
        "", options[["label"]])
  
  if (!is.null(check)) check(options) else options
}
