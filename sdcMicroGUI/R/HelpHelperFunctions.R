#  Original Function in R SVN src/library/tools/R/Rd.R
#  Original Function in R SVN src/library/utils/R/help.R
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

.getHelpFile_sdcMicroGUI <- function(file)
{
  path <- dirname(file)
  dirpath <- dirname(path)
  if(!file.exists(dirpath))
    stop(gettextf("invalid %s argument", sQuote("file")), domain = NA)
  pkgname <- basename(dirpath)
  RdDB <- file.path(path, pkgname)
  if(!file.exists(paste(RdDB, "rdx", sep = ".")))
    stop(gettextf("package %s exists but was not installed under R >= 2.10.0 so help cannot be accessed", sQuote(pkgname)), domain = NA)
  fetchRdDB_sdcMicroGUI(RdDB, basename(file))
}

fetchRdDB_sdcMicroGUI <- function(filebase, key = NULL)
{
  fun <- function(db) {
    vals <- db$vals
    vars <- db$vars
    datafile <- db$datafile
    compressed <- db$compressed
    envhook <- db$envhook
    
    fetch <- function(key)
      lazyLoadDBfetch(vals[key][[1L]], datafile, compressed, envhook)
    
    if(length(key)) {
      if(! key %in% vars)
        stop(gettextf("No help on %s found in RdDB %s",
                sQuote(key), sQuote(filebase)),
            domain = NA)
      fetch(key)
    } else {
      res <- lapply(vars, fetch)
      names(res) <- vars
      res
    }
  }
  res <- lazyLoadDBexec(filebase, fun)
  if (length(key))
    res
  else
    invisible(res)
}