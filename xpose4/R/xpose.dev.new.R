# Xpose 4
# An R-based population pharmacokinetic/
# pharmacodynamic model building aid for NONMEM.
# Copyright (C) 1998-2004 E. Niclas Jonsson and Mats Karlsson.
# Copyright (C) 2005-2008 Andrew C. Hooker, Justin J. Wilkins, 
# Mats O. Karlsson and E. Niclas Jonsson.
# Copyright (C) 2009-2010 Andrew C. Hooker, Mats O. Karlsson and 
# E. Niclas Jonsson.

# This file is a part of Xpose 4.
# Xpose 4 is free software; you can redistribute it and/or
# modify it under the terms of the GNU Lesser General Public License
# as published by the Free Software Foundation, either version 3
# of the License, or (at your option) any later version.

# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU Lesser General Public License for more details.

# You should have received a copy of the GNU Lesser General Public License
# along with this program.  A copy can be cound in the R installation
# directory under \share\licenses. If not, see http://www.gnu.org/licenses/.

xpose.dev.new <- function (...) 
{
  if(getRversion()>="2.8.0"){
    dev.new(...)
  } else {
    dev <- getOption("device")
    if (!is.character(dev) && !is.function(dev)) 
      stop("invalid setting for 'getOption(\"device\")'")
    if (is.character(dev)) {
      dev <- if (exists(dev, .GlobalEnv)) 
        get(dev, .GlobalEnv)
      else if (exists(dev, asNamespace("grDevices"))) 
        get(dev, asNamespace("grDevices"))
      else stop(gettextf("device '%s' not found", dev), domain = NA)
    }
    a <- list(...)
    a2 <- names(formals(dev))
    a <- a[names(a) %in% a2]
    if (identical(dev, pdf)) {
      if (is.null(a[["file"]]) && file.exists("Rplots.pdf")) {
        fe <- file.exists(tmp <- paste("Rplots", 1:999, ".pdf", 
                                       sep = ""))
        if (all(fe)) 
          stop("no suitable unused file name for pdf()")
        message(gettextf("dev.new(): using pdf(file=\"%s\")", 
                         tmp[!fe][1]), domain = NA)
        a$file <- tmp[!fe][1]
      }
    }
    else if (identical(dev, postscript)) {
      if (is.null(a[["file"]]) && file.exists("Rplots.ps")) {
        fe <- file.exists(tmp <- paste("Rplots", 1:999, ".ps", 
                                       sep = ""))
        if (all(fe)) 
          stop("no suitable unused file name for postscript()")
        message(gettextf("dev.new(): using postscript(file=\"%s\")", 
                         tmp[!fe][1]), domain = NA)
        a$file <- tmp[!fe][1]
      }
    }
    else if (!is.null(a[["width"]]) && !is.null(a[["height"]]) && 
             (identical(dev, png) || identical(dev, jpeg) || identical(dev, 
                                                                       bmp) || identical(dev, tiff))) {
      if (is.null(a[["units"]]) && is.null(a[["res"]])) {
        a$units <- "in"
        a$res <- 72
      }
    }
    do.call(dev, a)
  }
}
#<environment: namespace:grDevices>
