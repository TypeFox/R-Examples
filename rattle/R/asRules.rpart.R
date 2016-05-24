# Rattle: A GUI for Data Mining in R
#
# RPART RULES
#
# Time-stamp: <2014-09-05 21:27:43 gjw>
#
# Copyright (c) 2009-2014 Togaware Pty Ltd
#
# This files is part of Rattle.
#
# Rattle is free software: you can redistribute it and/or modify it
# under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 2 of the License, or
# (at your option) any later version.
#
# Rattle is distributed in the hope that it will be useful, but
# WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
# General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with Rattle. If not, see <http://www.gnu.org/licenses/>.

asRules <- function(model, compact=FALSE, ...) UseMethod("asRules")

asRules.rpart <- function(model, compact=FALSE, ...)
{
  if (!inherits(model, "rpart")) stop(Rtxt("Not a legitimate rpart tree"))
  # if (model$method != "class")) stop("Model method needs to be class")
  #
  # Get some information.
  #
  rtree <- length(attr(model, "ylevels")) == 0
  target <- as.character(attr(model$terms, "variables")[2])
  frm <- model$frame
  names <- row.names(frm)
  ylevels <- attr(model, "ylevels")
  ds.size <-  model$frame[1,]$n
  #
  # Print each leaf node as a rule.
  #
  if (rtree)
    # Sort rules by coverage
    ordered <- rev(sort(frm$n, index=TRUE)$ix)
  else
    # Sort rules by probabilty of second class (usually the last in binary class)
    ordered <- rev(sort(frm$yval2[,5], index=TRUE)$ix)
  for (i in ordered)
  {
    if (frm[i,1] == "<leaf>")
    {
      # The following [,5] is hardwired and works on one example....
      if (rtree)
        yval <- frm[i,]$yval
      else
        yval <- ylevels[frm[i,]$yval]
      cover <- frm[i,]$n
      pcover <- round(100*cover/ds.size)
      if (! rtree) prob <- frm[i,]$yval2[,5]
      cat("\n")
      pth <- rpart::path.rpart(model, nodes=as.numeric(names[i]), print.it=FALSE)
      pth <- unlist(pth)[-1]
      if (! length(pth)) pth <- "True"
      if (compact)
      {
        cat(sprintf("R%03s ", names[i]))
        if (rtree)
          cat(sprintf("[%2.0f%%,%0.2f]", pcover, prob))
        else
          cat(sprintf("[%2.0f%%,%0.2f]", pcover, prob))
        cat(sprintf(" %s", pth), sep="")
      }
      else
      {
        cat(sprintf(Rtxt(" Rule number: %s "), names[i]))
        if (rtree)
          cat(sprintf("[%s=%s cover=%d (%.0f%%)]\n",
                      target, yval, cover, pcover))
        else
          cat(sprintf("[%s=%s cover=%d (%.0f%%) prob=%0.2f]\n",
                      target, yval, cover, pcover, prob))
        cat(sprintf("   %s\n", pth), sep="")
      }
    }
  }
  cat("\n")
  invisible(ordered)
}

