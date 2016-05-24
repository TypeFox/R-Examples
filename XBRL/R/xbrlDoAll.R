## -*- ess-indent-level: 2; ess-basic-offset: 2; tab-width: 8 -*-
##
## Copyright (C) 2014-2015 Roberto Bertolusso and Marek Kimmel
##
## This file is part of XBRL.
##
## XBRL is free software: you can redistribute it and/or modify it
## under the terms of the GNU General Public License as published by
## the Free Software Foundation, either version 2 of the License, or
## (at your option) any later version.
##
## XBRL is distributed in the hope that it will be useful, but
## WITHOUT ANY WARRANTY; without even the implied warranty of
## MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
## GNU General Public License for more details.
##
## You should have received a copy of the GNU General Public License
## along with XBRL. If not, see <http://www.gnu.org/licenses/>.


xbrlDoAll <- function(file.inst, cache.dir=NULL, prefix.out=NULL, verbose=FALSE) {
  xbrl <- XBRL()
  xbrl$setVerbose(verbose)
  if (!is.null(cache.dir)) { 
    xbrl$setCacheDir(cache.dir)
  }
  xbrl$openInstance(file.inst)
  xbrl$processSchema(xbrl$getSchemaName())
  xbrl$processContexts()
  xbrl$processFacts()
  xbrl$processUnits()
  xbrl$processFootnotes()
  xbrl$closeInstance()
  
  xbrl.vars <- xbrl$getResults()
  if (!is.null(prefix.out)) {
    if (verbose) {
      cat("Saving data\n")
    }
    write.csv(xbrl.vars$role, file=paste0(prefix.out, "_roles.csv"))
    write.csv(xbrl.vars$element, file=paste0(prefix.out, "_elements.csv"))
    write.csv(xbrl.vars$label, file=paste0(prefix.out, "_labels.csv"))
    write.csv(xbrl.vars$presentation, file=paste0(prefix.out, "_presentations.csv"))
    write.csv(xbrl.vars$definition, file=paste0(prefix.out, "_definitions.csv"))
    write.csv(xbrl.vars$calculation, file=paste0(prefix.out, "_calculations.csv"))
    write.csv(xbrl.vars$context, file=paste0(prefix.out, "_contexts.csv"))
    write.csv(xbrl.vars$fact, file=paste0(prefix.out, "_facts.csv"))
    write.csv(xbrl.vars$unit, file=paste0(prefix.out, "_units.csv"))
    write.csv(xbrl.vars$footnote, file=paste0(prefix.out, "_footnotes.csv"))
  }
  invisible(xbrl.vars)
}
