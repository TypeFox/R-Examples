# Gnome R Data Miner: GNOME interface to R for Data Mining
#
# Time-stamp: <2015-05-17 08:58:38 gjw>
#
# Reporting support
#
# Copyright (c) 2009 Togaware Pty Ltd
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

# TODO
# 100307 Consider moving to using reporttools and using Sweave instead

on_report_toolbutton_clicked <- function(action, window)
{
  # Wrap the actual call with a "try" so that the watch cursor turns
  # off even on error.
  
  setStatusBar("Generating report.")
  set.cursor("watch")
  try(dispatchReportButton())
  set.cursor()
}

dispatchReportButton <- function()
{
  # Prerequisites: Can not report on data if there is no dataset.

  if (noDatasetLoaded()) return(FALSE)

  if (! questionDialog("The Report button is very experimental.",
                       "Please report issues and updates to",
                       "support@togaware.com.",
                       "\n\nKnown issues:",
                       "\n\n\tAlways saves to the same fixed file",
                       "- need chooser.",
                       "\n\tA plot is displayed on screen - need to suppress.",
                       "\n\tToo much generated to the console - how remove?",
                       "\n\nOtherwise it is safe to use!",
                       "\n\nDo you wish to continue?"))
    return(FALSE)
  
  startLog("GENERATE A REPORT")
  
  if (! packageIsAvailable("odfWeave", "generate a report")) return(FALSE)
  lib.cmd <- "library(odfWeave, quietly=TRUE)"
  appendLog("The odfWeave package processes ODT document templates.", lib.cmd)
  eval(parse(text=lib.cmd))
  
  # Check which tab of the notebook is active and dispatch to the
  # appropriate execute action.

  ct <- getCurrentPageLabel(crv$NOTEBOOK)

  if (ct == crv$NOTEBOOK.DATA.NAME ||
      ct == crv$NOTEBOOK.EXPLORE.NAME)
  {
    # For the DATA or EXPLORE tabs generate a dataset summary.
    
    reportDataTab()
  }
  else if (ct == crv$NOTEBOOK.MODEL.NAME )
  {
    if (! is.null(crs$rpart))
      reportTreeModel(crs$rpart)
    else
    {
      infoDialog("Report functionality is only available for the Tree",
                 ct, " and no Tree model found.")
      return(FALSE)
    }
  }
  else
      
  {
    infoDialog("No report functionality is available for the",
               ct, "tab as yet. Nothing done.")
    return(FALSE)
  }
}

#-----------------------------------------------------------------------

reportDataTab <- function()
{
  if (file.exists("../odf/data_summary.odt"))
  {
    summary <- "../odf/data_summary.odt" # For Testing
    warning(Rtxt("Rattle Report is using local template ../odf"), immediate.=TRUE)
  } 
  else
    summary <- system.file("odt", "data_summary.odt", package="rattle")

  ofile <- paste(getwd(), "data_summary_rattle.odt", sep="/")

  odf.cmd <- sprintf(paste('odfWeave("%s",',
                           '\n         "%s",',
                           '\n         control = odfWeaveControl(verbose = FALSE))'),
                     summary, ofile)

  appendLog(Rtxt("Generate a data report."), odf.cmd)

  eval(parse(text=odf.cmd))

  setStatusBar(sprintf("Report written to %s.", ofile))

  system(paste("oowriter", ofile), wait = FALSE)
}

  
#-----------------------------------------------------------------------

reportTreeModel <- function(model)
{
  model <<- model
  if (file.exists("../odf/model_rpart_summary.odt"))
  {
    summary <- "../odf/model_rpart_summary.odt" # For Testing
    warning("Rattle Report is using local template ../odf", immediate.=TRUE)
  } 
  else
    summary <- system.file("odt", "mode_rpart_summary.odt", package="rattle")

  ofile <- paste(getwd(), "model_rpart_summary_rattle.odt", sep="/")
  
  odfWeave::odfWeave(summary, ofile, control=odfWeave::odfWeaveControl(verbose=FALSE))

  if (! is.null(crv$rattleGUI)) setStatusBar(sprintf("Report written to %s.", ofile))
}
