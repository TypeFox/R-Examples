# Gnome R Data Miner: GNOME interface to R for Data Mining
#
# Time-stamp: <2015-05-17 08:55:47 gjw>
#
# Implement functionality associated with the Execute button and Menu.
#
# Copyright (c) 2009-2013 Togaware Pty Ltd
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

on_execute_button_clicked <- function(action, window)
{
  # 100402 Allow Execute to be running just once, irrespective of the
  # number of times the Execute button is clicked. Otherwise we get a
  # second load of a CSV dataset whilst still loading the first.
  
  if (! is.null(crv$executing) && crv$executing) return()
  crv$executing <- TRUE
  on.exit(crv$executing <- FALSE)

  # Wrap up the actual call with a "try" so that the watch cursor
  # turns off even on error.

  setStatusBar()

  # 081117 This ensures spinbuttons, for example, lose focus and hence
  # their current value is properly noted. Otherwise I was finding the
  # user had to either press Enter or click somewhere else to ensure
  # the value is noted.

  theWidget("rattle_window")$setFocus()

  # 090102 Set the cursor to busy, and make sure on failure or
  # interrupt we set it back. TODO Currently, I can interrupt with
  # Ctrl-C in the console, and that does interrupt the Rattle process,
  # but I can't work out how to get a Ctrl-C (or perhaps an ESC) in
  # the Rattle GUI to cause an interrupt.

  set.cursor("watch")
  tryCatch(dispatchExecuteButton(),
           interrupt=function(m) setStatusBar(Rtxt("Processing interrupted by user.")),
           finally=set.cursor())

#  library(multicore)
#  set.cursor("watch")
#  crs$process <- parallel(dispatchExecuteButton())
#                                   interrupt=function(m)
#                                   setStatusBar("Processing interrupted by user."),
#                                   finally=set.cursor()))
#  collect()

  # 090103 Return nothing, otherwise we get the results from the
  # tryCatch above.

  return()
}

dispatchExecuteButton <- function()
{
  # Check which tab of notebook and dispatch to appropriate execute action

  ct <- getCurrentPageLabel(crv$NOTEBOOK)
  # REMOVE 100424 No longer required here - this is done earlier now.
  # Encoding(ct) <- "UTF-8" # 100408 For French, but see if it's okay always!

  if (ct == crv$NOTEBOOK.DATA.NAME)
  {
    executeDataTab()
  }
  else if (ct == crv$NOTEBOOK.EXPLORE.NAME)
  {
    executeExploreTab()
  }
  else if (ct == crv$NOTEBOOK.TEST.NAME)
  {
    executeTestTab()
  }
  else if (ct == crv$NOTEBOOK.TRANSFORM.NAME)
  {
    executeTransformTab()
  }
  else if (ct == crv$NOTEBOOK.CLUSTER.NAME)
  {
    executeClusterTab()
  }
  else if (ct == crv$NOTEBOOK.ASSOCIATE.NAME)
  {
    executeAssociateTab()
  }
  else if (ct == crv$NOTEBOOK.MODEL.NAME)
  {
    executeModelTab()
  }
  else if (ct == crv$NOTEBOOK.EVALUATE.NAME)
  {

    # The wrap mode of the confusion_textview may have been set to
    # word wrap when a model was Executed if it had more than 2
    # classes, since a message is printed about ROCR etc not handling
    # any more than 2 classes.

    theWidget("confusion_textview")$setWrapMode("none")
    executeEvaluateTab()
  }
  else if (ct == crv$NOTEBOOK.LOG.NAME)
  {
    executeLogTab()
  }
  else
  {
    errorDialog(Rtxt("'dispatchExecuteButton' has been called with an unknown tab."),
                "\n\n", ct,
                "\n\n", crv$support.msg)
    return()
  }
}
