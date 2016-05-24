# Gnome R Data Miner: GNOME interface to R for Data Mining
#
# Time-stamp: <2014-09-07 05:23:16 gjw>
#
# Textview widget support
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
#
# TODO
#
#	We are in the middle of migrating from using a string to
#	identify the widget, and assuming theWidget will get the
#	actual widget (which only works for the main crv$rattleGUI
#	window) to passing the actual widget itself, which is much
#	more general. For now, allow both, through the use of
#	getTextview.
#
#       071128 New version of R (since about 6.2.0) have a default
#       value for useFancyFonts that does not work with textviews in
#       RGtk2 (and probably other things in RGtk2). This affects only
#       MSWindows. A quick fix, which has been done for setTextview,
#       but not the rest yet, is to wrap textview displays with:
#
#         oldopt <- options(useFancyQuotes="utf8")
#         ...
#         options(oldopt)
#

allTextviews <- function()
{
  # 100407 The log_textview is not included here. This is probably
  # because we don't want to reset it anytime, so that it is a
  # copmlete log. But then it's font is not being set to monospace
  # like all the rest. So that needs a special initialisation.
  
  return(c("summary_textview", "interactive_textview",
           "correlation_textview",
           "prcomp_textview", "test_textview", "kmeans_textview",
           "hclust_textview", "associate_textview", "rpart_textview",
           "glm_textview", "ada_textview", "rf_textview",
           "esvm_textview", "ksvm_textview", "nnet_textview",
           "model_survival_textview",
           "confusion_textview", "risk_textview", "roc_textview"))
}


getTextview <- function(tv)
{
  # From either a string or the actual object, return the textview
  # object. This is usful internally in this file whilst we migrate
  # to not using the string to name the textview, but passing the
  # object itself.
  
  wid <- FALSE
  if (inherits(tv, "GtkTextView"))
    wid <- tv
  else if (inherits(tv, "character"))
    wid <- theWidget(tv)
  return(wid)
}

resetTextview <- function(tv, ..., sep="", tvsep=TRUE)
{
  # We call this if we simply want to blank the textview, and
  # optionally, add new text to the textview. We use a default
  # separator of nothing, rather than a single space, since this
  # suits the usual list of strings being passed into the textview,
  # where I have a string introducing what follows, then there is the
  # output of some command. By default we also add a Rattle separator
  # line, as the usual usage of resetTexview will be to clear the
  # textview and add a message, in one go. There are times when we
  # don't want the separator though.
  
  if (is.null(wid <- getTextview(tv)))
  {
    errorDialog("E138: Should not be here.",
                "The textview object supplied to resetTextview",
                "is neither a GtkTextView nor a string.",
                "We found a", class(tv)[1], crv$support.msg)
    return(FALSE)
  }
  if (! isJapanese()) wid$modifyFont(RGtk2::pangoFontDescriptionFromString(crv$textview.font))
  msg <- paste(sep=sep, ...)
  if (length(msg) > 0)
  {
    wid$getBuffer()$setText(msg)
    if (tvsep) appendTextview(tv)
  }
  else
    wid$getBuffer()$setText("")

}

appendTextview <- function(tv, ..., sep="", tvsep=TRUE)
{
  # Append a message to the given textview. Optionally add a Rattle
  # separator to the textview. By default, paste the strings of the
  # message together without a speartor.
  
  if (is.null(wid <- getTextview(tv)))
  {
    errorDialog("E140: Should not be here.",
                "The textview object supplied to appendTextview",
                "is neither a GtkTextView nor a string.",
                "We found a", class(tv)[1], crv$support.msg)
    return(FALSE)
  }
  msg <- paste(sep=sep, ...)
  if (tvsep) msg <- paste(msg, textviewSeparator(), sep="")
  buf <- wid$getBuffer()
  location <- buf$getEndIter()$iter
  buf$insert(location, msg)
}

textviewSeparator <- function()
{
  return(paste("\n\n",
               if (not.null(crv$show.timestamp) && crv$show.timestamp)
                 sprintf("%s %s %s %s\n", crv$appname, Rtxt("timestamp:"), Sys.time(),
                         Sys.info()["user"]),
               paste(rep("=", 70), collapse=""), "\n", sep=""))
}

getTextviewContent <- function(TV)
{
  # Extract text contents of specified textview and return
  # it. Designed for use in saveProject.
  
  log.buf <- theWidget(TV)$getBuffer()
  start <- log.buf$getStartIter()$iter
  end <- log.buf$getEndIter()$iter
  return(log.buf$getText(start, end))
}

# STOP USING THE FOLLOWING

setTextview <- function(tv, ..., sep="")
{
  # Stop using this - use resetTextview instead

  if (is.null(wid <- getTextview(tv)))
  {
    errorDialog("E137: Should not be here.",
                "The textview object supplied to setTextview",
                "is neither a GtkTextView nor a string.",
                "We found a", class(tv)[1], crv$support.msg)
    return(FALSE)
  }
  oldopt <- options(useFancyQuotes="utf8") # Bug fix for MSWindows [071128]
  msg <- paste(sep=sep, ...)
  if (length(msg) == 0) msg <-""
  wid$getBuffer()$setText(msg)
  options(oldopt) # Bug fix for MSWindows [071128]
}

addTextview <- function(tv, ..., sep="")
{
  if (is.null(wid <- getTextview(tv)))
  {
    errorDialog("E139: Should not be here.",
                "The textview object supplied to addTextview",
                "is neither a GtkTextView nor a string.",
                "We found a", class(tv)[1], crv$support.msg)
    return(FALSE)
  }
  msg <- paste(sep=sep, ...)
  if (length(msg) == 0) msg <-""
  tv.buf <- wid$getBuffer()
  loc <- tv.buf$getEndIter()$iter
  tv.buf$insert(loc, msg)
}
  
setTextviewContents <- function(TV, text)
{
  # Set the text contents of the specified textview to the supplied
  # text. Designed for use in loadProject.

  resetTextview(TV)
  theWidget(TV)$setWrapMode("none")
  if (is.null(text))
    theWidget(TV)$getBuffer()$setText("")
  else
    theWidget(TV)$getBuffer()$setText(text)
}

resetTextviews <- function(tv=NULL)
{
  # 090202 Reset all text views to default content, as when Rattle
  # starts up or the user has selected New Project, or has loaded a
  # new dataset. The text for the texviews come from textviews.xml.
  
  # 090214 I probably actually just want to go through a clear each of
  # the textviews first, and then populate with the text from
  # textviews.xml if there is one and XML package is available!

  if (is.null(tv))
    sapply(allTextviews(), resetTextview)
  
  if (! packageIsAvailable("XML", "load textview texts"))
  {
    warning("The XML package is not available. Textview texts will not be available.")
    return(FALSE)
  }

  result <- try(etc <- file.path(path.package(package="rattle")[1], "etc"),
                silent=TRUE)
  if (inherits(result, "try-error"))
    doc <- XML::xmlTreeParse("textviews.xml", useInternalNodes=TRUE)
  else
    doc <- XML::xmlTreeParse(file.path(etc, "textviews.xml"), useInternalNodes=TRUE)

  if (is.null(tv))
    sapply(XML::getNodeSet(doc, "//textview"),
           function(tt)
           {
             wd <- XML::xmlGetAttr(tt, 'widget')
             resetTextview(wd, Rtxt(XML::xmlValue(tt)), tvsep=FALSE)
           })
  else
    sapply(XML::getNodeSet(doc, "//textview"),
           function(tt)
           {
             wd <- XML::xmlGetAttr(tt, 'widget')
             if (wd %in% tv) resetTextview(wd, Rtxt(XML::xmlValue(tt)), tvsep=FALSE)
           })
  invisible()
}
