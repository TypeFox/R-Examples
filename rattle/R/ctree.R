# Gnome R Data Miner: GNOME interface to R for Data Mining
#
# Time-stamp: <2015-05-17 08:54:45 gjw>
#
# CTREE OPTION OF THE TREE TAB
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

########################################################################
#
# Model -> Tree -> Conditional
#

# 100815 TODO The "partykit" package from R-Forge (only for now)
# includes .list.rules.party() to convert tree into rules:
#
# install.packages("partykit", repos = "http://R-Forge.R-project.org")
# library("partykit")
# Rebuild the ctree as partykit provides new ctree.
# partykit:::.list.rules.party(crs$rpart)

executeModelCTree <- function()
{
  # 080815 This is currently just copied from rpart.R, and slowly
  # being tuned for ctree specifically.
  
  # Initial setup 

  TV <- "rpart_textview"

  num.classes <- length(levels(as.factor(crs$dataset[[crs$target]])))
  control <- NULL
  parms <- NULL

  # Scrape the value of the tuning controls

  tune.controls <- theWidget("rpart_tune_entry")$getText()
  
  # Retrieve the Priors, and check there is the right number and that
  # they add up to 1.
  
  priors <- theWidget("model_tree_priors_entry")$getText()
  if (nchar(priors) > 0)
  {
    pr <- as.numeric(unlist(strsplit(priors, ",")))
    if (length(pr) != num.classes)
      {
        errorDialog(sprintf(Rtxt("The supplied priors (%s)",
                                 "need to correspond to the number of classes",
                                 "found in the target variable '%s'.",
                                 "Please supply exactly %d priors."),
                            priors, crs$target, num.classes))
        return(FALSE)
      }
    if (sum(pr) != 1)
      {
        errorDialog(sprintf(Rtxt("The supplied priors (%s)",
                                 "add up to %0.2f whereas",
                                 "they need to add up 1.00.",
                                 "Please provide appropriate priors."),
                            priors, sum(pr)))
        return(FALSE)
      }
    if (is.null(parms))
      parms <- sprintf(", parms=list(prior=c(%s))", priors)
    else
      parms <- gsub(")$", sprintf(", prior=c(%s)", priors), parms)
  }

  # Retrieve the Min Split and check if it is different from the
  # default, and if so then use it.

  minsplit <- theWidget("rpart_minsplit_spinbutton")$getValue()
  if (minsplit != crv$rpart.minsplit.default)
  {
    if (is.null(control))
      control <- sprintf(", control=ctree_control(minsplit=%d)", minsplit)
    else
      control <- gsub(")$", sprintf(", minsplit=%d)", minsplit), control)
  }

  # Retrieve the Min Bucket and check if it is different from the
  # default, and if so then use it.

  minbucket <- theWidget("rpart_minbucket_spinbutton")$getValue()
  if (minbucket != crv$rpart.minbucket.default)
  {
    if (is.null(control))
      control <- sprintf(", control=ctree_control(minbucket=%d)", minbucket)
    else
      control <- gsub(")$", sprintf(", minbucket=%d)", minbucket), control)
  }

  # Retrieve the Max Depth and check if it is different from the
  # default, and if so then use it.

  maxdepth <- theWidget("rpart_maxdepth_spinbutton")$getValue()
  if (maxdepth != crv$rpart.maxdepth.default)
  {
    if (is.null(control))
      control <- sprintf(", control=ctree_control(maxdepth=%d)", maxdepth)
    else
      control <- gsub(")$", sprintf(", maxdepth=%d)", maxdepth), control)
  }

  # Build the formula for the model.

  frml <- paste(crs$target, "~ .")

  # Variables to be included --- a string of indicies.
  
  # included <- getIncludedVariables()
  included <- "c(crs$input, crs$target)" # 20110102
  
  # Some convenience booleans

  sampling  <- not.null(crs$sample)
  including <- not.null(included)
  subsetting <- sampling || including
  
  # Commands.
  
  lib.cmd <- "library(party, quietly=TRUE)"
  if (! packageIsAvailable("party", Rtxt("build conditional trees"))) return(FALSE)

  fit.cmd <- paste("crs$rpart <- ctree(", frml, ", data=crs$dataset",
                   if (subsetting) "[",
                   if (sampling) "crs$sample",
                   if (subsetting) ",",
                   if (including) included,
                   if (subsetting) "]",
                   if (! is.null(crs$weights))
                   sprintf(",\n    weights=as.integer(%s)%s",
                           crs$weights,
                           ifelse(sampling, "[crs$train]", "")),
                   ifelse(is.null(control), "", control),
                   ")", sep="")

  print.cmd <- "print(crs$rpart)"
                               
  # Load the required library.

  startLog(Rtxt("Conditional inference tree."))
  appendLog(Rtxt("Build a conditional tree using the party package."), lib.cmd)

  eval(parse(text=lib.cmd))

  # Build the model.

  appendLog(Rtxt("Build a ctree model."), fit.cmd)
  start.time <- Sys.time()
  result <- try(eval(parse(text=fit.cmd)), silent=TRUE)
  time.taken <- Sys.time()-start.time
  if (inherits(result, "try-error"))
  {
    errorDialog(errorMessageFun("ctree", result))
    return(FALSE)
  }

  # Display the resulting model.

  appendLog(Rtxt("Generate summary of the ctree model."), print.cmd)

  resetTextview(TV)
  setTextview(TV,
              sprintf(Rtxt("Summary of the %s model for %s (built using '%s'):\n"),
                      commonName("ctree"),
                      Rtxt("Classification"), # 080604 TODO put the right type
                      "ctree"),
              collectOutput(print.cmd), "\n")

  if (sampling) crs$smodel <- union(crs$smodel, crv$RPART)

  # Now that we have a model, make sure the rules and plot buttons are
  # not visible.
  
  showModelRPartExists()

  # Finish up.

  reportTimeTaken(TV, time.taken, model=commonName(crv$RPART))

  return(TRUE)
}
