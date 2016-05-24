# Gnome R Data Miner: GNOME interface to R for Data Mining
#
# Time-stamp: <2015-05-17 08:54:03 gjw>
#
# Implement associations functionality.
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
# TODO
#
# 100308 Implement arulesSequences

########################################################################
#
# CALLBACKS
#

on_tools_associate_activate <- function(action, window)
{
  crv$NOTEBOOK$setCurrentPage(getNotebookPage(crv$NOTEBOOK,
                                              crv$NOTEBOOK.ASSOCIATE.NAME))
  switchToPage(crv$NOTEBOOK.ASSOCIATE.NAME)
}

on_associate_plot_frequency_button_clicked <-  function(action, window)
{
  plotAssociateFrequencies()
}

on_associate_rules_button_clicked <-  function(action, window)
{
  listAssociateRules()
}

on_associate_plot_button_clicked <-  function(action, window)
{
  plotAssociateRules()
}

########################################################################
# SUPPORT

generateAprioriSummary <- function(ap)
{
 result <- sprintf(paste(Rtxt("Number of Rules: %d"), "\n\n"), length(ap))
 # 080726 Started failing.... no size for ap@lhs
 ## result <- paste(result, "Rule Length Distribution (LHS + RHS):\n\n", sep="")
 ## pp <- table(size(ap@lhs)+size(ap@rhs))
 ## result <- paste(result, paste(sprintf("\tRules of length %s:  %-4d",
 ##                                          names(pp), pp), collapse="\n"),
 ##                "\n\n", sep="")
 return(result)
}

########################################################################
# EXECUTION

executeAssociateTab <- function()
{
  # We require a dataset.

  if (noDatasetLoaded()) return()

  # If it looks like the DATA page has not been executed, complain.

  if (variablesHaveChanged(Rtxt("identifying association rules"))) return()

  # Check if sampling needs executing.

  if (sampleNeedsExecute()) return()

  # Determine whether we want basket analysis of transactions or
  # rules from the categorical variables. The former is indicated if
  # there is a single IDENT variable and a TARGET with multiple
  # values. Perhaps we just see if the GUI checkbutton is set and if
  # so, check that the variables meet these criteria, and if not the
  # return. Also, on  executeSelectTab, if there is one ID and the
  # TARGET is factor or integer and there are no inputs then set the
  # default to Baskets.

  baskets <- theWidget("associate_baskets_checkbutton")$getActive()
  if (baskets && length(crs$ident) != 1)
  {
    errorDialog(sprintf(Rtxt("Exactly one variable must be identified as an Ident",
                             "in the Data tab to be used as",
                             "the identifier of the transactions.",
                             "There are %d Ident variables.",
                             "The observations need to be aggregated by the Ident to",
                             "create the baskets for association analysis."),
                        length(crs$ident)))
    return()
  }
  if (baskets && length(crs$target) != 1)
  {
    errorDialog(Rtxt("A Target variable must be identified in the Data tab.",
                     "This variable then identifies the items associated with each",
                     "basket or transaction in the analysis. Each basket or",
                     "transaction is uniquely identified using the Ident variable."))
    return()
  }
      
  # Check that we have only categorical attributes.

  include <- getCategoricVariables("names")
  
  if (!baskets && length(include) == 0)
  {
    errorDialog(Rtxt("Associations are calculated only for categoric data.",
                     "\n\nNo categoric variables were found in the dataset",
                     "from amongst those having an input role.",
                     "\n\nIf you wanted a basket analysis with the Target variable",
                     "listing the items, and the Ident variable identifying",
                     "the baskets, then please select the Baskets option."))
    return()
  }

  # 080724 Check that we have more than just one categoric
  # variable. With only one we get an error: no method or default for
  # coercing "factor" to "transactions". Makes sense to need at least
  # two categorics for associating.

  if (!baskets && length(include) == 1)
  {
    errorDialog(sprintf(Rtxt("Associations, when not using the Baskets option,",
                             "can only be identified when there are",
                             "multiple categoric variables.",
                             "\n\nOnly one categoric variable was found (%s)",
                             "from amongst those having an input role.",
                             "\n\nIf you wanted a basket analysis with the Target",
                             "variable listing the items, and the Ident variable",
                             "identifying the baskets, then please deselect the",
                             "Baskets button."),
                        include[1]))
    return()
  }

  # Ensure the arules library is available and loaded.

  if (! packageIsAvailable("arules", Rtxt("generate associations"))) return()

  startLog(commonName("arules"))
  lib.cmd <- "library(arules, quietly=TRUE)"
  appendLog(packageProvides("arules", "arules"), lib.cmd)
  eval(parse(text=lib.cmd))
 
  # Initialise the textview.
  
  TV <- "associate_textview"
  
  # Required information.
  
  sampling   <- not.null(crs$sample)
  support    <- theWidget("associate_support_spinbutton")$getValue()
  confidence <- theWidget("associate_confidence_spinbutton")$getValue()
  minlen     <- theWidget("associate_minlen_spinbutton")$getValue()

  # Transform data into a transactions dataset for arules.

  include <- "crs$categoric" # 20110102 getCategoricVariables()
  
  if (baskets)
    transaction.cmd <- paste("crs$transactions <- as(split(",
                             sprintf('crs$dataset[%s, %s],\n%*scrs$dataset[%s, %s]',
                                     ifelse(sampling, "crs$sample", ""),
                                     "crs$target",
                                     29, "",
                                     ifelse(sampling, "crs$sample", ""),
                                     "crs$ident"),
                             '),',
                             sprintf('\n%*s"transactions")', 23, ""),
                             sep="") 
  else
    transaction.cmd <- paste("crs$transactions <- as(",
                             sprintf('crs$dataset[%s, %s], "transactions")',
                                     ifelse(sampling, "crs$sample", ""),
                                     include), sep="")
  appendLog(Rtxt("Generate a transactions dataset."), transaction.cmd)

  start.time <- Sys.time()
  result <- try(eval(parse(text=transaction.cmd)))
  time.taken <- Sys.time() - start.time
  if (inherits(result, "try-error"))
  {
    if (any(grep("slot i is not \\*strictly\\* increasing inside a column", result)) ||
        any(grep("transactions with duplicated items", result)))
      errorDialog(sprintf(Rtxt("One or more baskets appear to have repeated items.",
                               "This is not handled by the internal conversion to transactions.",
                               "Please pre-process your data to remove repeated items",
                               "from any baskets. The actual error message was:\n\n%s"),
                          result))
    else
      errorDialog(errorMessageFun("as", result))
    return(FALSE)
  }

  # Now generate the association rules.

  apriori.cmd <- paste("crs$apriori <- apriori(crs$transactions, ",
                       "parameter = list(",
                       sprintf("support=%.3f, confidence=%.3f, minlen=%d",
                               support, confidence, minlen),
                       "))", sep="")
  appendLog(Rtxt("Generate the association rules."), apriori.cmd)
  cmd.output <- collectOutput(apriori.cmd)

  # Add a summary of the transactions and the rules.

  tr.summary.cmd <- "summary(crs$transactions)"
  
  mysummary.cmd <- "generateAprioriSummary(crs$apriori)"
  appendLog(Rtxt("Summarise the resulting rule set."), mysummary.cmd)
  summary.cmd <- "summary(crs$apriori@quality)"

  resetTextview(TV)
  setTextview(TV,
              Rtxt("Summary of the Transactions:"), "\n\n",
              paste(capture.output(eval(parse(text=tr.summary.cmd))), collapse="\n"),
              "\n\n",
              Rtxt("Summary of the Apriori Association Rules:"), "\n\n",
              collectOutput(mysummary.cmd, use.cat=TRUE),
              "\n",
              Rtxt("Summary of the Measures of Interestingness:"), "\n\n",
              collectOutput(summary.cmd, use.print=TRUE),
              "\n\n",
              Rtxt("Summary of the Execution of the Apriori Command:"),
              "\n",
              cmd.output,
              "\n")
  
  reportTimeTaken(TV, time.taken, model=commonName("arules"))
}

plotAssociateFrequencies <- function()
{
  ## We require a dataset

  if (noDatasetLoaded()) return()

  ## If it looks like the VARIABLES page has not been executed, complain..

  if (variablesHaveChanged(Rtxt("building and then plotting association rules")))
    return(FALSE)

  ## Check if sampling needs executing.

  if (sampleNeedsExecute()) return()
    
  baskets <- theWidget("associate_baskets_checkbutton")$getActive()
  if (baskets && length(crs$ident) != 1)
  {
    errorDialog(sprintf(Rtxt("Exactly one variable must be identified as an Ident",
                             "in the Data tab to be used as the identifier of the transactions.",
                             "There were %s variables found.",
                             "The entities need to be aggregated by the Ident to",
                             "create the baskets for association analysis."),
                        length(crs$ident)))
    return()
  }
  if (baskets && length(crs$target) != 1)
  {
    errorDialog(Rtxt("You need to specify a Target variable in the Data tab.",
                     "This vairable then identifies the items associated with each",
                     "basket or transaction in the analysis. Each basket or",
                     "transaction is uniquely identified using the Ident variable."))
    return()
  }

  # Check that we have categorical attributes.

  include <- getCategoricVariables()
  if (! baskets && length(include) == 0)
  {
    errorDialog(Rtxt("Associations are calculated only for categoric variables.",
                     "No categorical variables were found in the dataset",
                     "from amongst those having an Input role.\n\n",
                     "If instead you wanted a basket analysis with the Target variable",
                     "listing the items, and the Ident variable identifying",
                     "the baskets, then please click the Baskets button."))
    return()
  }

  # Ensure the arules library is available and loaded.

  if (! packageIsAvailable("arules", Rtxt("generate associations"))) return()
  startLog(Rtxt("Relative Frequencies Plot"))
  lib.cmd <- "library(arules, quietly=TRUE)"
  appendLog(Rtxt("Association rules are implemented in the 'arules' package."), lib.cmd)
  eval(parse(text=lib.cmd))
 
  # Required information
  
  sampling  <- not.null(crs$sample)
  support <- theWidget("associate_support_spinbutton")$getValue()

  # Transform data into a transactions dataset for arules.

  # TODO 080425 Note that aprior does "as(..., 'transactions')" so we
  # don't need this here do we?
  
  if (baskets)
    transaction.cmd <- paste("crs$transactions <- as(split(",
                             sprintf('crs$dataset%s$%s, crs$dataset%s$%s',
                                     ifelse(sampling, "[crs$sample,]", ""),
                                     crs$target,
                                     ifelse(sampling, "[crs$sample,]", ""),
                                     crs$ident),
                             '), "transactions")', sep="")
  else
    transaction.cmd <- paste("crs$transactions <- as(",
                             sprintf('crs$dataset[%s,%s], "transactions")',
                                     ifelse(sampling, "crs$sample", ""),
                                     include), sep="")
  appendLog(Rtxt("Generate a transactions dataset."), transaction.cmd)
  eval(parse(text=transaction.cmd))

  # Now plot the relative frequencies.

  plot.cmd <- paste("itemFrequencyPlot(crs$transactions, support=",
                    support, ", cex=0.8)", sep="")
  newPlot()
  appendLog(Rtxt("Plot the relative frequencies."), plot.cmd)
  eval(parse(text=plot.cmd))

  setStatusBar(Rtxt("Generated the relative frequency plot."))
}

plotAssociateRules <- function()
{
  # We require a dataset

  if (noDatasetLoaded()) return()

  # Also make sure we have already generated the association rules.
  
  if (is.null(crs$apriori))
  {
    errorDialog(Rtxt("You first need to generate the association rules.",
                     "Perhaps you need to click the Execute button."))
    return()
  }

  # Ensure the required package is available

  if (! packageIsAvailable("arulesViz", Rtxt("generate a plot of association rules"))) return()
  startLog(Rtxt("Graph of Rules"))
  lib.cmd <- "library(arulesViz, quietly=TRUE)"
  appendLog(Rtxt("Use 'arulesViz' to plot a graph representing the rules."), lib.cmd)
  eval(parse(text=lib.cmd))

  plot.cmd <- paste('plot(crs$apriori, method="graph")\n',
#                    genPlotTitleCmd(Rtxt("Association Rules"), commonName(mtype),
#                                    testname, risk),
                    sep="")

  appendLog(Rtxt("Display a graph representing the rules."), plot.cmd)
  eval(parse(text=plot.cmd))
}

listAssociateRules <- function()
{
  # We require a dataset

  if (noDatasetLoaded()) return()

  # Also make sure we have already generated the association rules.
  
  if (is.null(crs$apriori))
  {
    errorDialog(Rtxt("You first need to generate the association rules.",
                     "Perhaps you need to click the Execute button."))
    return()
  }

  # Note the textview.
  
  TV <- "associate_textview"
  
  # Required information
  
#  lift    <- theWidget("associate_lift_spinbutton")$getValue()

#  appendTextview(TV, "Top Rules\n\n",
#                 "For now, run the following command in the console:\n\n",
#                 paste('inspect(sort(subset(crs$apriori, lift >',
#                       lift, '), by="confidence"))'))

  # I wanted to use the subset function to list just the top rules,
  # but the use of "lift" has issues when run within a library, even
  # though it is fine when loaded using source. The problems seems to
  # have been fixed by 080419 so refine this to do the right
  # thing. Some more testing 080421 indicates I'm back to the old
  # problem...

#  if (lift == 0)
  sby <- ""
  if (! isMac())
      sby <- tolower(theWidget("associate_sort_comboboxtext")$getActiveText())
  sort.ds <- ifelse(sby=="",
                    "crs$apriori",
                    sprintf('sort(crs$apriori, by="%s")', sby))
  summary1.cmd <- sprintf('inspect(%s)', sort.ds)

  
#  else
#lift<-1
#      summary1.cmd <- paste('sort(subset(crs$apriori, lift > ',
#                          lift, '),  by="confidence")')
#    summary1.cmd <- paste('inspect(sort(subset(crs$apriori, lift > ',
#                          lift, '),  by="confidence"))')
  appendLog(Rtxt("List rules."), summary1.cmd)
  # print(summary1.cmd)
  ## This returns "" 080429 when "lift > 1.3" is included in the
  ## subset command.

  ## 111224 collectOutput was failing with:
  ##
  ## Warning message:
  ## In is.na(x) : is.na() applied to non-(list or vector) of type 'S4'
  ##
  ## Error in inspect(sort(crs$apriori, by = "confidence")) : 
  ## error in evaluating the argument 'x' in selecting a method for
  ## function 'inspect': Error in
  ## x[order(x, na.last = na.last, decreasing = decreasing)] : 
  ## error in evaluating the argument 'i' in selecting a method for
  ## function '[': Error in slot(x, s)[i] : subscript out of bounds
  ##
  ## So try capture.output.
  ##
  ## result <- collectOutput(summary1.cmd)
  ##

  result <- paste(capture.output(eval(parse(text=summary1.cmd), envir=globalenv())), collapse="\n")

  ## It did not dump an error, but produces no output until I added
  ## the paste. But then I'm back to this working when I source, but
  ## failing when I run it from the package.
  
  ## This 
##   zz <- textConnection("commandsink", "w", TRUE)
##   sink(zz)
##   cat(eval(parse(text=summary1.cmd)))
##   sink()
##   close(zz)
##   result <- paste(commandsink, collapse="\n")
  # print(result) # DEBUG
  appendTextview(TV, Rtxt("All Rules"), "\n\n", result, "\n")
                 #"\n\nKnown Bug: If nothing appears above, ",
                 #"set the Lift to 0.0\n")
#                 paste('inspect(sort(subset(crs$apriori, lift >',
#                       lift, '), by="confidence"))'))
  
  # This works but it lists all rules.

#  summary.cmd <- 'inspect(sort(crs$apriori, by="confidence"))'
#  appendLog("List all rules.", summary.cmd)
#  appendTextview(TV, "All Rules\n\n", collectOutput(summary.cmd))

  # Measures of interestingness

  measures <- paste('c("chiSquare", "hyperLift", "hyperConfidence",',
                    '"leverage", "oddsRatio", "phi")')
  interest.cmd <- sprintf('interestMeasure(%s, %s, crs$transactions)',
                          sort.ds, measures)
  appendLog(Rtxt("Interesting Measures."), interest.cmd)
  result <- paste(capture.output(eval(parse(text=interest.cmd),
                                      envir=globalenv())), collapse="\n")
  appendTextview(TV, Rtxt("Interestng Measures"), "\n\n", result)

  setStatusBar(Rtxt("Finished listing the rules",
                    "- scroll the text window to view the rules."))
}

########################################################################
#
# EXPORT
#

exportAssociateTab <- function()
{
  # Make sure we have already done something in Rattle.
  
  if (noDatasetLoaded()) return()

  # Make sure we have a model first!
  
  if (is.null(crs$apriori))
  {
    errorDialog(Rtxt("No association rules model is available. Be sure to build",
                     "the model before trying to export it! You will need",
                     "to press the Execute button (F2) in order to build the",
                     "model."))
    return()
  }

  # Require the pmml package
  
  lib.cmd <- "library(pmml, quietly=TRUE)"
  if (! packageIsAvailable("pmml", Rtxt("export association rules"))) return(FALSE)
  appendLog(Rtxt("Load the PMML package to export association rules."), lib.cmd)
  # Load the package unless we already have a pmml defined (through source).
  if (! exists("pmml")) eval(parse(text=lib.cmd))
  
  # Obtain filename to write the PMML to.
  
  dialog <- RGtk2::gtkFileChooserDialog(Rtxt("Export PMML"), NULL, "save",
                                 "gtk-cancel", RGtk2::GtkResponseType["cancel"],
                                 "gtk-save", RGtk2::GtkResponseType["accept"])
  dialog$setDoOverwriteConfirmation(TRUE)

  if(not.null(crs$dataname))
    dialog$setCurrentName(paste(get.stem(crs$dataname), "_arules.xml", sep=""))

  ff <- RGtk2::gtkFileFilterNew()
  ff$setName(Rtxt("PMML Files"))
  ff$addPattern("*.xml")
  dialog$addFilter(ff)

  ff <- RGtk2::gtkFileFilterNew()
  ff$setName(Rtxt("All Files"))
  ff$addPattern("*")
  dialog$addFilter(ff)
  
  if (dialog$run() == RGtk2::GtkResponseType["accept"])
  {
    save.name <- dialog$getFilename()
    dialog$destroy()
  }
  else
  {
    dialog$destroy()
    return()
  }

  # if (get.extension(save.name) == "") save.name <- sprintf("%s.xml", save.name)
    
  # We can't pass "\" in a filename to the parse command in
  # MS/Windows so we have to run the save/write command separately,
  # i.e., not inside the string thaat is being parsed.

  pmml.cmd <- 'pmml(crs$apriori)'
  appendLog(Rtxt("Export association rules as PMML."),
            sprintf('saveXML(%s, "%s")', pmml.cmd, save.name))
  XML::saveXML(eval(parse(text=pmml.cmd)), save.name)

  setStatusBar(sprintf(Rtxt("The PMML file '%s' has been written."), save.name))
}
