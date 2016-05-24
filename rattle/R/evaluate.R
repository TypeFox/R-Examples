# Gnome R Data Miner: GNOME interface to R for Data Mining
#
# Time-stamp: <2016-01-26 10:16:11 gjw>
#
# Implement evaluate functionality.
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

########################################################################
# CALLBACKS

on_evaluate_csv_radiobutton_toggled <- function(button)
{
  if (button$getActive())
    theWidget("evaluate_filechooserbutton")$setSensitive(TRUE)
  else
    theWidget("evaluate_filechooserbutton")$setSensitive(FALSE)
  setStatusBar()
}

on_evaluate_rdataset_radiobutton_toggled <- function(button)
{
  if (button$getActive())
  {
    theWidget("evaluate_rdataset_combobox")$setSensitive(TRUE)
    updateRDatasets(current=theWidget("evaluate_rdataset_combobox")$getActiveText(),
                    cbox.name="evaluate_rdataset_combobox")
  }
  else
    theWidget("evaluate_rdataset_combobox")$setSensitive(FALSE)
  setStatusBar()
}

on_evaluate_confusion_radiobutton_toggled <- function(button)
{
  if (button$getActive())
    crv$EVALUATE$setCurrentPage(crv$EVALUATE.CONFUSION.TAB)
  setStatusBar()
}

on_evaluate_risk_radiobutton_toggled <- function(button)
{
  if (button$getActive())
  {
    crv$EVALUATE$setCurrentPage(crv$EVALUATE.RISK.TAB)
    theWidget("evaluate_risk_variable_label")$setSensitive(TRUE)
    theWidget("evaluate_risk_label")$setSensitive(TRUE)
  }
  else
  {
    theWidget("evaluate_risk_variable_label")$setSensitive(FALSE)
    theWidget("evaluate_risk_label")$setSensitive(FALSE)
  }
  setStatusBar()
}

on_evaluate_lift_radiobutton_toggled <- function(button)
{
  if (button$getActive())
    crv$EVALUATE$setCurrentPage(crv$EVALUATE.LIFT.TAB)
  setStatusBar()
}

on_evaluate_roc_radiobutton_toggled <- function(button)
{
  if (button$getActive())
    crv$EVALUATE$setCurrentPage(crv$EVALUATE.ROC.TAB)
  setStatusBar()
}

on_evaluate_precision_radiobutton_toggled <- function(button)
{
  if (button$getActive())
    crv$EVALUATE$setCurrentPage(crv$EVALUATE.PRECISION.TAB)
  setStatusBar()
}

on_evaluate_sensitivity_radiobutton_toggled <- function(button)
{
  if (button$getActive())
    crv$EVALUATE$setCurrentPage(crv$EVALUATE.SENSITIVITY.TAB)
  setStatusBar()
}

on_evaluate_score_radiobutton_toggled <- function(button)
{
  # Configure the Report and Include options that relate only to
  # scoring. For non-Survival models the report is either Class or
  # Probability. For Survival models (coxph at least) it is Time or
  # Risk. We change the labels appropriately.

  if (button$getActive())
  {
    make.sensitive <- (#length(active.models) == 1 &&
                       !(theWidget("evaluate_kmeans_checkbutton")$getActive() ||
                         theWidget("evaluate_hclust_checkbutton")$getActive())
                       && ! numericTarget())

    # Show the Score textview.

    crv$EVALUATE$setCurrentPage(crv$EVALUATE.SCORE.TAB)

    # Configure the Report/Include options

    theWidget("score_report_label")$setSensitive(make.sensitive)
    theWidget("score_class_radiobutton")$setSensitive(make.sensitive)
    theWidget("score_probability_radiobutton")$setSensitive(make.sensitive)
    if (length(crs$survival) && class(crs$survival) == "survreg")
      theWidget("score_probability_radiobutton")$setSensitive(FALSE)
    theWidget("score_include_label")$setSensitive(TRUE)
    theWidget("score_idents_radiobutton")$setSensitive(TRUE)
    theWidget("score_all_radiobutton")$setSensitive(TRUE)
    theWidget("evaluate_enterdata_radiobutton")$setSensitive(TRUE)

    #091112 I don't think this belongs here anymore.
#    if (not.null(crs$kmeans))
#      theWidget("evaluate_kmeans_checkbutton")$setSensitive(TRUE)
#    if (not.null(crs$hclust))
#      theWidget("evaluate_hclust_checkbutton")$setSensitive(TRUE)
  }
  else
  {
    theWidget("score_report_label")$setSensitive(FALSE)
    theWidget("score_class_radiobutton")$setSensitive(FALSE)
    theWidget("score_probability_radiobutton")$setSensitive(FALSE)
    theWidget("score_include_label")$setSensitive(FALSE)
    theWidget("score_idents_radiobutton")$setSensitive(FALSE)
    theWidget("score_all_radiobutton")$setSensitive(FALSE)
    theWidget("evaluate_enterdata_radiobutton")$setSensitive(FALSE)
  }
  setStatusBar()
}

on_evaluate_pvo_radiobutton_toggled <- function(button)
{
  if (button$getActive())
    crv$EVALUATE$setCurrentPage(crv$EVALUATE.PVO.TAB)
  setStatusBar()
}

on_evaluate_costcurve_radiobutton_toggled <- function(button)
{
  if (button$getActive())
    crv$EVALUATE$setCurrentPage(crv$EVALUATE.COSTCURVE.TAB)
  setStatusBar()
}

on_evaluate_radiobutton_group_changed <- function(action, window)
{
  setStatusBar()
}

on_risk_comboboxentry_changed <- function(action, window)
{
  setStatusBar()
}

########################################################################
# UTILITIES

configureEvaluateTab <- function()
{
  # 091112 Redesign. Configure based on the available models.  Don't
  # change the status of any model checkboxes.  This is the one place
  # where all this is done, and this is then called from the callbacks
  # whenever any of the model buttons within the Evaluate tab is
  # toggled.

  # Check if any models are checked for evaluation.

  active.models <- NULL

  MODEL <- union(crv$PREDICT, setdiff(crv$DESCRIBE, crv$APRIORI))

  for (m in MODEL)
  {
    avail <- ! is.null(eval(parse(text=paste("crs$", m, sep=""))))
    theWidget(paste("evaluate", m, "checkbutton", sep="_"))$setSensitive(avail)
    if (theWidget(paste("evaluate", m, "checkbutton", sep="_"))$getActive())
      active.models <- c(active.models, m)
  }

  active.models.exist <- length(active.models) > 0

  # Automatically work out what needs to be sensistive, based on data
  # type of the target plus whether kmeans or hclust is active and
  # selected.

  TYPE <- c("confusion", "hand", "risk", "costcurve", "lift", "roc",
            "precision", "sensitivity", "pvo", "score")

  if (! active.models.exist)
    buttons <- NULL
  else if (is.null(crs$target))
    buttons <- c("score")
  else if (theWidget("evaluate_kmeans_checkbutton")$getActive() ||
           theWidget("evaluate_hclust_checkbutton")$getActive() ||
           theWidget("evaluate_survival_checkbutton")$getActive()) # 090929 For now.
    buttons <- c("score")
  else if (multinomialTarget())
    buttons <- c("confusion", "score")
  else if (numericTarget())
    buttons <- c("risk", "pvo", "score")
  else if (binomialTarget() && is.factor(crs$dataset[,crs$target]))
    # 090222 If it is a binomial target then do not allow Pr v Ob
    # since the target might be categoric and the display makes no
    # sense (and fails with the jitter function).
    buttons <- setdiff(TYPE, "pvo")
  else
    buttons <- TYPE

  # Now enable each of the identified evaluate buttons.

  for (b in buttons)
    theWidget(paste("evaluate", b, "radiobutton", sep="_"))$setSensitive(TRUE)

  # Disable each of the buttons not specified and make sure none are
  # set as active.

  for (b in setdiff(TYPE, buttons))
  {
    theWidget(paste("evaluate", b, "radiobutton", sep="_"))$setSensitive(FALSE)
    theWidget(paste("evaluate", b, "radiobutton", sep="_"))$setActive(FALSE)
  }

  # Need a button to be set as default. Use the first in the
  # list. 091114 TODO We also need to check if there is a current
  # sensitive button that is active, and if so then leave it as it is.

  sensitive.active <- NULL
  for (b in TYPE)
    if (theWidget(paste("evaluate", b, "radiobutton", sep="_"))$isSensitive()
        && theWidget(paste("evaluate", b, "radiobutton", sep="_"))$getActive())
      sensitive.active <- c(sensitive.active, b)

  if (length(buttons) > 0 && length(sensitive.active) == 0)
    theWidget(paste("evaluate", buttons[1], "radiobutton", sep="_"))$setActive(TRUE)

  # Set the Data options of the Evaluate tab appropraitely. 101116 But
  # only do this if there are active models.

  if (active.models.exist)
  {
    for (b in c("training", "csv", "rdataset"))
      theWidget(paste("evaluate", b, "radiobutton", sep="_"))$setSensitive(TRUE)

    # When we have partitioning enabled, select the appropriate default.

    if (theWidget("data_sample_checkbutton")$getActive())
    {
      theWidget("evaluate_validation_radiobutton")$setSensitive(length(crs$validate))
      # 100328 Until RStat catches up with Rattle's train/validate/test,
      # we need to make the Testing option available for RStat when
      # sampling is active.
      theWidget("evaluate_testing_radiobutton")$setSensitive(length(crs$test) ||
                                                             crv$appname == "RStat")
      theWidget("evaluate_fulldata_radiobutton")$setSensitive(TRUE)
      if (length(crs$validate))
        theWidget("evaluate_validation_radiobutton")$setActive(TRUE)
      else
        theWidget("evaluate_testing_radiobutton")$setActive(TRUE)
    }
    else
    {
      theWidget("evaluate_validation_radiobutton")$setSensitive(FALSE)
      theWidget("evaluate_testing_radiobutton")$setSensitive(FALSE)
      theWidget("evaluate_fulldata_radiobutton")$setSensitive(FALSE)
      theWidget("evaluate_training_radiobutton")$setActive(TRUE)
    }

    # 101116 This is set to FALSE here since it is only available when
    # the Score option is chosen.

    theWidget("evaluate_enterdata_radiobutton")$setSensitive(FALSE)
  }

  #----------------------------------------------------------------------

  # 081206 Handle the sensitivity of the Report options: Class/Time
  # and Probability/Risk. These are only available if one of the
  # non-cluster models is active but not if it is a multinomial
  # target.

  #091114 Now done in the Score callback?
#  predictive.model <- (theWidget("evaluate_rpart_checkbutton")$getActive() ||
#                       theWidget("evaluate_ada_checkbutton")$getActive() ||
#                       theWidget("evaluate_rf_checkbutton")$getActive() ||
#                       theWidget("evaluate_ksvm_checkbutton")$getActive() ||
#                       theWidget("evaluate_glm_checkbutton")$getActive() ||
#                       theWidget("evaluate_nnet_checkbutton")$getActive())
#
#  make.sensitive <- (! is.null(crs$survival) ||
#                     (existsCategoricModel()
#                      && predictive.model
#                      && ! multinomialTarget()))

  make.sensitive <- (#length(active.models) == 1 &&
                     !(theWidget("evaluate_kmeans_checkbutton")$getActive() ||
                       theWidget("evaluate_hclust_checkbutton")$getActive())
                     && ! numericTarget())

  theWidget("score_report_label")$setSensitive(make.sensitive)
  theWidget("score_class_radiobutton")$setSensitive(make.sensitive)
  theWidget("score_probability_radiobutton")$setSensitive(make.sensitive)

  default.to.class <- (theWidget("evaluate_kmeans_checkbutton")$getActive() ||
                       theWidget("evaluate_hclust_checkbutton")$getActive() ||
                       theWidget("evaluate_rpart_checkbutton")$getActive() ||
                       theWidget("evaluate_ada_checkbutton")$getActive() ||
                       theWidget("evaluate_rf_checkbutton")$getActive() ||
                       theWidget("evaluate_ksvm_checkbutton")$getActive() ||
                       (theWidget("evaluate_survival_checkbutton")$getActive() &&
                        class(crs$survival) %in% "survreg"))

  if (default.to.class)
    theWidget("score_class_radiobutton")$setActive(TRUE)
  else
    theWidget("score_probability_radiobutton")$setActive(TRUE)

  # Change the labels for Report if Survival is set.

  if (theWidget("evaluate_survival_checkbutton")$getActive())
  {
    if (length(active.models) == 1) # Only Survival is selected
    {
      theWidget("score_class_radiobutton")$setLabel(Rtxt("Time"))
      theWidget("score_probability_radiobutton")$setLabel(Rtxt("Risk"))
    }
    else
    {
      theWidget("score_class_radiobutton")$setLabel(Rtxt("Class/Time"))
      theWidget("score_probability_radiobutton")$setLabel(Rtxt("Prob/Risk"))
    }
  }
  else if (length(active.models))
  {
    theWidget("score_class_radiobutton")$setLabel(Rtxt("Class"))
    theWidget("score_probability_radiobutton")$setLabel(Rtxt("Probability"))
  }
}

resetEvaluateTab <- function()
{
  # 091112 Simply deactivate and desensitise everything back to
  # default.

  TYPE <- c("confusion", "hand", "risk", "costcurve", "lift", "roc",
            "precision", "sensitivity", "pvo", "score")

  for (b in TYPE)
    theWidget(paste("evaluate", b, "radiobutton", sep="_"))$setSensitive(FALSE)
  theWidget("evaluate_confusion_radiobutton")$setActive(TRUE)

  DATA <- c("training", "validation", "testing", "csv", "rdataset", "fulldata",
            "enterdata")

  for (b in DATA)
  {
    theWidget(paste("evaluate", b, "radiobutton", sep="_"))$setSensitive(FALSE)
    theWidget(paste("evaluate", b, "radiobutton", sep="_"))$setActive(FALSE)
  }

  MODEL <- union(crv$PREDICT, setdiff(crv$DESCRIBE, crv$APRIORI))

  for (m in MODEL)
  {
    theWidget(paste("evaluate", m, "checkbutton", sep="_"))$setSensitive(FALSE)
    theWidget(paste("evaluate", m, "checkbutton", sep="_"))$setActive(FALSE)
  }

  # Scoring options

  theWidget("score_report_label")$setSensitive(FALSE)
  theWidget("score_class_radiobutton")$setSensitive(FALSE)
  theWidget("score_probability_radiobutton")$setSensitive(FALSE)
  theWidget("score_include_label")$setSensitive(FALSE)
  theWidget("score_idents_radiobutton")$setSensitive(FALSE)
  theWidget("score_all_radiobutton")$setSensitive(FALSE)

  # 091215 Set default choice back to Confusion Matrix

  theWidget("evaluate_confusion_radiobutton")$setActive(TRUE)

}

getEvaluateModels <- function()
{
  # Return a list of models selected for evaluation

  models <- c()

  # First, check each of the traditional model checkboxes.

  for (m in crv$PREDICT)
    if (theWidget(paste("evaluate", m, "checkbutton", sep="_"))$getActive())
      models <- c(models, m)

  # Now add in the cluster models, which will eventually be part of
  # the Model tab.

  if (theWidget("evaluate_kmeans_checkbutton")$isSensitive() &&
      theWidget("evaluate_kmeans_checkbutton")$getActive())
    models <- c(models, "kmeans")

  if (theWidget("evaluate_hclust_checkbutton")$isSensitive() &&
      theWidget("evaluate_hclust_checkbutton")$getActive())
    models <- c(models, "hclust")

  return(models)
}

current.evaluate.tab <- function()
{
  cp <- crv$EVALUATE$getCurrentPage()
  return(crv$EVALUATE$getTabLabelText(crv$EVALUATE$getNthPage(cp)))
}

########################################################################
# EXECUTION

executeEvaluateTab <- function()
{
  # Perform the requested action from the Execute tab.

  # Obtain some background information.

  mtypes <- getEvaluateModels() # The chosen model types in the Evaluate tab.

  # Check any pre-conditions.

  # Ensure a dataset exists.

  if (noDatasetLoaded()) return()

  # Ensure we have at least one model to evaluate, otherwise warn the
  # user and do nothing.

  if (is.null(listBuiltModels(crs$APRIORI)))
  {
    warnDialog(Rtxt("No models suitable for evaluation have been built.\n\n",
                    "Please build a suitable model before evaluation."))
    return()
  }

  if (is.null(mtypes))
  {
    warnDialog(Rtxt("No model has been specified.",
                    "\n\nPlease select one or more from the",
                    "list of models available."))
    return()
  }

  # Ensure we recognise the model type.

  if (length(setdiff(mtypes, union(crv$PREDICT, c("kmeans", "hclust")))) > 0)
  {
    errorDialog(sprintf(Rtxt("E121: A model type is not recognised.",
                             "We found the model types to be: %s Known models: %s"),
                        mtypes, crv$PREDICT),
                "\n", crv$support.msg)
    return()
  }

  # Ensure there is a model for each model type that is selected.

  if (sum(sapply(mtypes, function(x) is.null(crs[[x]]))) > 0)
  {
    errorDialog(sprintf(Rtxt("E120: Some model has not been built?",
                             "We found the model types to be: %s",
                             "The models not built: %s",
                             "This is a Rattle bug."),
                        mtypes, sapply(mtypes, function(x) is.null(crs[[x]]))),
                "\n", crv$support.msg)
    return()
  }

  #   Ensure the appropriate package is loaded (in the case, for
  #   example, when loading a project and going straight to Evaluate,
  #   and wanting to run predict.svm on new data).

  if (crv$ADA %in%  mtypes &&
      ! packageIsAvailable("ada", sprintf(Rtxt("evaluate a %s model"),
                                          commonName(crv$ADA))))
    return()
  if (crv$KSVM %in%  mtypes &&
      ! packageIsAvailable("kernlab", sprintf(Rtxt("evaluate a %s model"),
                                              commonName(crv$KSVM))))
    return()
  if (crv$RF %in%  mtypes &&
      ! packageIsAvailable("randomForest", sprintf(Rtxt("evaluate a %s model"),
                                                   commonName(crv$RF))))
    return()
  if (crv$GLM %in%  mtypes && "multinom" %in% class(crs$glm) &&
      ! packageIsAvailable("nnet", sprintf(Rtxt("evaluate a %s model"),
                                           paste(" Multinomial", commonName(crv$GLM)))))
    return()
  if (crv$NNET %in% mtypes &&
      ! packageIsAvailable("nnet", sprintf(Rtxt("evaluate a %s model"),
                                           commonName(crv$NNET))))
    return()
  if (crv$SURVIVAL %in% mtypes &&
      ! packageIsAvailable("survival", sprintf(Rtxt("evaluate a %s model"),
                                               commonName(crv$SURVIVAL))))
    return()

  if(theWidget("evaluate_score_radiobutton")$getActive())
    startLog(Rtxt("Score a dataset."))
  else
    startLog(Rtxt("Evaluate model performance."))

  # Identify the data on which evaluation is to be performed.

  testset0 <- "crs$dataset"
  testname <- crs$dataname

  # 081028 For included we only need the input variables and perhaps
  # the risk variable. But after changing the definition of the
  # arguments to getIncludedVariables, where risk=FALSE by default, I
  # forgot to set it to TRUE here. However, it seems to be working so
  # far, at least for glm! 081029 However, we need the target variable
  # in the list for confusion matrix and risk chart, for
  # example. 100530 But we don't need the target for scoring, and so
  # we should remove it if we are scoring.

  #included <- getIncludedVariables(target=FALSE)
  if (theWidget("evaluate_score_radiobutton")$getActive())
    included <- "c(crs$input)" # 20110102 getIncludedVariables(target=FALSE)
  else
    included <- "c(crs$input, crs$target)" # 20110102 getIncludedVariables()

  if (theWidget("evaluate_training_radiobutton")$getActive())
  {
    # Evaluate on training data

    if (crv$show.warnings && theWidget("data_sample_checkbutton")$getActive())
      infoDialog(Rtxt("You are using the training dataset to evaluate your model.",
                      "This will give you an optimistic estimate",
                      "of the performance of your model.",
                      "\n\nYou may want to choose",
                      "to sample the dataset and evaluate the model on the",
                      "test dataset, or else",
                      "load a separate test dataset from a CSV File or a",
                      "pre-existing R Dataset here."))

    if (theWidget("data_sample_checkbutton")$getActive())
      if (is.null(included))
        testset0 <- "crs$dataset[crs$sample,]"
      else
        testset0 <- sprintf("crs$dataset[crs$sample, %s]", included)
    else
      if (is.null(included))
        testset0 <- "crs$dataset"
      else
        testset0 <- sprintf("crs$dataset[,%s]", included)

    testname <- sprintf("%s [**%s**]", crs$dataname, Rtxt("train"))
  }
  else if (theWidget("evaluate_validation_radiobutton")$getActive())
  {
    # Evaluate on validation data

    if (is.null(included))
      testset0 <- "crs$dataset[crs$validate,]"
    else
      testset0 <- sprintf("crs$dataset[crs$validate, %s]", included)
    testname <- sprintf("%s [%s]", crs$dataname, Rtxt("validate"))
  }
  else if (theWidget("evaluate_testing_radiobutton")$getActive())
  {
    # Evaluate on test data

    if (is.null(included))
      if (newSampling())
        testset0 <- "crs$dataset[crs$test,]"
      else
        testset0 <- "crs$dataset[-crs$sample,]"
    else
      if (newSampling())
        testset0 <- sprintf("crs$dataset[crs$test, %s]", included)
      else
        testset0 <- sprintf("crs$dataset[-crs$sample, %s]", included)
    testname <- sprintf("%s [%s]", crs$dataname, Rtxt("test"))
  }
  else if (theWidget("evaluate_csv_radiobutton")$getActive())
  {
    # Evaluate on CSV or TXT data. We identify which from the
    # seperator specified from the Data tab.

    # We need to allow for the case where the loaded csv data does not
    # have the risk and target variables when we are scoring the data
    # (i.e., not when we are generating confusion charts and other
    # evaluations). For scoring, it is only natural that we do not
    # have the risk and target variables.

    filename <- theWidget("evaluate_filechooserbutton")$getFilename()
    crs$dwd <- ifelse(length(filename), dirname(filename), "")
    crs$mtime <- urlModTime(filename)

    if (is.null(filename)
        || ! file.exists(filename)
        || (! isWindows() && file.info(filename)$isdir))
    {
      errorDialog(Rtxt("You have requested that a CSV file be used",
                       "as your testing dataset, but you have not",
                       "identified which file.\n\nPlease use the Spreadsheet",
                       "button to select the CSV file you wish",
                       "to use as your testset before you Execute."))
      return()
    }

    # Load the testset from file, but only load it if it is not
    # already loaded.

    if (is.null(crs$testset)
        || is.null(crs$testname)
        || (basename(filename) != crs$testname))
    {
      # Fix filename for MS/Windows - otherwise eval/parse strips the \\.

      if (isWindows()) filename <- gsub("\\\\", "/", filename)

      nastring <- ', na.strings=c(".", "NA", "", "?")'
      sep = theWidget("data_separator_entry")$getText()
      hdr = theWidget("data_header_checkbutton")$getActive()
      read.cmd <- sprintf(paste('crs$testset <- read.csv("%s"%s, header=%s,',
                                'sep="%s", encoding="%s", strip.white=TRUE)'),
                          filename, nastring,
                          ifelse(hdr, "TRUE", "FALSE"),
                          sep, crv$csv.encoding)

      appendLog(Rtxt("Read a dataset from file for testing the model."), read.cmd)
      eval(parse(text=read.cmd))

      testname <- basename(filename)
      crs$testname <- testname
    }

    # TODO The following case for included assumes the same column
    # orders. Should really check this to make sure.  For scoring a
    # dataset we do not include the target or the risk in the
    # variables, since they may not be present in the csv file that is
    # being loaded (if that option is active). Thus, in this case it
    # is best to simply use the whole dataset for scoring. But, for
    # the case where there are lots of columns that are ignored in the
    # model building, if they have lots of NAs then the scoring is
    # going to give NAs for RF, etc. (Pointed out by Ed Cox 9 Feb
    # 2008.) In general, not sure how to handle this, except for now
    # say that the schema must be identical in the scoring dataset to
    # the training dataset (including the target, risk, and ignored
    # columns). In fact, if the target etc are the last columns then
    # we can get away with it.

    if (is.null(included)) # || theWidget("score_radiobutton")$getActive())
      testset0 <- "crs$testset"
    else
      testset0 <- sprintf("crs$testset[,%s]", included)
  }
  else if (theWidget("evaluate_rdataset_radiobutton")$getActive())
  {
    dataset <- theWidget("evaluate_rdataset_combobox")$
               getActiveText()

    if (is.null(dataset) || nchar(dataset) == 0)
    {
      errorDialog(Rtxt("The R Dataset is active but",
                       "no dataset name has been specified.",
                       "Please identify the name of the R dataset",
                       "on which you would like to evaluate the model.",
                       "This dataset will be one that has been defined",
                       "in the R Console."))
      return()
    }

    testset0 <- 'crs$testset'
    testname <- dataset
    crs$testname <- testname

    assign.cmd <- sprintf("crs$testset <- %s", dataset)
    appendLog(Rtxt("Assign the R dataset to be used as the test set."), assign.cmd)
    eval(parse(text=assign.cmd))
  }

  # Ensure the test dataset has the same levels for each variable of
  # the training dataset. This can arise when we externally split a
  # dataset into a training and testing dataset, and the smaller
  # testing dataset may not have representatives of all of the
  # variables. Be sure to add any new levels to the end, otherwise
  # you'll end up renaming some of the other levels! This won't help a
  # model that uses the variable and does not find the particular
  # level, although it is okay if it is missing levels. TODO this
  # might need to check for the former and error out if it is the
  # case. TODO 081204 I don't need to do this for every factor, just
  # those with different levels.

  if (not.null(crs$testname) && crs$testname != crs$dataname)
    # 130128 Why do we need to add the target variable back in?
    # Especially when adding it in and setting it to NA, na.omit for
    # RF has a particular issue!!!! So don't add target back in if it
    # is missing. TODO The code around here has evolved poorly and
    # could really do with quite a cleanup.
    for (cn in setdiff(colnames(crs$dataset), crs$target))
      if (is.factor(crs$dataset[[cn]]))
      {
        # 100701 If the categoric variable is missing (like it might
        # be for the target variable) then be sure to add it in. TODO
        # 100701 Does the order of the variables matter in scoring? I
        # think it does (since we use indicies to subset), so be sure
        # to insert the column in the right place.

        if (! cn %in% colnames(crs$testset))
        {
          place <- which(cn == colnames(crs$dataset))
          cmd <- sprintf(paste("crs$testset <- cbind(crs$testset[1:%d], ",
                               "%s=rep(NA, nrow(crs$testset))",
                               ifelse(place < ncol(crs$testset),
                                      ", crs$testset[%d:ncol(crs$testset)]",
                                      "%s"),
                               ")", sep=""),
                         place-1, cn,
                         ifelse(place < ncol(crs$testset), place, ""))
          appendLog(sprintf(Rtxt("Add missing column `%s'",
                                 "to the testing dataset"), cn),
                          cmd)
          eval(parse(text=cmd))
        }

        # 090808 Be sure to expose this trick to the log file since
        # the user will otherwise be unable to repeat the scoring for
        # the case where the levels are not the same as the training
        # dataset.

        cmd <- sprintf(paste('levels(crs$testset[["%s"]]) <-',
                             '\n  c(levels(crs$testset[["%s"]]),',
                             '\n    setdiff(levels(crs$dataset[["%s"]]),',
                             '\n               levels(crs$testset[["%s"]])))'),
                       cn, cn, cn, cn)
        appendLog(sprintf(Rtxt("Ensure the levels are the same as the training",
                               "data for variable `%s'."), cn),
                          cmd)
        eval(parse(text=cmd))
      }

  ## The default command for prediction from any model is
  ## predict(model, data). Here we tune the predict command to
  ## particular types of models where they have specific general
  ## requirements. We then modify the default predict command to
  ## generate either a prediction of the response or a probability of
  ## the class, as appropriate to the particular evaluator.
  ##
  ## PREDICT: crs$pr <- predict(crs$model, crs$testset[crs$sample, c(...)])

  ## PROBABILITY: this predicts a matrix, each column a probability
  ## for that class.

  ## We want to obtain the probablity of class 1 (i.e., the second of
  ## a two level class). Start with the default predict.cmd.

  # Now build model specific strings for each model

  testset <- list() # The string representing the test dataset
  predcmd <- list() # Command string for predictions
  respcmd <- list() # Command string for response - class of entities
  probcmd <- list() # Command string for probability

  # Why a Predict and Response command? Need better documentation.
  #
  # modeller    pred		resp		prob
  # ada				=pred
  # kmeans			=pred		=pred
  # hclust			=pred		=pred
  # survival			=pred

  if (crv$ADA %in%  mtypes)
  {
    testset[[crv$ADA]] <- testset0

    predcmd[[crv$ADA]] <- genPredictAda(testset[[crv$ADA]])
    respcmd[[crv$ADA]] <- genResponseAda(testset[[crv$ADA]])
    probcmd[[crv$ADA]] <- genProbabilityAda(testset[[crv$ADA]])
  }

  if (crv$KMEANS %in% mtypes)
  {
    testset[[crv$KMEANS]] <- testset0

    # These are all the same!

    predcmd[[crv$KMEANS]] <- genPredictKmeans(testset[[crv$KMEANS]])
    respcmd[[crv$KMEANS]] <- genResponseKmeans(testset[[crv$KMEANS]])
    probcmd[[crv$KMEANS]] <- genProbabilityKmeans(testset[[crv$KMEANS]])
  }

  if (crv$HCLUST %in% mtypes)
  {
    testset[[crv$HCLUST]] <- testset0

    # These are all the same!

    predcmd[[crv$HCLUST]] <- genPredictHclust(testset[[crv$HCLUST]])
    respcmd[[crv$HCLUST]] <- genResponseHclust(testset[[crv$HCLUST]])
    probcmd[[crv$HCLUST]] <- genProbabilityHclust(testset[[crv$HCLUST]])
  }

  if (crv$SURVIVAL %in%  mtypes)
  {
    testset[[crv$SURVIVAL]] <- sprintf("na.omit(%s)", testset0)

    predcmd[[crv$SURVIVAL]] <- genPredictSurvival(testset[[crv$SURVIVAL]])
    respcmd[[crv$SURVIVAL]] <- genResponseSurvival(testset[[crv$SURVIVAL]])
    probcmd[[crv$SURVIVAL]] <- genProbabilitySurvival(testset[[crv$SURVIVAL]])
  }

  if (crv$NNET %in%  mtypes)
  {
    testset[[crv$NNET]] <- testset0

    # 090316 For a binomial target the output node is built with
    # linout=TRUE, thus the value is trying to be close to 1 or 0. Use
    # 0.5 threshold to predict it as 1 or 0. 090820 Move to doing a
    # logistic for nnet, so that we don't use linout in the model
    # building.

    predcmd[[crv$NNET]] <- sprintf("crs$pr <- predict(crs$nnet, newdata=%s)",
                                   testset[[crv$NNET]])

    if (binomialTarget())
      respcmd[[crv$NNET]] <- sub(")$", ', type="class")', predcmd[[crv$NNET]])

#090820 REMOVE the old commented code from here once nnet is stable.
#      predcmd[[crv$NNET]] <- sprintf("crs$pr <- as.integer(predict(crs$nnet, %s)>=0.5)",
#                                     testset[[crv$NNET]])

    else
      respcmd[[crv$NNET]] <- predcmd[[crv$NNET]]

    if (binomialTarget())
#      # 090809 Change from using the pred command to using the raw
#      # predicticed value, which is what it originally probably meant
#      # to be.  090820 TODO Combine these two that both set to
#      # basecmd. In fact, currently nnet is either only binomial or
#      # numeric?
      probcmd[[crv$NNET]] <- predcmd[[crv$NNET]]
#      ## gsub(")$", ', type="raw")', predcmd[[crv$NNET]])
    else if (numericTarget())
      probcmd[[crv$NNET]] <- predcmd[[crv$NNET]]
    else
      probcmd[[crv$NNET]] <- gsub(")$", ', type="prob")', predcmd[[crv$NNET]])
  }

  if (crv$RPART %in%  mtypes)
  {
    cond.tree <- attr(class(crs$rpart), "package") %in% "party"
    if (! length(cond.tree)) cond.tree <- FALSE
    
    testset[[crv$RPART]] <- testset0
    predcmd[[crv$RPART]] <- sprintf("crs$pr <- predict(crs$rpart, newdata=%s)",
                                testset[[crv$RPART]])

    # For crv$RPART, the default is to generate class probabilities for
    # each output class, so ensure we instead generate the response.

    respcmd[[crv$RPART]] <- gsub(")$",
                                 ifelse(cond.tree,
                                        ', type="response")',
                                        ', type="class")'),
                                 predcmd[[crv$RPART]])

    # For RPART the default predict command generates the probabilities
    # for each class and we assume we are interested in the final class
    # (i.e., for binary classification we are interested in the 1's).

    if (cond.tree)
      probcmd[[crv$RPART]] <- sub(')$', '), function(x) x[2])',
                                  sub("predict", "sapply(treeresponse",
                                      predcmd[[crv$RPART]]))
    else
      if (binomialTarget())
        probcmd[[crv$RPART]] <- sprintf("%s[,2]", predcmd[[crv$RPART]])
      else
        probcmd[[crv$RPART]] <- sprintf("%s", predcmd[[crv$RPART]])

    if (multinomialTarget())
    {
      # 081226 Add on the actual class also. This is useful for Score
      # but may be a problem for other types of evaluations (of which
      # there are currently none that use probcmd for multinom).

      probcmd[[crv$RPART]] <- sub("<- ", "<- data.frame(",
                               sub(")$",
                                   sprintf(paste("), rpart=predict(crs$rpart,",
                                                 "newdata=%s, type='class'))"),
                                           testset[[crv$RPART]]),
                                   probcmd[[crv$RPART]]))
    }
  }

  if (crv$RF %in%  mtypes)
  {
    # 090301 Having added support for random forest regression seems
    # like we need to take into acocunt missing for PvO and scoring
    # with numeric targets. In fact, we can probably add na.omit also
    # for categoric targets, since randomForest also does na.omit
    # internally. So it won't help, and will keep in line with other
    # algorithms that actually need the na.omit to be done here.

    cond.rf <- "RandomForest" %in% class(crs$rf) # party conditional rf

    # 090301 testset[[crv$RF]] <- testset0
    if (cond.rf)
      testset[[crv$RF]] <- testset0 # 130323 cforest handles NA in predict
    else
      testset[[crv$RF]] <- sprintf("na.omit(%s)", testset0)

    predcmd[[crv$RF]] <- sprintf("crs$pr <- predict(crs$rf, newdata=%s)",
                             testset[[crv$RF]])

    # The default for crv$RF is to predict the class, so no
    # modification of the predict command is required.

    respcmd[[crv$RF]] <- predcmd[[crv$RF]]

    # For RF we request a probability with the type argument, and as
    # with RPART we extract the column of interest (the last column).

    if (numericTarget())
      probcmd[[crv$RF]] <- predcmd[[crv$RF]]
    else
      if (cond.rf) 
        probcmd[[crv$RF]] <- sub(')$', '), function(x) x[2])',
                                 sub("predict", "sapply(treeresponse",
                                     predcmd[[crv$RF]]))
      else
        probcmd[[crv$RF]] <- sprintf("%s[,2]",
                                     gsub(")$", ', type="prob")', predcmd[[crv$RF]]))

  }

  if (crv$KSVM %in%  mtypes)
  {

    ## For SVM and KSVM, we need to deal with NA's. The predict seems to
    ## need to have NA's removed from the testset, (unlike rpart and rf
    ## where perhaps the NAs don't appear in columns that are used in
    ## the model? An SVM will use all the columns. But in the way we
    ## construct the evaluate command we add extra columns in the third
    ## argument to make sure we get the risk variable in the dataset.
    ## So we need to ensure we get the same subset. It might be smaller
    ## otherwise since the extra columns may have NAs.
    ##
    ## 060527 Comment this out since I was ending up with different
    ## lengths in the 2nd and 3rd arguments in the call to evaluateRisk
    ## in the svm stuff? Using survery-2k the 2nd arg length was 600 and
    ## the 3rd 561, both using na.omit. Perhaps the 3rd should not be
    ## using na.omit, but I haven't investigated this.
    ##
    ## 060603 Put this back in!!! I was again getting the 600 in the
    ## testdata and 561 in the result from predict. Doing a
    ## na.omit(testset) resulted in 561, from the original 600., which
    ## matches the number output from the predict. Try out
    ## survey-training with 10% training to see that that also works! I
    ## suspect that if it does not work, then the issue is missing
    ## levels.
    ## 060622 Seems like the problem is that the na.omit is working on
    ## different subsets of the columns:
    ##   na.omit(crs$dataset[-crs$sample, c(2:22,25)]) versus
    ##   na.omit(crs$dataset[-crs$sample,])
    ## because in the second one we want to retrieve the Risk variable,
    ## which is
    ## not in the first! Instead, let's always extract the list of omitted
    ## rows, and use that here.
    ##
    ##  romit <- attr(na.omit(testset), "na.action")[]
    ## Then use testset[-romit,]
    ## Note that predict automatically removes NAs.
    ## I.e.
    ## crs$eval <- evaluateRisk(crs$pr, na.omit(crs$dataset[-crs$sample,
    ## c(2:22,25)])$Target1,crs$dataset[-crs$sample,][-romit,]$NETADJ_AS_LBLTY)
    ##
    ## 060623 For now in the risk chart function we add the risk
    ## variable back into the testset to ensure it can be accessed,
    ## whereas previously we added in all columns to ensure the risk
    ## variable was included, and this latter resulted in much more
    ## potential for a row to include an NA, and hence to be omitted,
    ## leading to different sized vetors being passed to evaluateRisk. I
    ## hope this now solves the problem and we don't need the top
    ## solution for now.

    testset[[crv$KSVM]] <- sprintf("na.omit(%s)", testset0)

    predcmd[[crv$KSVM]] <- sprintf("crs$pr <- kernlab::predict(crs$ksvm, newdata=%s)",
                               testset[[crv$KSVM]])

    ## The default for KSVM is to predict the class, so no
    ## modification of the predict command is required.

    respcmd[[crv$KSVM]] <- predcmd[[crv$KSVM]]

    ## For KSVM we request a probability with the type argument set to
    ## probability (but need prob.model=TRUE in model building). For SVM
    ## we request probabilities by setting probablity=TRUE and don't
    ## need the second column stuff (and in building the model we needed
    ## probability=TRUE).

    probcmd[[crv$KSVM]] <- sprintf("%s[,2]",
                               gsub(")$",
                                    ', type="probabilities")',
                                    predcmd[[crv$KSVM]]))
    ## For SVM:
    ## probability.cmd <- sprintf("%s",
    ##                             gsub(")$",
    ##                                  ', probability=TRUE)',
    ##                                  probability.cmd))
  }

  if (crv$GLM %in%  mtypes)
  {
    # 080716 The multinom model has been moved to GLM, even though it
    # is using the nnet library. So we need to do the nnet predict
    # here.

    if ("multinom" %in% class(crs$glm))
    {
      testset[[crv$GLM]] <- testset0
      predcmd[[crv$GLM]] <- sprintf("crs$pr <- predict(crs$glm, newdata=%s)",
                                     testset[[crv$GLM]])
      respcmd[[crv$GLM]] <- predcmd[[crv$GLM]]
      probcmd[[crv$GLM]] <- sub(")$", ', type="prob")', predcmd[[crv$GLM]])

      # Add on the actual class also. This is useful for Score but may
      # be a problem for other types of evaluations (of which there
      # are currently none that use probcmd for multinom).

      probcmd[[crv$GLM]] <- sub("<- ", "<- cbind(",
                                sub(")$",
                                    sprintf(paste("), crs$glm$lab[predict(crs$glm,",
                                                  "newdata=%s)])"),
                                            testset[[crv$GLM]]),
                                    probcmd[[crv$GLM]]))

    }
    else
    {

      # GLM's predict removes rows with missing values, so we also need
      # to ensure we remove rows with missing values here.

      # 081029 Try without na.omit since if the target has missing
      # values the record won't be scored, yet there is no reason not
      # to score it. Example is w_reg_logistic.

      ## testset[[crv$GLM]] <- sprintf("na.omit(%s)", testset0)

      testset[[crv$GLM]] <- testset0

      predcmd[[crv$GLM]] <- sprintf(paste("crs$pr <- predict(crs$glm,",
                                          'type="response", newdata=%s)'),
                                    testset[[crv$GLM]])

      # For GLM, a response is a figure close to the class, either close
      # to 1 or close to 0, so threshold it to be either 1 or 0. TODO
      # Simplify this like?
      #    response.cmd <- gsub("predict", "(predict",
      #                         gsub(")$", ")>0.5)*1", response.cmd))

      # 081025 Why do the as.factor? Try just the 0/1 instead. In fact
      # we have now modified this to use the actual levels.

##      respcmd[[crv$GLM]] <- gsub("predict", "as.factor(as.vector(ifelse(predict",
##                                  gsub(")$", ', type="response") > 0.5, 1, 0)))',
##                                       predcmd[[crv$GLM]]))

      # 081029 No longer need response - already there?
      # lvls <- sprintf(', type="response") > 0.5, "%s", "%s"))',
      lvls <- sprintf(') > 0.5, "%s", "%s"))',
                      levels(as.factor(crs$dataset[[crs$target]]))[2],
                      levels(as.factor(crs$dataset[[crs$target]]))[1])
      respcmd[[crv$GLM]] <- gsub("predict", "as.vector(ifelse(predict",
                                 gsub(")$", lvls, predcmd[[crv$GLM]]))

      # For GLM, the response is a probability of the class.

      # 081029 No longer need response - already there?
      # probcmd[[crv$GLM]] <- gsub(")$", ', type="response")', predcmd[[crv$GLM]])
      probcmd[[crv$GLM]] <- predcmd[[crv$GLM]]
    }
  }

##   if (GBM %in%  mtypes)
##   {
##     testset[[GBM]] <- testset0

##     ## For GBM the default needs to know the number of trees to include.

##     predcmd[[GBM]] <- sprintf(paste("crs$pr <- predict(crs$gbm, %s,",
##                                     "n.trees=length(crs$gbm$trees))"),
##                               testset[[GBM]])
##     respcmd[[GBM]] <- predcmd[[GBM]]
##     probcmd[[GBM]] <- predcmd[[GBM]]
##   }

  # Currently (and perhaps permanently) the ROCR package deals only
  # with binary classification, as does my own Risk Chart.

  if (!(theWidget("evaluate_confusion_radiobutton")$getActive()

        # 090506 I had a note here that pvo was not working for
        # multiclass targets for the PrvOb plot, but uncommenting the
        # following and having a numeric target from a categoric
        # variable we get a decent looking plot, but we get a warning
        # about only fitting to first two points and the plotted lines
        # may be relating only to the first two etc. So until I review
        # what happens leave this out. Make sure pvo is not available
        # under such circumstances. Add the check of the Categoric
        # target choice to cover the case of a Numeric override of a
        # Categoric variable.

        || theWidget("evaluate_pvo_radiobutton")$getActive()
        || theWidget("evaluate_score_radiobutton")$getActive())
      && !theWidget("data_target_categoric_radiobutton")$getActive() # See Note Above
      && is.factor(crs$dataset[[crs$target]])
      && length(levels(crs$dataset[[crs$target]])) > 2)
  {
    errorDialog(Rtxt("The number of levels in the target variable is greater",
                     "than 2. Currently, Risk charts and the ROCR package",
                     "(which implements the Lift, ROC, Precision, and Specificity",
                     "charts) apply only to binary classification."))
    return()
  }

  # Dispatch to the appropriate function.

  if (theWidget("evaluate_confusion_radiobutton")$getActive())
    msg <- executeEvaluateConfusion(respcmd, testset, testname)
  else if (theWidget("evaluate_risk_radiobutton")$getActive())
    msg <- executeEvaluateRisk(probcmd, testset, testname)
  else if (theWidget("evaluate_costcurve_radiobutton")$getActive())
    msg <- executeEvaluateCostCurve(probcmd, testset, testname)
  else if (theWidget("evaluate_roc_radiobutton")$getActive())
    msg <- executeEvaluateROC(probcmd, testset, testname)
  else if (theWidget("evaluate_lift_radiobutton")$getActive())
    msg <- executeEvaluateLift(probcmd, testset, testname)
  else if (theWidget("evaluate_precision_radiobutton")$getActive())
    msg <- executeEvaluatePrecision(probcmd, testset, testname)
  else if (theWidget("evaluate_sensitivity_radiobutton")$getActive())
    msg <- executeEvaluateSensitivity(probcmd, testset, testname)
  else if (theWidget("evaluate_hand_radiobutton")$getActive())
    msg <- executeEvaluateHand(probcmd, testset, testname)
  else if (theWidget("evaluate_pvo_radiobutton")$getActive())
  {
    if (categoricTarget())
      msg <- executeEvaluatePvOplot(probcmd, testset, testname)
    else if (numericTarget())
      msg <- executeEvaluatePvOplot(predcmd, testset, testname)
  }

  else if (theWidget("evaluate_score_radiobutton")$getActive())
  {
    if (categoricTarget() || crv$SURVIVAL %in%  mtypes)

      # 081025 Which is best? For trees, traditionally we return the
      # class, but for logistic regression we might return the
      # probability. 081204 So we pass both to the function and decide in
      # there based on a radiobutton setting.

      # 091115 Add the survival option. For a survival model we want
      # access to both the prob and resp (i.e., pred) commands, yet we
      # do not have a categoric target. The prob will be the risk
      # prediction and the resp will be the time to event
      # prediction. Unfortunately, if other models are also selected,
      # this will fail since normally they would go down the other
      # branch here - should be testing that case also I think and put
      # up a warning.

      msg <- executeEvaluateScore(probcmd, respcmd, testset, testname)

    else

      msg <- executeEvaluateScore(predcmd, predcmd, testset, testname)
  }
  else
    msg <- Rtxt("No appropriate evaluator found.")

  if (not.null(msg)) setStatusBar(msg)
}

#----------------------------------------------------------------------
# Evaluate Confusion Table

executeEvaluateConfusion <- function(respcmd, testset, testname)
{
  TV <- "confusion_textview"

  resetTextview(TV)

  for (mtype in getEvaluateModels())
  {

    setStatusBar(sprintf(Rtxt("Applying the %s model to the dataset",
                              "to generate a confusion matrix..."),
                         commonName(mtype)))

    # Generate the command to show the confusion matrix.

    # To ensure we get the target variable in the list, remove any
    # column limits. Doesnot work for na.omit.....

    #ts <- sub(',.*\\]', ', ]', testset[[mtype]])
    ts <- testset[[mtype]]

    confuse.cmd <- paste(sprintf("table(%s$%s, crs$pr,",
                                 ts, crs$target),
                         '\n        useNA="ifany",',
                         '\n        dnn=c("',
                         Rtxt("Actual"), '", "',
                         Rtxt("Predicted"), '"))', sep="")

    percentage.cmd <- paste('pcme <- function(actual, cl)',
                            '{',
                            '  x <- table(actual, cl)',
                            '  nc <- nrow(x) # Number of classes.',
                            '  nv <- length(actual) - sum(is.na(actual) | is.na(cl)) # Number of values.',
                            '  tbl <- cbind(x/nv,',
                            '               Error=sapply(1:nc,',
                            '                 function(r) round(sum(x[r,-r])/sum(x[r,]), 2)))',
                            '  names(attr(tbl, "dimnames")) <- c("Actual", "Predicted")',
                            '  return(tbl)',
                            '}',
                            sprintf('per <- pcme(%s$%s, crs$pr)', ts, crs$target),
                            'round(per, 2)',
                            sep="\n")

    if (categoricTarget())
    {
      # 080528 TODO generalise to categoricTarget. 091023 Handle the
      # case where there is only one value predicted from the two
      # possible values. 150715 TODO Need to update the formulation
      # to account for a single row in the percentage table.

      error.cmd <- 'cat(100*round(1-sum(diag(per), na.rm=TRUE), 2))'
      avgerr.cmd <- 'cat(100*round(mean(per[,"Error"], na.rm=TRUE), 2))'
      
    }

    # Log the R commands and execute them.

    appendLog(sprintf(Rtxt("Generate an Error Matrix for the %s model."),
                      commonName(mtype)))
    appendLog(sprintf(Rtxt("Obtain the response from the %s model."),
                      commonName(mtype)), respcmd[[mtype]])

    result <- try(eval(parse(text=respcmd[[mtype]])), TRUE)

    # Check for errors - in particular, new levels in the test dataset.

    if (inherits(result, "try-error"))
    {
      if (any(grep("has new level", result)) || any(grep("New levels",result)))
        infoDialog(sprintf(Rtxt("It seems that the dataset on which the predictions",
                                "from the %s model are required has a categoric",
                                "variable with levels not found in the training",
                                "dataset. The predictions can not be made in",
                                "this situation. You may need to either ensure",
                                "the training dataset has representatives of all levels",
                                "or else remove them from the testing dataset.",
                                "Alternatively, do not include that variable in the",
                                "modelling. \n\n The actual error message was:\n\n%s"),
                           mtype, paste(result, "\n")))
      else if (any(grep("undefined columns", result)))
        infoDialog(sprintf(Rtxt("It seems that the dataset on which the predictions",
                                "from the %s model are required has some variables",
                                "missing. This is often the case when your CSV",
                                "dataset does not have the target",
                                "variable included (e.g., when your test dataset",
                                "is meant to be used as a scoring dataset, in which case",
                                "we can't perform an evaluation).",
                                "For producing a confusion matrix we need",
                                "to include the target variable.",
                                "Please load a CSV file which has",
                                "the risk and target variables included.",
                                "\n\n The actual error message was:\n\n%s"),
                           mtype, paste(result, "\n")))
      else
        errorReport(respcmd, result)
      next()
    }

    appendLog(Rtxt("Generate the confusion matrix showing counts."), confuse.cmd)

    confuse.output <- collectOutput(confuse.cmd, TRUE)

    appendLog(Rtxt("Generate the confusion matrix showing proportions."), percentage.cmd)
    percentage.output <- collectOutput(percentage.cmd)

    if (categoricTarget())
    {
      appendLog(Rtxt("Calculate the overall error percentage."), error.cmd)
      error.output <- collectOutput(paste(percentage.cmd, ";", error.cmd))
      appendLog(Rtxt("Calculate the averaged class error percentage."), avgerr.cmd)
      avgerr.output <- collectOutput(paste(percentage.cmd, ";", avgerr.cmd))
    }

    appendTextview(TV,
                   sprintf(Rtxt("Error matrix for the %s model",
                                "on %s (counts):"),
                           commonName(mtype), testname),
                   "\n\n",
                   confuse.output,
                   "\n\n",
                   sprintf(Rtxt("Error matrix for the %s model",
                                "on %s (proportions):"),
                           commonName(mtype), testname),
                   "\n\n",
                   percentage.output,
                   if (categoricTarget())
                   {
                     paste("\n\n", sprintf(Rtxt("Overall error: %s%%"),
                                           format(error.output)),
                           ", ",   sprintf(Rtxt("Averaged class error: %s%%"),
                                           format(avgerr.output)), sep="")
                   })

  }

  return(sprintf(Rtxt("Generated confusion matrix."), mtype, testname))
}

#----------------------------------------------------------------------
#
# EVALUATE RISK CHART
#

executeEvaluateRisk <- function(probcmd, testset, testname)
{
  # Initial setup.

  TV <- "risk_textview"
  resetTextview(TV)

  advanced.graphics <- theWidget("use_ggplot2")$getActive()

  # Ensure a risk variable has been specified. 081002 The plotRisk
  # function still works if no risk variable has been specified, so
  # let's go with it.

  risk <- crs$risk

  # Put 1 or 2 charts onto their own plots. Otherwise, put the
  # multiple charts onto one plot, keeping them all the same size
  # (thus if numplots is odd, leave a cell of the plot empty.

  numplots <- length(getEvaluateModels())
  if (numplots == 1)
    newPlot(1)
  else if (numplots == 2)
    newPlot(1)
  else if (numplots %% 2 == 0)
    newPlot(numplots)
  else
    newPlot(numplots + 1)

  if (numplots <= 2 )
    cex <- 1.0
  else if (numplots <= 4)
    cex <- 0.5
  else
    cex <- 0.5

  opar <- par(cex=cex)

  model.list <- getEvaluateModels()

  for (mtype in model.list)
  {

    setStatusBar(sprintf(Rtxt("Applying %s model to the dataset",
                              "to generate a risk chart ..."),
                         commonName(mtype)))

    # We need the base testset name here to get the risk variable, which
    # is not usually in the list of included columns.

    # testbase <- gsub(", ?c\\(.*\\]", ",]", testset)

    # Instead, obtain the column list, and if it exists, add the risk
    # variable to it, to avoid listing all columns, since this can
    # affect the na.omit function which will omit more rows if these
    # extra columns have NAs.

    # If using new graphics option then check ggplot2 is installed.
      
    if (advanced.graphics)
    {
      lib.cmd <- "library(ggplot2)"
      if (! packageIsAvailable("ggplot2", "plot a riskchart")) return()
    }
      
    if (length(crs$risk))
    {
      # Extract the columns selected from the test dataset as we will
      # augment this with the risk variable name as we may need to get
      # the same rows removed through NAs, before then extracting the
      # relevant risk variable.

      testcols <- sub("]$", "", sub("[^,]*, ", "", testset[[mtype]]))
      if (testcols != "")
      {
        newcols <- gsub(")", sprintf(", %s)", "crs$risk"), testcols)
        testsetr <- gsub(testcols, newcols, testset[[mtype]], fixed=TRUE)
      }

      evaluate.cmd <- paste("crs$eval <- evaluateRisk(crs$pr,",
                            sprintf("\n    %s$%s,", testset[[mtype]], crs$target),
                            sprintf("\n    %s$%s)", testsetr, risk))


      if (advanced.graphics)
        plot.cmd <- paste("print(riskchart(crs$pr,",
                          '\n                ',
                          sprintf("%s$%s, ", testset[[mtype]], crs$target),
                          '\n                ',
                          sprintf("%s$%s, ", testsetr, risk),
                          '\n                ',
                          'title="',
                          paste("Risk Chart", commonName(mtype), testname, crs$target,
                                '", '),
                          '\n                ',
                          'risk.name="', risk, '", recall.name="', crs$target, '",',
                          '\n                ',
                          'show.lift=', ifelse(numericTarget(), "FALSE", "TRUE"),
                          ', show.precision=', ifelse(numericTarget(), "FALSE", "TRUE"),
                          ', legend.horiz=FALSE',
                          '))\n',
                          sep="")
      else
        plot.cmd <- paste("plotRisk(crs$eval$Caseload, ",
                          "crs$eval$Precision, crs$eval$Recall, crs$eval$Risk,",
                          '\n    risk.name="', risk,
                          '", recall.name="', crs$target, '"',
                          ', show.lift=', ifelse(numericTarget(), "FALSE", "TRUE"),
                          ', show.precision=', ifelse(numericTarget(), "FALSE", "TRUE"),
                          ")\n",
                          genPlotTitleCmd(Rtxt("Risk Chart"), commonName(mtype),
                                        testname, risk),
                        sep="")
    }
    else
    {
      evaluate.cmd <- paste("crs$eval <- evaluateRisk(crs$pr,",
                            sprintf("%s$%s)", testset[[mtype]], crs$target))

      if (advanced.graphics)
        plot.cmd <- paste("print(riskchart(crs$pr, ",
                          '\n                ',
                          sprintf("%s$%s, ", testset[[mtype]], crs$target),
                          '\n                ',
                          'title="',
                          paste("Performance Chart", commonName(mtype), testname, '"'),
                          ', show.lift=', ifelse(numericTarget(), "FALSE", "TRUE"),
                          ', show.precision=', ifelse(numericTarget(), "FALSE", "TRUE"),
                          ', legend.horiz=FALSE',
                          '))\n',
                          sep="")
      else
        plot.cmd <- paste("plotRisk(crs$eval$Caseload, ",
                          "crs$eval$Precision, crs$eval$Recall",
                          ', show.lift=', ifelse(numericTarget(), "FALSE", "TRUE"),
                          ', show.precision=', ifelse(numericTarget(), "FALSE", "TRUE"),
                          ")\n",
                          genPlotTitleCmd("Performance Chart", commonName(mtype),
                                          testname),
                          sep="")
    }

    if (advanced.graphics)
    {
      appendLog(Rtxt("Risk Chart: requires the ggplot2 package."), lib.cmd)
      # 140906 Move to using namespaces within the code, though still
      # expose the interactive commands.
      eval(parse(text=lib.cmd))
    }

    appendLog(Rtxt("Generate a risk chart."),
              ifelse(advanced.graphics,
                     Rtxt("# Rattle provides evaluateRisk() and riskchart()."),
                     Rtxt("# Rattle provides evaluateRisk() and plotRisk().")),
              "\n\n",
              probcmd[[mtype]], "\n",
              evaluate.cmd, "\n",
              plot.cmd, sep="")

    result <- try(eval(parse(text=probcmd[[mtype]])), TRUE)

    # Check for errors - in particular, new levels in the test dataset.

    if (inherits(result, "try-error"))
    {
      if (any(grep("has new level", result)) || any(grep("New levels",result)))
        infoDialog("It seems that the dataset on which the probabilities",
                   "from the", mtype, "model are required has a categoric",
                   "variable with levels not found in the training",
                   "dataset. The probabilities can not be determined in",
                   "this situation. You may need to either ensure",
                   "the training dataset has representatives of all levels",
                   "or else remove them from the testing dataset.",
                   "Alternatively, do not include that variable in the",
                   "modelling. \n\n The actual error message was:",
                   "\n\n", paste(result, "\n"))
      else if (any(grep("undefined columns", result)))
        infoDialog("It seems that the dataset on which the predictions",
                   "from the", mtype, "model are required has some variables",
                   "missing. This is often the case when your CSV",
                   "dataset does not have the risk or target",
                   "variables included (e.g., when your test dataset",
                   "is meant to be used as a scoring dataset, in which case",
                   "we can't perform an evaluation).",
                   "For producing a risk chart we need",
                   "to include the risk and target variables.",
                   "Please load a CSV file which has",
                   "the risk and target variables included.",
                   "\n\n The actual error message was:",
                   "\n\n", paste(result, "\n"))
      else
        errorReport(probcmd, result)
      next()
    }

    # Check for all results the same.

    if (length(levels(as.factor(crs$pr))) == 1)
    {
      errorDialog("The model predicts the same result for all records,",
                  "so there is nothing to plot!")
      return()
    }

    # Now generate a summary of the performance at various probability
    # cutoffs, with the result being stored in crs$eval. We really
    # only do this for the older style plot as now riskchart() does
    # the evaluateRisk itself. But we also want to display the textual
    # results, so keep it here.

    eval(parse(text=evaluate.cmd))

    # We want to display the numeric results of the evaluation. But if
    # there are too many rows, as produced by KSVM for example, it will
    # be too much, so limit it to 100 row, which need to be selected
    # every Nth.

    ne <- nrow(crs$eval)
    maxev <- 100
    if (ne > maxev)
    {
      id <- round(seq(1, ne, length=maxev))
      msg <- sprintf("The sequence has been truncated to just %d from %d.\n\n",
                     maxev, ne)
    }
    else
    {
      id <- seq_len(ne)
      msg <- ""
    }
    id <- sprintf("c(%s)", paste(id, collapse=","))
    msg <- paste("Summary ", commonName(mtype), " model ",
                 sprintf("(built using %s)", mtype), " on ",
                 testname,
                 " by ",
                 ifelse(numericTarget(), "predicted value", "probability"),
                 " cutoffs.\n\n", msg, sep="")
    appendTextview(TV, msg, collectOutput(sprintf("crs$eval[%s,]", id), TRUE))

    # Display the AUC measures.

    #auc <- calculateRiskAUC(crs$eval)
    #print(auc)
    if (length(crs$risk))
    {
      aucRisk <- with(crs$eval, calculateAUC(Caseload, Risk))
      aucRisk <- (2*aucRisk)/(2-crs$eval$Precision[1])
    }
      
    aucRecall <- with(crs$eval, calculateAUC(Caseload, Recall))
    aucRecall <- (2*aucRecall)/(2-crs$eval$Precision[1])
    appendTextview(TV, paste("The area under the ",
                             if (length(crs$risk))
                             "Recall curve " else "Risk and Recall curves ",
                             "for ", commonName(mtype), " model\n\n",
                             "Area under the Recall (green) curve: ",
                             sprintf("%d%% (%0.3f)\n",
                                     round(100*aucRecall), aucRecall),
                             if (length(crs$risk))
                             "Area under the Risk   (red)   curve: ",
                             if (length(crs$risk))
                             sprintf("%d%% (%0.3f)",
                                     round(100*aucRisk), aucRisk),
                             sep=""))

    # Display the Risk Chart itself now.

    # For 2 plots, so as not to overwrite the first plot, if we are
    # about to plot the second plot, initiate a new plot.

    if (advanced.graphics && dev.cur() != 2) # 121219 Not sure why 2?
      newPlot()
    else if (numplots == 2 && mtype == model.list[length(model.list)])
      newPlot(1)

    eval(parse(text=plot.cmd))

  }

  # Restore par

  par(opar)

  return(sprintf(Rtxt("Generated %d risk chart%s."),
                 numplots, ifelse(numplots>1, "s", "")))

}

grid.plot <- function (colour="gray", tops=100)
{
  opar = par(lwd=1)
  abline(v=seq(0,tops,tops/10), col=colour, lty="dotted")
  abline(h=seq(0,tops,tops/10), col=colour, lty="dotted")
  par(opar)
}

plotOptimalLine <- function(x, y1, y2, pr=NULL, colour="plum", label=NULL)
{
  lines(c(x, x), c(-13, max(y1, y2)), lty=6, col=colour)
  lines(c(-13, x), c(y1, y1), lty=6, col=colour)
  lines(c(-13, x), c(y2, y2), lty=6, col=colour)
  if (not.null(label))
  {
    text(x, 0, label, pos=4)
    text(x, 0, sprintf("%2.0f%%", x), pos=2)
    text(0, y2, sprintf("%2.0f%%", y2), pos=3, offset=0.2)
    text(0, y1, sprintf("%2.0f%%", y1), pos=3, offset=0.2)
    if (not.null(pr))
      text(x, pr+4, sprintf("%2.0f%%", pr), pos=2)
  }
}

evaluateRisk <- function(predicted, actual, risks=NULL)
{
  # 081002 We allow risk to be not specified.

  if (is.factor(actual))
    actual <- as.integer(actual)-1

  # With na.rm=TRUE we cater for the case when the actual data has
  # missing values for the target.

  # 090802 Try allowing any range of numbers

#  if (min(actual, na.rm=TRUE) != 0 || max(actual, na.rm=TRUE) !=1 )
#    stop("actual must be binary (0,1) but found (",
#         min(actual, na.rm=TRUE), ",", max(actual, na.rm=TRUE), ").")

  # For KSVMs, and perhaps other modellers, the predictied values are
  # probabilites, which may be a very high level of precision (e.g.,
  # 0.999999999999996 or 2.58015830922886e-13), and thus, when
  # converted to a factor, we have almost a one-to-one match from an
  # entity to a probability. When converted to a data frame the
  # resulting row names (these probablities of being a 1) have
  # caseloads of 1, 2, or 3, thus there are very many, and sometimes,
  # the probablities are the same! We then get duplicate row names and
  # the assigning of new names to the columns below gives an error
  # about duplicate row names! We should aggregate up to three
  # significant figures in the probabilities to make everything much
  # easier. BUT this then lumps all of the 0.9999999.... together, and
  # leaves a very large jump at the right end of the plot! We really
  # might want to instead aggregate on caseload! But rounding it to 13
  # digits seems okay! We get a good plot.

  predicted <- as.factor(round(predicted, 13))

  if (is.null(risks))
    ds.actual <- data.frame(Actual=actual,
                            Predict=predicted)
  else
    ds.actual <- data.frame(Actual=actual,
                            Risk=as.numeric(risks), # Avoid integer overflow
                            Predict=predicted)
  #Predict=as.factor(ds.predict[,2]))

  # With na.rm=TRUE in the first sum here we cater for the case when
  # the actual data has missing values for the target.

  ds.evaluation <- as.data.frame(t(rbind(tapply(ds.actual$Actual,
                                                ds.actual$Predict,
                                                sum, na.rm=TRUE),
                                         if (is.null(risks))
                                         NULL
                                         else
                                         tapply(ds.actual$Risk,
                                                ds.actual$Predict,
                                                sum, na.rm=TRUE),
                                         tapply(ds.actual$Actual,
                                              ds.actual$Predict, length))))

  if (is.null(risks))
    colnames(ds.evaluation) <- c("Recall", "Caseload")
  else
    colnames(ds.evaluation) <- c("Recall", "Risk", "Caseload")

  last <- nrow(ds.evaluation)
  ds.evaluation$Precision[last] <- ds.evaluation$Recall[last]/
    ds.evaluation$Caseload[last]

  for (i in (nrow(ds.evaluation)-1):1)
    {
      ds.evaluation$Recall[i] <- ds.evaluation$Recall[i+1] +
        ds.evaluation$Recall[i]
      if (not.null(risks))
        ds.evaluation$Risk[i] <- ds.evaluation$Risk[i+1] +
          ds.evaluation$Risk[i]
      ds.evaluation$Caseload[i] <- ds.evaluation$Caseload[i+1] +
        ds.evaluation$Caseload[i]
      ds.evaluation$Precision[i] <- ds.evaluation$Recall[i] /
        ds.evaluation$Caseload[i]
    }
  ds.evaluation$Recall <- ds.evaluation$Recall/ds.evaluation$Recall[1]
  if (not.null(risks)) ds.evaluation$Risk <- ds.evaluation$Risk/ds.evaluation$Risk[1]
  ds.evaluation$Caseload <- ds.evaluation$Caseload/ds.evaluation$Caseload[1]
  # This is Michael's measure of performance.
  if (not.null(risks))
    ds.evaluation$Measure <- abs(ds.evaluation$Recall - ds.evaluation$Caseload) +
      abs(ds.evaluation$Risk - ds.evaluation$Caseload)

  # 120502 Be sure we include the 0, 0 point. 140503 If Caseload=0 is
  # not present in the data, then add one. Previously used a check on
  # whether row.names included a 1 but that does not acutally catch
  # the case where the risk scores include 1. Recall, and Risk are 0
  # and Strike Rate is 1. That is, no caseload means no recall, no
  # risk, no measure and let's say 100% strike rate.

  if (! 0 %in% ds.evaluation$Caseload)
  {
    if (is.null(risks))
      ds.evaluation <- rbind(ds.evaluation, c(0.00, 0.00, 1.00))
    else
      ds.evaluation <- rbind(ds.evaluation, c(0.00, 0.00, 0.00, 1.00, 0.00))
  
    row.names(ds.evaluation)[nrow(ds.evaluation)] <- "1.0"
  }
  
  return(ds.evaluation)
}

calculateAUC <- function(x, y)
{
  len <- length(x)
  ria <- x[len] * y[len] / 2

  for (i in (len-1):1)
  {
    ria <- ria +
      (x[i] - x[i+1]) * y[i+1] + (x[i] - x[i+1]) * (y[i] - y[i+1]) / 2
  }
  return(ria)
}

openMyDevice <- function(dev, filename)
{
  if (dev == "" && filename != "")
  {
    fn <- unlist(strsplit(filename, "\\."))
    dev=fn[length(fn)]
  }

  if (dev == "wmf")
    eval(parse(text=sprintf("win.metafile(%s)", filename)))
  else if (dev == "png")
    png(filename)
  else if (dev == "pdf")
    pdf(filename)

  return(dev)

}

plotRisk <- function (cl, pr, re, ri=NULL,
                      title=NULL,
                      show.legend=TRUE,
                      xleg=60, yleg=55,
                      # xleg=20, yleg=16, # gets lost with  multiplots.
                      optimal=NULL, optimal.label="",
                      chosen=NULL, chosen.label="",
                      include.baseline=TRUE,
                      dev="", filename="",
                      show.knots=NULL,
                      show.lift=TRUE,
                      show.precision=TRUE,
                      risk.name="Risk",
                      recall.name="Recall",
                      precision.name="Precision")
{
  openMyDevice(dev, filename)

  ## If proportions, convert to percentages

  if (all(cl <= 1)) cl <- cl * 100
  if (all(re <= 1)) re <- re * 100
  if (not.null(ri) && all(ri <= 1.5)) ri <- ri * 100 # Can sometimes be just >1
  if (all(pr <= 1)) pr <- pr * 100
  #
  # If list is in min to max order then reverse
  #
  if (cl[1] < cl[length(cl)])
  {
    cl <- rev(cl)
    pr <- rev(pr)
    re <- rev(re)
    ri <- rev(ri)
  }
  #
  # Add a zero point for the display
  #
  cl <- c(cl, 0)
  re <- c(re, 0)
  if (not.null(ri)) ri <- c(ri, 0)
  pr <- c(pr, NA)
  #
  # Also add the 100 point just in case?
  #
  if (cl[1] != 100)
  {
    cl <- c(100, cl)
    re <- c(100, re)
    if (not.null(ri)) ri <- c(100, ri)
    pr <- c(min( pr[!is.na(pr)]), pr)
  }
  #
  # Now plot
  #
  if (show.lift)
    opar <- par(lwd=2, mar=c(5.1, 4.1, 4.1, 4.1))
  else
    opar <- par(lwd=2)
  plot(c(0,100), c(0,100), type='l', col=1,
       xlab=Rtxt("Caseload (%)"), ylab=Rtxt("Performance (%)"),
       ylim=c(0,100), xlim=c(0,100))
  grid.plot()
  if (not.null(title))
    title(main=title, sub=paste("Rattle", Sys.time(), Sys.info()["user"]))
  points(re ~ cl, type='l', col=3, lty=5)
  if (show.precision) points(pr ~ cl, type='l', col=4, lty=4)
  if (not.null(ri)) points(ri ~ cl, type='l', col=2, lty=1)
  if (include.baseline && show.precision) text(100, pr[1]+4, sprintf("%0.0f%%", pr[1]))
  # Optimal
  if (not.null(optimal))
  {
    optimal.index <- which(abs(cl-optimal) == min(abs(cl-optimal)))
    if (length(optimal.index) > 1) optimal.index <- optimal.index[1]
    plotOptimalLine(optimal, ri[optimal.index], re[optimal.index],
                      pr[optimal.index], label=optimal.label)
  }
  # Chosen
  if (not.null(chosen))
  {
    chosen.index <- which(abs(cl-chosen) == min(abs(cl-chosen)))
    if (length(chosen.index) > 1) chosen.index <- chosen.index[1]
    plotOptimalLine(chosen, ri[chosen.index], re[chosen.index],
                      label=chosen.label, colour="grey")
  }

  legend <- c()
  lty <- c()
  col <- c()
  if (not.null(ri))
  {
    auc <- calculateAUC(cl/100, ri/100)
    legend <- c(legend, sprintf("%s (%d%%)", risk.name, round(100*auc)))
    lty <- c(lty, 1)
    col <- c(col, 2)
  }
  auc <- calculateAUC(cl/100, re/100)
  legend <- c(legend, sprintf("%s (%d%%)", recall.name, round(100*auc)))
  if (show.precision) legend <- c(legend, precision.name)
  lty <- c(lty,5,4)
  col <- c(col,3,4)
  if (not.null(optimal))
  {
    legend <- c(legend, Rtxt("Optimal"))
    lty <- c(lty,6)
    col <- c(col,"plum")
  }
  if (not.null(chosen))
  {
    legend <- c(legend, Rtxt("Chosen"))
    lty <- c(lty,6)
    col <- c(col,"grey")
  }
  if (show.legend)
    legend(xleg, yleg, legend, lty=lty, lwd=2, col=col, bty="n")

  if (show.lift)
  {
    lifts <- seq(pr[1], 100, pr[1])
    axis(4, at=lifts, labels=seq(1, length(lifts)))
    mtext(Rtxt("Lift"), side=4, line=3)
  }

  #
  #
  # Add in knot labels
  #
  if (not.null(show.knots))
  {
    len <- length(cl)
    text(cl[c(-1,-len)]-2, ri[c(-1,-len)]+3, rev(show.knots)[-1])
  }
  if (dev != "") dev.off()
  par(opar)
}

handleMissingValues <- function(testset, mtype)
{
  # 091205 A shared function to generate predictions in the variable
  # pred that do not include any missing values (from the target
  # variable).

  return(paste("",
               Rtxt("# Remove observations with missing target."),
               '',
               sprintf('no.miss   <- na.omit(%s$%s)', testset[[mtype]], crs$target),
               'miss.list <- attr(no.miss, "na.action")',
               # 110320 Avoid prediction complaining about invalid labels.
               'attributes(no.miss) <- NULL',
               '',
               'if (length(miss.list))',
               '{',
               '  pred <- prediction(crs$pr[-miss.list], no.miss)',
               '} else',
               '{',
               '  pred <- prediction(crs$pr, no.miss)',
               '}',
               sep="\n"))
}

#doRiskChart <- function(pr, data, test, target, risk, main, thresholds=NULL)
#{
    # TODO 20130727 REMOVE this function sine now riskchart take pr,
    # ac, ri itself so this wrapper is no longer required.
#    riskchart(pr, data[test, target], data[test, risk],
#              risk.name=risk, recall.name=target, show.lift=TRUE,
#            show.precision=TRUE, title=main, thresholds=thresholds)
  
#  plotRisk(eval$Caseload, eval$Precision, eval$Recall, eval$Risk,
#           risk.name=risk,
#           recall.name=target,
#           show.lift=TRUE,
#           show.precision=TRUE)
  
#  title(main=main,
#        sub=paste("Rattle", format(Sys.time(), "%Y-%b-%d %H:%M:%S"),
#          Sys.info()["user"]))
#}

#----------------------------------------------------------------------
# EVALUATE COST CURVE 080524

executeEvaluateCostCurve <- function(probcmd, testset, testname)
{
  # 080524 Display Cost Curves (Drummond and Holte)

  lib.cmd <- "library(ROCR)"
  if (! packageIsAvailable("ROCR", "plot a cost curve")) return()

  # Put 1 or 2 charts onto their own plots. Otherwise, put the
  # multiple charts onto one plot, keeping them all the same size
  # (thus if numplots is odd, leave a cell of the plot empty.

  numplots <- length(getEvaluateModels())
  if (numplots == 1)
    newPlot(1)
  else if (numplots == 2)
    newPlot(1)
  else if (numplots %% 2 == 0)
    newPlot(numplots)
  else
    newPlot(numplots + 1)

  if (numplots <= 2 )
    cex <- 1.0
  else if (numplots <= 4)
    cex <- 0.5
  else
    cex <- 0.5

  opar <- par(cex=cex)

  nummodels <- length(probcmd)
  #140906 mcolors is not used? Why define it??? REMOVE.
  #  if (packageIsAvailable("colorspace"))
  #     mcolors <- colorspace::rainbow_hcl(nummodels) # 090524, start = 270, end = 150)
  #  else
  #    mcolors <- rainbow(nummodels, 1, .8)
  mcount <- 0

  model.list <- getEvaluateModels()

  for (mtype in model.list)
  {
    setStatusBar("Applying", commonName(mtype),
                 "model to the dataset to generate a cost curve ...")

    mcount <- mcount + 1
    plot.cmd <- paste(paste("plot(0, 0, xlim=c(0, 1), ylim=c(0, 1),",
                            sprintf('xlab="%s",', Rtxt("Probability cost function")),
                            sprintf('ylab="%s")', Rtxt("Normalized expected cost"))),
                      'lines(c(0,1),c(0,1))',
                      'lines(c(0,1),c(1,0))',
                      handleMissingValues(testset, mtype),
                      'perf1 <- performance(pred, "fpr", "fnr")',
                      'for (i in seq_along(perf1@x.values))\n{',
                      '\tfor (j in seq_along(perf1@x.values[[i]]))\n\t{',
                      '\t\tlines(c(0,1),c(perf1@y.values[[i]][j],',
                      '\t\t\t\tperf1@x.values[[i]][j]),',
                      '\t\t\t\tcol=terrain.colors(10)[i],lty=3)',
                      '\t}\n}',
                      'perf<-performance(pred, "ecost")',
                      '\n# Bug in ROCR 1.0-3 does not obey the add command.',
                      '# Calling the function directly does work.\n',
                      ".plot.performance(perf, lwd=1.5, xlim=c(0,1), ylim=c(0,1), add=T)",
                      "op <- par(xpd=TRUE)",
                      'text(0, 1.07, "FPR")',
                      'text(1, 1.07, "FNR")',
                      "par(op)",
                      'text(0.12, 1, "Predict +ve")',
                      'text(0.88, 1, "Predict -ve")',
                      # TODO 080810 Add text AUC=... to plot
                      genPlotTitleCmd(Rtxt("Cost Curve"), commonName(mtype),
                                      testname),
                      sep="\n")

    appendLog(Rtxt("Cost Curve: requires the ROCR package."), lib.cmd)
    # 140906 Move to using namespaces within the code, though still
    # expose the interactive commands.
    eval(parse(text=lib.cmd))

    appendLog(sprintf(Rtxt("Generate a Cost Curve for the %s model on %s."),
                     commonName(mtype), testname),
             probcmd[[mtype]], "\n", plot.cmd)

    result <- try(eval(parse(text=probcmd[[mtype]])), silent=TRUE)

    # Check for errors - in particular, new levels in the test dataset.

    if (inherits(result, "try-error"))
    {
      if (any(grep("has new level", result)) || any(grep("New levels",result)))
        infoDialog(sprintf(Rtxt("It seems that the dataset on which the probabilities",
                                "from the %s model are required has a categoric",
                                "variable with levels not found in the training",
                                "dataset. The probabilities can not be determined in",
                                "this situation. You may need to either ensure",
                                "the training dataset has representatives of all levels",
                                "or else remove them from the testing dataset.",
                                "Alternatively, do not include that variable in the",
                                "modelling. \n\nThe actual error message was:"),
                           mtype),
                   "\n\n", paste(result, "\n"))
      else
        errorReport(probcmd, result)
      next()
    }

    # Display the Cost Curve itself now.

    # For 2 plots, so as not to overwrite the first plot, if we are
    # about to plot the second plot, initiate a new plot.

    if (numplots == 2 && mtype == model.list[length(model.list)]) newPlot(1)

    eval(parse(text=plot.cmd))

  }
  return(sprintf(Rtxt("Generated Cost Curves on %s."), testname))
}


#----------------------------------------------------------------------
# LIFT CHART

executeEvaluateLift <- function(probcmd, testset, testname)
{
  lib.cmd <- "library(ROCR)"
  if (! packageIsAvailable("ROCR", Rtxt("plot a lift chart"))) return()

  newPlot()
  addplot <- "FALSE"
  xlab <- Rtxt("Caseload (%)")

  nummodels <- length(probcmd)
  mcolors <- rainbow(nummodels, 1, .8)
  mcount <- 0

  for (mtype in getEvaluateModels())
  {
    setStatusBar(sprintf(Rtxt("Applying %s model to the dataset to generate",
                              "a lift chart ..."), mtype))

    mcount <- mcount + 1
    plot.cmd <- paste(handleMissingValues(testset, mtype),
                      "",
                      Rtxt("# Convert rate of positive predictions to percentage."),
                      "",
                      'per <- performance(pred, "lift", "rpp")',
                      "per@x.values[[1]] <- per@x.values[[1]]*100\n",
                      Rtxt("# Plot the lift chart."),
                      paste("ROCR::plot(per,",
                            sprintf('col="%s", lty=%d,', mcolors[mcount], mcount),
                            sprintf('xlab="%s",', xlab),
                            sprintf("add=%s)", addplot)),
                      sep="\n")
    addplot <- "TRUE"

    appendLog(Rtxt("Lift Chart: requires the ROCR package."), lib.cmd)
    # 140906 Move to using namespaces within the code, though still
    # expose the interactive commands.
    eval(parse(text=lib.cmd))

    # print(mtype); print(testname)

    appendLog(sprintf(Rtxt("Obtain %s for the %s model on %s."),
                     Rtxt("predictions"), mtype, testname),
             probcmd[[mtype]], "\n", plot.cmd)

    result <- try(eval(parse(text=probcmd[[mtype]])), silent=TRUE)

    # Check for errors - in particular, new levels in the test dataset.

    if (inherits(result, "try-error"))
    {
      if (any(grep("has new level", result)) || any(grep("New levels",result)))
        infoDialog(sprintf(Rtxt("It seems that the dataset on which the probabilities",
                                "from the %s model are required has a categoric",
                                "variable with levels not found in the training",
                                "dataset. The probabilities can not be determined in",
                                "this situation. You may need to either ensure",
                                "the training dataset has representatives of all levels",
                                "or else remove them from the testing dataset.",
                                "Alternatively, do not include that variable in the",
                                "modelling. \n\n The actual error message was:"),
                           mtype),
                  "\n\n",  paste(result, "\n"))
      else
        errorReport(probcmd, result)
      next()
    }

    eval(parse(text=plot.cmd))

  }

  # If just one model, and we are plotting the test dataset, then
  # also plot the training dataset.

  if (nummodels==1 && length(grep("\\[test\\]", testname))>0)
  {
    mcount <- mcount + 1
    plot.cmd <- paste(Rtxt("\n# Also convert rate of positive predictions to percentage\n"),
                      "\nper <- performance(prediction(crs$pr, ",
                      sprintf("%s$%s),", sub("-crs\\$sample", "crs$sample",
                                  testset[[mtype]]), crs$target),
                      '"lift", "rpp")\n',
                      "per@x.values[[1]] <- per@x.values[[1]]*100\n\n",
                      Rtxt("# Now plot the lift.\n"),
                      Rtxt("\n# Bug in ROCR 1.0-3 plot does not obey the add command."),
                      Rtxt("# Calling the function directly (.plot.performance) does work.\n"),
                      "\n.plot.performance(per, ",
                      'col="#00CCCCFF", lty=2, ',
                      sprintf("add=%s)", addplot),
                      sep="")
    appendLog(sprintf(Rtxt("Generate a Lift Chart for the %s model on %s."),
                     mtype, sub('\\[test\\]', '[train]', testname)),
             sub("-crs\\$sample", "crs$sample",
                                   probcmd[[mtype]]), "\n", plot.cmd)

    result <- try(eval(parse(text=sub("-crs\\$sample",
                               "crs$sample", probcmd[[mtype]]))), silent=TRUE)
    eval(parse(text=plot.cmd))
    models <- c(Rtxt("Test"), Rtxt("Train"))
    nummodels <- 2
    legtitle <- sapply(getEvaluateModels(), commonName)
    title <- sub('\\[test\\]', '', testname)
  }
  else
  {
    models <- getEvaluateModels()
    legtitle <- Rtxt("Models")
    title <- testname
  }

  legendcmd <- paste('legend("topright",',
                     sprintf("c(%s),",
                             paste('"', models, '"',
                                   sep="", collapse=",")),
                     sprintf('col=rainbow(%d, 1, .8), lty=1:%d,',
                             nummodels, nummodels),
                     sprintf('title="%s", inset=c(0.05, 0.05))', legtitle))
  appendLog(Rtxt("Add a legend to the plot."), legendcmd)
  eval(parse(text=legendcmd))

  decorcmd <- paste(genPlotTitleCmd(Rtxt("Lift Chart"), "", title),
                    '\ngrid()', sep="")
  appendLog(Rtxt("Add decorations to the plot."), decorcmd)
  eval(parse(text=decorcmd))

  return(Rtxt("Generated Lift Charts."))
}

#----------------------------------------------------------------------
# ROC PLOT

executeEvaluateROC <- function(probcmd, testset, testname)
{
  TV <- "roc_textview"
  resetTextview(TV)
  advanced.graphics <- theWidget("use_ggplot2")$getActive()

  lib.cmd <- "library(ROCR)"
  if (! packageIsAvailable("ROCR", Rtxt("plot an ROC curve"))) return()

  if (advanced.graphics)
  {
    req.ggplot2.cmd <- "library(ggplot2, quietly=TRUE)"
    if (! packageIsAvailable("ggplot2", "plot an ROC curve")) return()
  }
  else
  {
    newPlot()
    addplot <- "FALSE"
  }

  nummodels <- length(probcmd)
  mcolors <- rainbow(nummodels, 1, .8)
  mcount <- 0


  for (mtype in getEvaluateModels())
  {
    # A new plot for each chart.
    
    if (advanced.graphics) newPlot()

    # Inform the user what is going on.
    
    setStatusBar(sprintf(Rtxt("Applying %s model to the dataset to generate",
                              "an ROC plot ..."),
                         mtype))

    mcount <- mcount + 1

    plot.cmd <- handleMissingValues(testset, mtype)
    if (advanced.graphics)
      plot.cmd <-
        paste(plot.cmd,
              '',
              'pe <- performance(pred, "tpr", "fpr")',
              'au <- performance(pred, "auc")@y.values[[1]]',
              'pd <- data.frame(fpr=unlist(pe@x.values), tpr=unlist(pe@y.values))',
              'p <- ggplot(pd, aes(x=fpr, y=tpr))',
              'p <- p + geom_line(colour="red")',
              'p <- p + xlab("False Positive Rate") + ylab("True Positive Rate")',
              sprintf('p <- p + ggtitle("ROC Curve %s %s %s")',
                      commonName(mtype), testname, crs$target),
              'p <- p + theme(plot.title=element_text(size=10))',
              paste('p <- p + geom_line(data=data.frame(), aes(x=c(0,1), y=c(0,1)),',
                    'colour="grey")'),
              'p <- p + annotate("text", x=0.50, y=0.00, hjust=0, vjust=0, size=5,',
              '                   label=paste("AUC =", round(au, 2)))',
              'print(p)',
              sep="\n")
    else
      plot.cmd <-
        paste(plot.cmd,
              '\n',
              'ROCR::plot(performance(pred, "tpr", "fpr"), ',
              sprintf('col="%s", lty=%d, ', mcolors[mcount], mcount),
              sprintf("add=%s)", addplot),
              sep="")
    
    addplot <- "TRUE"

    appendLog(Rtxt("ROC Curve: requires the ROCR package."), lib.cmd)
    # 140906 Move to using namespaces within the code, though still
    # expose the interactive commands.
    eval(parse(text=lib.cmd))

    if (advanced.graphics)
    {
      appendLog(Rtxt("ROC Curve: requires the ggplot2 package."), req.ggplot2.cmd)
      eval(parse(text=req.ggplot2.cmd))
    }
    
    appendLog(sprintf(Rtxt("Generate an ROC Curve for the %s model on %s."),
                     mtype, testname),
             probcmd[[mtype]], "\n", plot.cmd)

    result <- try(eval(parse(text=probcmd[[mtype]])), silent=TRUE)

    ## Check for errors - in particular, new levels in the test dataset.

    if (inherits(result, "try-error"))
    {
      if (any(grep("has new level", result)) || any(grep("New levels",result)))
        infoDialog(sprintf(Rtxt("It seems that the dataset on which the probabilities",
                                "from the %s model are required has a categoric",
                                "variable with levels not found in the training",
                                "dataset. The probabilities can not be determined in",
                                "this situation. You may need to either ensure",
                                "the training dataset has representatives of all levels",
                                "or else remove them from the testing dataset.",
                                "Alternatively, do not include that variable in the",
                                "modelling. \n\n The actual error message was:"),
                           mtype),
                   "\n\n",
                   paste(result, "\n"))
      else
        errorReport(probcmd, result)
      next()
    }

    eval(parse(text=plot.cmd))

    # Report the area under the curve.

    auc.cmd <- paste(handleMissingValues(testset, mtype),
                     'performance(pred, "auc")', sep="\n")
    appendLog(Rtxt("Calculate the area under the curve for the plot."), auc.cmd)
    auc <- eval(parse(text=auc.cmd))
    appendTextview(TV, sprintf(Rtxt("Area under the ROC curve for the",
                                    "%s model on %s is %0.4f"),
                               mtype, testname, attr(auc, "y.values")))
  }
  if (! advanced.graphics)
  {
    lines(c(0,1), c(0,1)) # Baseline

    # If just one model, and we are plotting the test dataset, then
    # also plot the training dataset.

    if (nummodels==1 && length(grep(sprintf("\\[%s\\]", Rtxt("test")), testname))>0)
    {
      mcount <- mcount + 1
      plot.cmd <- paste(Rtxt("\n# In ROCR (1.0-3) plot does not obey the add command.\n"),
                        Rtxt("# Calling the function directly works.\n\n"),
                        ".plot.performance(performance(prediction(crs$pr, ",
                        sprintf("%s$%s),",
                                sub("-crs\\$sample", "crs$sample",
                                    testset[[mtype]]), crs$target),
                        '"tpr", "fpr"), ',
                        'col="#00CCCCFF", lty=2, ',
                        sprintf("add=%s)\n", addplot),
                        sep="")
      appendLog(sprintf(Rtxt("Generate an ROC curve for the %s model on %s."),
                        mtype, sub(sprintf("\\[%s\\]", Rtxt("test")), '[train]', testname)),
                sub("-crs\\$sample", "crs$sample",
                    probcmd[[mtype]]), "\n", plot.cmd)
      
      result <- try(eval(parse(text=sub("-crs\\$sample",
                                   "crs$sample", probcmd[[mtype]]))), silent=TRUE)
      eval(parse(text=plot.cmd))
      models <- c(Rtxt("Test"), Rtxt("Train"))
      nummodels <- 2
      legtitle <- getEvaluateModels()
      title <- sub(sprintf("\\[%s\\]", Rtxt("test")), '', testname)
    }
    else
    {
      models <- getEvaluateModels()
      legtitle <- Rtxt("Models")
      title <- testname
    }
    
    legendcmd <- paste('legend("bottomright",',
                       sprintf("c(%s),",
                               paste('"', models, '"',
                                     sep="", collapse=",")),
                       sprintf('col=rainbow(%d, 1, .8), lty=1:%d,',
                               nummodels, nummodels),
                       sprintf('title="%s", inset=c(0.05, 0.05))', legtitle))
    appendLog(Rtxt("Add a legend to the plot."), legendcmd)
    eval(parse(text=legendcmd))

    decor.cmd <- paste(genPlotTitleCmd(Rtxt("ROC Curve"), "", title),
                       '\ngrid()', sep="")
    appendLog(Rtxt("Add decorations to the plot."), decor.cmd)
    eval(parse(text=decor.cmd))
  }

  return(sprintf(Rtxt("Generated ROC Curves on %s."), testname))
}

#----------------------------------------------------------------------
# EVALUATE PRECISION PLOT

executeEvaluatePrecision <- function(probcmd, testset, testname)
{
  lib.cmd <- "library(ROCR)"
  if (! packageIsAvailable("ROCR", Rtxt("plot a precision chart"))) return()

  newPlot()
  addplot <- "FALSE"

  nummodels <- length(probcmd)
  mcolors <- rainbow(nummodels, 1, .8)
  mcount <- 0

  for (mtype in getEvaluateModels())
  {
    setStatusBar(sprintf(Rtxt("Applying %s model to the dataset to generate",
                              "a Precision/Recall plot ..."), mtype))

    mcount <- mcount + 1

    plot.cmd <- paste(handleMissingValues(testset, mtype),
                      '\nROCR::plot(performance(pred, "prec", "rec"), ',
                      sprintf('col="%s", lty=%d, ', mcolors[mcount], mcount),
                      sprintf("add=%s)\n", addplot),
                      sep="")
    addplot <- "TRUE"

    appendLog(Rtxt("Precision/Recall Plot: requires the ROCR package"), lib.cmd)
    # 140906 Move to using namespaces within the code, though still
    # expose the interactive commands.
    eval(parse(text=lib.cmd))

    appendLog(sprintf(Rtxt("Generate a Precision/Recall Plot for the %s model on %s."),
                     mtype, testname),
             probcmd[[mtype]], "\n", plot.cmd)

    result <- try(eval(parse(text=probcmd[[mtype]])), silent=TRUE)

    ## Check for errors - in particular, new levels in the test dataset.

    if (inherits(result, "try-error"))
    {
      if (any(grep("has new level", result)) || any(grep("New levels",result)))
        infoDialog(sprintf(Rtxt("It seems that the dataset on which the probabilities",
                                "from the %s model are required has a categoric",
                                "variable with levels not found in the training",
                                "dataset. The probabilities can not be determined in",
                                "this situation. You may need to either ensure",
                                "the training dataset has representatives of all levels",
                                "or else remove them from the testing dataset.",
                                "Alternatively, do not include that variable in the",
                                "modelling. \n\n The actual error message was:"),
                           mtype), "\n\n", paste(result, "\n"))
      else
        errorReport(probcmd, result)
      next()
    }

    eval(parse(text=plot.cmd))
  }

  ## If just one model, and we are plotting the test dataset, then
  ## also plot the training dataset.

  if (nummodels==1 && length(grep("\\[test\\]", testname))>0)
  {
    mcount <- mcount + 1
    plot.cmd <- paste(Rtxt("\n# In ROCR (1.0-3) plot does not obey the add command.\n"),
                      Rtxt("# Calling the function directly (.plot.performance) does work.\n\n"),
                      ".plot.performance(performance(prediction(crs$pr, ",
                      sprintf("%s$%s),",
                              sub("-crs\\$sample", "crs$sample",
                                  testset[[mtype]]), crs$target),
                      '"prec", "rec"), ',
                      'col="#00CCCCFF", lty=2, ',
                      sprintf("add=%s)\n", addplot),
                      sep="")
    appendLog(sprintf(Rtxt("Generate a Precision/Recall Plot for the %s model on %s."),
                     mtype, sub('\\[test\\]', '[train]', testname)),
             sub("-crs\\$sample", "crs$sample",
                                   probcmd[[mtype]]), "\n", plot.cmd)

    result <- try(eval(parse(text=sub("-crs\\$sample",
                               "crs$sample", probcmd[[mtype]]))), silent=TRUE)
    eval(parse(text=plot.cmd))
    models <- c(Rtxt("Test"), Rtxt("Train"))
    nummodels <- 2
    legtitle <- getEvaluateModels()
    title <- sub('\\[test\\]', '', testname)
  }
  else
  {
    models <- getEvaluateModels()
    legtitle <- Rtxt("Models")
    title <- testname
  }

  legendcmd <- paste('legend("bottomleft",',
                     sprintf("c(%s),",
                             paste('"', models, '"',
                                   sep="", collapse=",")),
                     sprintf('col=rainbow(%d, 1, .8), lty=1:%d,',
                             nummodels, nummodels),
                     sprintf('title="%s", inset=c(0.05, 0.05))', legtitle))
  appendLog(Rtxt("Add a legend to the plot."), legendcmd)
  eval(parse(text=legendcmd))

  decor.cmd <- paste(genPlotTitleCmd(Rtxt("Precision/Recall Plot"), "", title),
                    '\ngrid()', sep="")
  appendLog(Rtxt("Add decorations to the plot."), decor.cmd)
  eval(parse(text=decor.cmd))

  return(sprintf(Rtxt("Generated Precision/Recall Plot on %s."), title))
}

#----------------------------------------------------------------------
# EVALUATE SENSITIVITY PLOT

executeEvaluateSensitivity <- function(probcmd, testset, testname)
{
  lib.cmd <- "library(ROCR)"
  if (! packageIsAvailable("ROCR", Rtxt("plot a sensitivity chart"))) return()

  newPlot()
  addplot <- "FALSE"

  nummodels <- length(probcmd)
  mcolors <- rainbow(nummodels, 1, .8)
  mcount <- 0

  for (mtype in getEvaluateModels())
  {
    setStatusBar(sprintf(Rtxt("Applying %s model to the dataset to generate",
                              "a Sensitivity plot ..."),
                         mtype))

    mcount <- mcount + 1
    plot.cmd <- paste(handleMissingValues(testset, mtype),
                      '\nROCR::plot(performance(pred, "sens", "spec"), ',
                      sprintf('col="%s", lty=%d, ', mcolors[mcount], mcount),
                      sprintf("add=%s)\n", addplot),
                      sep="")
     addplot <- "TRUE"

    appendLog(Rtxt("Sensitivity/Specificity Plot: requires the ROCR package"),
             lib.cmd)
    eval(parse(text=lib.cmd))

    appendLog(sprintf(Rtxt("Generate Sensitivity/Specificity Plot for %s model on %s."),
                     mtype, testname),
             probcmd[[mtype]], "\n", plot.cmd)

    result <- try(eval(parse(text=probcmd[[mtype]])), silent=TRUE)

    ## Check for errors - in particular, new levels in the test dataset.

    if (inherits(result, "try-error"))
    {
      if (any(grep("has new level", result)) || any(grep("New levels",result)))
        infoDialog(sprintf(Rtxt("It seems that the dataset on which the probabilities",
                                "from the %s model are required has a categoric",
                                "variable with levels not found in the training",
                                "dataset. The probabilities can not be determined in",
                                "this situation. You may need to either ensure",
                                "the training dataset has representatives of all levels",
                                "or else remove them from the testing dataset.",
                                "Alternatively, do not include that variable in the",
                                "modelling.\n\nThe actual error message was:"),
                           mtype),
                   "\n\n", paste(result, "\n"))
      else
        errorReport(probcmd, result)
      next()
    }

    eval(parse(text=plot.cmd))
  }
  ## If just one model, and we are plotting the test dataset, then
  ## also plot the training dataset.

  if (nummodels==1 && length(grep("\\[test\\]", testname))>0)
  {
    mcount <- mcount + 1
    plot.cmd <- paste(Rtxt("\n#In ROCR (1.0-3) plot does not obey the add command.\n"),
                      Rtxt("# Calling the function directly (.plot.performance) does work.\n\n"),
                      ".plot.performance(performance(prediction(crs$pr, ",
                      sprintf("%s$%s),",
                              sub("-crs\\$sample", "crs$sample",
                                  testset[[mtype]]), crs$target),
                      '"sens", "spec"), ',
                      'col="#00CCCCFF", lty=2, ',
                      sprintf("add=%s)\n", addplot),
                      sep="")
    appendLog(sprintf(Rtxt("Generate a Lift Chart for the %s model on %s."),
                     mtype, sub('\\[test\\]', '[train]', testname)),
             sub("-crs\\$sample", "crs$sample",
                                   probcmd[[mtype]]), "\n", plot.cmd)

    result <- try(eval(parse(text=sub("-crs\\$sample",
                               "crs$sample", probcmd[[mtype]]))), silent=TRUE)
    eval(parse(text=plot.cmd))
    models <- c(Rtxt("Test"), Rtxt("Train"))
    nummodels <- 2
    legtitle <- getEvaluateModels()
    title <- sub('\\[test\\]', '', testname)
  }
  else
  {
    models <- getEvaluateModels()
    legtitle <- Rtxt("Models")
    title <- testname
  }

  legendcmd <- paste('legend("bottomleft",',
                     sprintf("c(%s),",
                             paste('"', models, '"',
                                   sep="", collapse=",")),
                     sprintf('col=rainbow(%d, 1, .8), lty=1:%d,',
                             nummodels, nummodels),
                     sprintf('title="%s", inset=c(0.05, 0.05))', legtitle))
  appendLog(Rtxt("Add a legend to the plot."), legendcmd)
  eval(parse(text=legendcmd))

  decor.cmd <- paste(genPlotTitleCmd(Rtxt("Sensitivity/Specificity (tpr/tnr)"), "",
                                    title),
                    '\ngrid()', sep="")
  appendLog(Rtxt("Add decorations to the plot."), decor.cmd)
  eval(parse(text=decor.cmd))

  return(sprintf(Rtxt("Generated Sensitivity/Specificity Plot on %s."), testname))
}

#-----------------------------------------------------------------------
# SCORE

executeEvaluateScore <- function(probcmd, respcmd, testset, testname, dfedit.done=FALSE)
{
  # Apply each selected model to the selected dataset and save the
  # results to a file with columns containing the score (or scores in
  # the case of a multinomial model) from a specific model. Other
  # columns depend on the radio button options, and will either be
  # just the identifiers, or a copy of the full data, or else, the
  # score columns are written to the original file (assuming CSV).
  # TODO: Would this be better as the Export functionality for the
  # Evaluate tab? 081227 Add cluster export in here.

  # 100306 Allow data to be entered manually, and score that.

  entered <- theWidget("evaluate_enterdata_radiobutton")$getActive()

  if (entered & ! dfedit.done)
  {
    if (packageIsAvailable("RGtk2Extras"))
    {
      # 100307 Not quite ready yet - needs to know when to continue
      # after the data has been editted. Tom is fixing this up for
      # RGtk2Extras so I will be able to start using it then.
#      infoDialog(Rtxt ("RGtk2Extras will be used to edit",
#                      "a data frame called 'rattle.entered.dataset'. Once you have",
#                      "edited the dataset and the window is closed the dataset",
#                      "will be scored."))
      dsname <- "rattle.entered.dataset"
      if (exists(dsname))
        rattle.edit.obj <- RGtk2Extras::dfedit(rattle.entered.dataset, size=c(800, 400))
      else
        rattle.edit.obj <- RGtk2Extras::dfedit(crs$dataset[nrow(crs$dataset),
                                              c(crs$ident, crs$input, crs$target)],
                                  size=c(800, 400), dataset.name=dsname)

      probcmd <- lapply(probcmd, function(x) sub("crs\\$dataset", dsname, x))
      respcmd <- lapply(respcmd, function(x) sub("crs\\$dataset", dsname, x))
      testset <- lapply(testset, function(x) sub("crs\\$dataset", dsname, x))
      testname <- "manually entered data"

      RGtk2::gSignalConnect(rattle.edit.obj, "unrealize", data=rattle.edit.obj,
                     function(obj, data)
                     {
                       #print(paste("Exited", data$getDatasetName()))
                       #print(data$getDataFrame())

                       # 121210 As of now, remove the assign to global
                       # env - it is a bad idea and against CRAN
                       # policy, and iritates nasty riples. No
                       # solution for now - remove the functionality -
                       # rather minor anyhow.

                       # assign(dsname, data$getDataFrame(), envir=.GlobalEnv)
                       executeEvaluateScore(probcmd, respcmd, testset, testname, TRUE)
                       setStatusBar(Rtxt("Scored manually entered data."))
                     })

      setStatusBar(Rtxt("Enter the data in the editor and then close",
                        "it to have the data scored."))
      return()
    }
    else
    {
      if (! is.null(crs$entered))
        crs$entered <- edit(crs$entered) # Use previously manually entered data.
      else
        crs$entered <- edit(crs$dataset[nrow(crs$dataset),
                                        c(crs$ident, crs$input, crs$target)])
      probcmd <- lapply(probcmd, function(x) sub("crs\\$dataset", "crs$entered", x))
      respcmd <- lapply(respcmd, function(x) sub("crs\\$dataset", "crs$entered", x))
      testset <- lapply(testset, function(x) sub("crs\\$dataset", "crs$entered", x))
      testname <- "manually entered data"
    }
  }

  # Obtain information from the interface: what other data is to be
  # included with the scores.

  sinclude <- NULL
  if (theWidget("score_idents_radiobutton")$getActive())
    sinclude <- "idents"
  else if (theWidget("score_all_radiobutton")$getActive())
    sinclude <- "all"

  if (! entered)
  {

    # Obtain the filename to write the scores to.  We ask the user for a
    # filename if RATTLE_SCORE and RATTLE.SCORE.OUT are not provided.
    # TODO should we add getwd() to the RATTLE_SCORE or
    # RATTLE.SCORE.OUT if a relative path.

    fname <- Sys.getenv("RATTLE_SCORE")
    if (fname == "" && exists("RATTLE.SCORE.OUT") && not.null(RATTLE.SCORE.OUT))
      fname <- RATTLE.SCORE.OUT

    if (fname == "")
    {
      # The default filename is the testname with spaces replaced by
      # "_", etc., and then "_score" is appended, and then "_all" or
      # "_idents" to indicate what other columns are included.

      default <- sprintf("%s_score_%s.csv",
                         gsub(" ", "_",
                              gsub("\\.[[:alnum:]]*", "",
                                   gsub("(\\[|\\])", "",
                                        gsub("\\*", "", testname)))),
                         sinclude)
      # fname <- paste(getwd(), default, sep="/")

      dialog <- RGtk2::gtkFileChooserDialog(Rtxt("Score Files"), NULL, "save",
                                     "gtk-cancel", RGtk2::GtkResponseType["cancel"],
                                     "gtk-save", RGtk2::GtkResponseType["accept"])
      dialog$setDoOverwriteConfirmation(TRUE)

      if(not.null(testname)) dialog$setCurrentName(default)

      #dialog$setCurrentFolder(crs$dwd) Generates errors.

      ff <- RGtk2::gtkFileFilterNew()
      ff$setName(Rtxt("CSV Files"))
      ff$addPattern("*.csv")
      dialog$addFilter(ff)

      ff <- RGtk2::gtkFileFilterNew()
      ff$setName(Rtxt("All Files"))
      ff$addPattern("*")
      dialog$addFilter(ff)

      if (dialog$run() == RGtk2::GtkResponseType["accept"])
      {
        fname <- dialog$getFilename()
        dialog$destroy()
      }
      else
      {
        dialog$destroy()
        return()
      }
    }
  }

  # Score the data with each model, collect the outputs, and then
  # write them all at once to file.
  #
  # Note that there is at least one testset (hence, below we look at
  # just the first testset), but possibly others, and there is an
  # assumption that they are all of the forms:
  #
  # crs$dataset[-crs$sample, c(...)]
  # na.omit(crs$dataset[-crs$sample, c(...)])
  #
  # or else they are all of the forms:
  #
  # crs$testset[,c(...)]
  # na.omit(crs$testset[,c(...)])
  #
  # 080425 TODO It would be best to test they are all of the same form
  # to make sure the assumption is not breached.
  #
  # We first remove the na.omit so we can get all row names. The
  # na.omit is there for those models, like glm and ksvm, which do not
  # handle NA's themselves.

  ts <- testset[[1]]
  if (substr(ts, 1, 7) == "na.omit") ts <- sub('na.omit\\((.*)\\)$', '\\1', ts)

  # Create the data frame to hold the scores, initialised to NA in
  # every cell.

  the.names <- eval(parse(text=sprintf("row.names(%s)", ts)))
  the.models <- getEvaluateModels()
  scores <- as.data.frame(matrix(nrow=length(the.names),
                                 ncol=length(the.models)))
  row.names(scores) <- the.names
  names(scores) <- the.models

  # Obtain a list of the identity variables to output. 080713 Include
  # the target to output. 100531 Don't put the target in if it is not
  # in the dataset (like when we read a CSV file that does not contain
  # the target variable).

  if (length(grep("\\.csv$", testname)) &&
      ! getSelectedVariables("target") %in% names(crs$testset))
    idents <- getSelectedVariables("ident")
  else
    idents <- union(getSelectedVariables("ident"), getSelectedVariables("target"))

  for (mtype in the.models)
  {
    setStatusBar(sprintf(Rtxt("Scoring dataset using %s ..."), mtype))

    # Determine whether we want the respcmd (for trees and multinom)
    # or the probcmd (for logistic regression). 081204 Originally we
    # returned probabilities for glm and class for everything
    # else. But users want one or the other, so add a radiobutton
    # option to choose Class or Probability.
    if (theWidget("score_probability_radiobutton")$getActive())
    # if (mtype == crv$GLM) # This was RStat original approach.
      thecmd <- probcmd
    else
      thecmd <- respcmd

    # 110306 For Japanese, when we have a Japanese filename (in UTF-8)
    # if I don't do the following then the sprintf in the following
    # appendLog fails. I don't really know why this fixes it?

    if (isJapanese()) testname <- iconv(testname, from="UTF-8")

    # Apply the model to the dataset.

    appendLog(sprintf(Rtxt("Obtain %s for the %s model on %s."),
                      ifelse(mtype %in% c("kmeans", "hclust"),
                             Rtxt("cluster number"),
                             ifelse(categoricTarget(),
                                    Rtxt("probability scores"),
                                    ifelse(numericTarget(),
                                           Rtxt("predictions"),
                                           Rtxt("class")))),
                      commonName(mtype), testname),
              thecmd[[mtype]])

    result <- try(eval(parse(text=thecmd[[mtype]])), silent=TRUE)

    # Check for errors - in particular, new levels in the testset. If
    # an error is found we skip this mtype and proceed to the
    # next. This will leave NA's in the score file for this mtype.

    if (inherits(result, "try-error"))
    {
      if (any(grep("has new level", result)) || any(grep("New levels",result)))
        infoDialog(sprintf(Rtxt("The dataset on which the %s",
                                "model is applied to has a categoric",
                                "variable with levels not found in the training",
                                "dataset. The model can not be applied in",
                                "this situation. You may need to either ensure",
                                "the training dataset has representatives of all levels",
                                "or else remove them from the testing dataset.",
                                "Alternatively, do not include that variable in the",
                                "modelling.\n\nThe actual error message was:\n\n"),
                           mtype),
                   paste(result, "\n"))
      else
        errorReport(thecmd, result)
      next()
    }

    # 080417 Communicate the score file name. Note that originally I
    # intended to export the user's choice as an environment variable
    # to communicate that back to a calling process. But setenv
    # unfortunately does not export the name outside the R process so
    # it is of no use. TODO We could get a bit more sophisticated
    # here and add getwd() to the RATTLE_SCORE if it is a relative
    # path.

    # Transform the dataset expression into what we need to extract
    # the relevant columns.
    #
    # Various formats include:
    #
    #    train	crs$dataset[crs$sample, c(3:12,14)]
    #    test	crs$dataset[-crs$sample, c(3:12,14)]
    #    csv	crs$testset[,c(3:12,14)]
    #    df	crs$testset
    #
    # Want
    #    subset(crs$dataset[-crs$sample,], select=Idents) + crs$pr
    #

    scoreset <- testset[[mtype]]

    # If no comma in scoreset, leave as is, else find first comma,
    # remove everything after, and replace with "]". PROBLEM TODO If
    # the testset[[crv$MODEL]] includes na.omit, we need to do something
    # different because after the following step of replacing the
    # column list with nothing, it is very likely that new columns
    # are included that have NAs, and hence the na.omit will remove
    # even more rows for the subset command than it does for the
    # predict command. Yet we still want to ensure we have all the
    # appropriate columns available. So use
    # na.omit(crs$dataset[-crs$sample, c(2:4,6:10,13)])@na.action to
    # remove the rows from crs$dataset[-crs$sample,] that have
    # missing values with regard the columns c(2:4,6:10,13). Thus if
    # we have scoreset as:
    #
    #  na.omit(crs$dataset[-crs$sample, c(2:4,6:10,13)])
    #
    # we want to:
    #
    #  omitted <- na.omit(crs$dataset[-crs$sample, c(2:4,6:10,13)])@na.action
    #
    # and then scoreset should become:
    #
    #  crs$dataset[-crs$sample,][-omitted,]

    # First deal with the na.omit case, to capture the list of rows
    # omitted.

    # Mod by Ed Cox (080301) to fix error when no NAs in test set. The
    # error was:
    #
    #  Error in eval(expr, envir, enclos) :
    #    no slot of name "na.action" for this object of class "data.frame"
    #
    # The issue appears to be that it complains if there are no NAs to
    # remove so bypass the omission code.  This is not a general fix
    # because Ed has run into a case where it failed, but until we
    # have an insight as to what the real problem is we go with
    # this. It may be related to regression versus classification. Ed
    # was doing a regression (and testing the Pr v Ob plots).

    #scorevarlist <- c(getSelectedVariables("ident"),
    #                  getSelectedVariables("target"),
    #                  getSelectedVariables("input"),
    #                  getSelectedVariables("risk"))

    omitted <- NULL
    if (substr(scoreset, 1, 7) == "na.omit")
    {
      narm.dim <- eval(parse(text=sprintf("dim(%s)", scoreset)))[1]
      orig.dim <- eval(parse(text=sub('na.omit', 'dim', scoreset)))[1]
      if (narm.dim != orig.dim)

        # Ed had: && !dim(tmpset)[1]==dim(na.omit(tmpset))[1])

        # End of Ed's modification.
        # if (substr(scoreset, 1, 7) == "na.omit")
      {
        omit.cmd <- paste("omitted <- attr(", scoreset, ', "na.action")',
                          sep="")
        appendLog(Rtxt("Record rows omitted from predict command."), omit.cmd)
        eval(parse(text=omit.cmd))
      }
    }

    # Add the scores into the scores variable.

    if (is.null(omitted))
      scores[[mtype]] <- result
    else
      scores[[mtype]][-omitted] <- result
##  }
  }

  # Generate the other columns to be included in the score file.

  # Ensure we have all columns available in the dataset to start with,
  # so remove the " c(....)" selector from ts. We are then going to
  # include the identifiers or all columns in the output (depending on
  # the value of sinclude) so select those columns.

  if (length(grep(",", ts)) > 0) ts <- gsub(",.*]", ",]", ts)

  if (sinclude == "all")
    scoreset <- ts
  else if (sinclude == "idents")
    scoreset <- sprintf('subset(%s, select=c(%s))', ts,
                        ifelse(!length(idents), "",
                               sprintf('"%s"',
                                       paste(idents, collapse='", "'))))
  else
  {
    errorDialog(sprintf(Rtxt("We should not be here! The value of sinclude should have",
                             "been one of all or idents. We found: %s\n\n"),
                        sinclude), crv$support.msg)
    return()
  }

  appendLog(Rtxt("Extract the relevant variables from the dataset."),
            sprintf("sdata <- %s", scoreset))

  sdata <- eval(parse(text=scoreset))

  # 081107 Special case: for multinom, multiple probs are saved, plus
  # the decision. But the decision as a dataframe used to become the
  # column "glm.glm". I used to change that to "glm" but using a cbind
  # results in a column with no name, so things work out for
  # multinom. The cbinding seems to play a role when we have multiple
  # models? 081226 For rpart needed to use a data frame, and so need
  # to hanle the rpart.rpart column label problem.

  scores <- as.matrix(scores)
  gcol <- which(colnames(scores) == "glm.")
  if (length(gcol)) colnames(scores)[gcol] <- "glm"
  rcol <- grep("rpart.", colnames(scores))
  if (length(rcol)) colnames(scores)[rcol[length(rcol)]] <- "rpart"

  if (entered)
  {
    resetTextview("score_textview",
                  Rtxt("Scores for the manually entered data.\n\n"),
                  collectOutput("cbind(sdata, scores)", envir=environment()))
    return.msg <- Rtxt("The manually entered data has been scored.")
  }
  else
  {
    # 110306 Fix for Japanese
    if (isJapanese()) fname <- iconv(fname, from="UTF-8")

    appendLog(Rtxt("Output the combined data."),
              "write.csv(cbind(sdata, crs$pr), ",
              sprintf('file="%s", ', fname),
              "row.names=FALSE)")
    writeCSV(cbind(sdata, scores), file=fname)
    return.msg <- sprintf(Rtxt("Scores have been saved to the file %s"), fname)
  }

  # StatusBar is enough so don't pop up a dialog?
  # infoDialog("The scores have been saved into the file", fname)

  return(return.msg)
}

#-----------------------------------------------------------------------
# PR VERSUS OB PLOT

executeEvaluatePvOplot <- function(probcmd, testset, testname)
{
  ## print(probcmd)

  # This modification to executeEvaluateSave was provided by Ed Cox
  # (080201) to plot predictions vs. observed values. Graham added the
  # logging and some fine tuning. It includes a pseudo R-squared.

  # Put 1 or 2 charts onto their own plots. Otherwise, put the
  # multiple charts onto one plot, keeping them all the same size
  # (thus if numplots is odd, leave a cell of the plot empty.

  # TODO [140623] This needs to be updaed to work wiht crs$train, etc
  # rather than crs$sample.
  
  model.list <- getEvaluateModels()
  numplots <- length(model.list)

  if (numplots == 1)
    newPlot(1)
  else if (numplots == 2)
    newPlot(1)
  else if (numplots %% 2 == 0)
    newPlot(numplots)
  else
    newPlot(numplots + 1)

  if (numplots <= 2 )
    cex <- 1.0
  else if (numplots <= 4)
    cex <- 0.5
  else
    cex <- 0.5

  opar <- par(cex=cex)

  for (mtype in model.list)
  {

    appendLog(sprintf(Rtxt("%s: Generate a Predicted v Observed plot",
                            "for %s model on %s."),
                      toupper(mtype), mtype, testname),
              probcmd[[mtype]])

    result <- try(eval(parse(text=probcmd[[mtype]])), silent=TRUE)

    # Check for errors - in particular, new levels in the test
    # dataset. TODO This should be factored into a separate function,
    # since it is used in a number of places, including
    # executeEvaluateSave.

    if (inherits(result, "try-error"))
    {
      if (any(grep("has new level", result)) || any(grep("New levels",result)))
        infoDialog(sprintf(Rtxt("It seems that the dataset on which the probabilities",
                                "from the %s model are required has a categoric",
                                "variable with levels not found in the training",
                                "dataset. The probabilities can not be determined in",
                                "this situation. You may need to either ensure",
                                "the training dataset has representatives of all levels",
                                "or else remove them from the testing dataset.",
                                "Alternatively, do not include that variable in the",
                                "modelling.\n\nThe actual error message was:"),
                           mtype),
                   "\n\n", paste(result, "\n"))
      else
        errorReport(probcmd, result)
      next()
    }

    # Obtain a list of the identity variables.

    idents <- getSelectedVariables("ident")

    # Transform the dataset expression into what we need to extract
    # the relevant columns.
    #
    # TODO This should be factored into a separate function, since it
    # is used in a number of places, including executeEvaluateSave.
    #
    #
    # Various formats include:
    #
    #    train	crs$dataset[crs$sample, c(3:12,14)]
    #    test	crs$dataset[-crs$sample, c(3:12,14)]
    #    csv	crs$testset[,c(3:12,14)]
    #    df	crs$testset
    #
    # Want
    #    subset(crs$dataset[-crs$sample,], select=Idents) + crs$pr
    #

    scoreset <- testset[[mtype]]

    # If no comma in scoreset, leave as is, else find first comma,
    # remove everything after, and replace with "]". PROBLEM TODO If
    # the testset[[crv$MODEL]] includes na.omit, we need to do something
    # different because after the following step of replacing the
    # column list with nothing, it is very likely that new columns are
    # included that have NAs, and hence the na.omit will remove even
    # more rows for the subset command than it does for the predict
    # command. Yet we still want to ensure we have all the appropriate
    # columns available. So use na.omit(crs$dataset[-crs$sample,
    # c(2:4,6:10,13)])@na.action to remove the rows from
    # crs$dataset[-crs$sample,] that have missing values with regard
    # the columns c(2:4,6:10,13). Thus if we have scoreset as:
    #
    #  na.omit(crs$dataset[-crs$sample, c(2:4,6:10,13)])
    #
    # we want to:
    #
    #  omitted <- na.omit(crs$dataset[-crs$sample, c(2:4,6:10,13)])@na.action
    #
    # and then scoreset should become:
    #
    #  crs$dataset[-crs$sample,][-omitted,]

    # First deal with the na.omit case, to capture the list of rows
    # omitted.

    # Mod by Ed Cox (080301) to fix error when no NAs in test set. The
    # error was:
    #
    #  Error in eval(expr, envir, enclos) :
    #    no slot of name "na.action" for this object of class "data.frame"
    #
    # The issue appears to be that it complains if there are no NAs to
    # remove so bypass the omission code.  This is not a general fix
    # because Ed has run into a case where it failed, but until we
    # have an insight as to what the real problem is we go with
    # this. It may be related to regression versus classification. Ed
    # was doing a regression (and testing the Pr v Ob plots).

    scorevarlist <- c(getSelectedVariables("ident"),
                      getSelectedVariables("target"),
                      getSelectedVariables("input"),
                      getSelectedVariables("risk"))
    if (is.null(crs$sample))
      tmpset <- crs$dataset[, scorevarlist]
    else
      tmpset <- crs$dataset[-crs$sample, scorevarlist]

    if (substr(scoreset, 1, 7) == "na.omit" &&
        !dim(tmpset)[1]==dim(na.omit(tmpset))[1])
      # End of Ed's modification.
      # if (substr(scoreset, 1, 7) == "na.omit")
    {
      omit.cmd <- paste("omitted <- attr(", scoreset, ', "na.action")', sep="")
      appendLog(Rtxt("Record rows omitted from predict command."), omit.cmd)
      eval(parse(text=omit.cmd))
    }
    else
      omitted <- NULL

    # Now clean out the column subsets. [140623] Big fixing the case
    # where there are no missing amongst the crs$input but thre are
    # amongst the other variables. I need to check again why we remove
    # the variable list here - it seems correct to be
    # c(crs$input,crs$target) when we run the obs <- commands
    # below... So store the old value for now. Otherwise we get a
    # mismatch in the number of crs$pr and the target values. However,
    # the modified scoreset is not used for anything else so there
    # must have been a reason to remove it?

    orig.scoreset <- scoreset

    if (length(grep(",", scoreset)) > 0)
      scoreset <- gsub(",.*]", ",]", scoreset)

    # And finally, remove the na.omit if there is one, replacing it
    # with specifically removing just the rows that were removed in
    # the predict command.

    if (not.null(omitted))
      scoreset <- sub(")", "[-omitted,]", sub("na.omit\\(", "", scoreset))
    else
      # [140623] Restore the original version here to see if that
      # fixes the bug, without interferring with the complicated code
      # to handle the omitted.
      scoreset <- orig.scoreset

    # Extract the actual (i.e., observed) values that are to be
    # compared to the probabilities (i.e., predictions) from the
    # model.

    obsset <- sprintf('subset(%s, select=crs$target)', scoreset)
    appendLog(Rtxt("Obtain the observed output for the dataset."),
              paste("obs <-", obsset))
    obs <- eval(parse(text=obsset))

    # 090506 Account for the case when a categoric is being treated as
    # a numeric.

    obsfix <- sprintf(paste("obs.rownames <- rownames(obs)\n",
                            "obs <- as.numeric(obs[[1]])\n",
                            "obs <- data.frame(%s=obs)\n",
                            "rownames(obs) <- obs.rownames", sep=""),
                      crs$target)
    appendLog(Rtxt("Handle in case categoric target treated as numeric."),
              obsfix)
    eval(parse(text=obsfix))

    # fitcorr is the so called psuedo-R-square. It has a maximum less
    # than 1 and is often used in either binary or multinomial
    # logistic regression. This is to be interpreted differently to
    # the standard R-square.

    fit.cmd <- sprintf("na.omit(cbind(obs, %s=crs$pr))", Rtxt("Predicted"))
    appendLog(Rtxt("Combine the observed values with the predicted."),
              paste("fitpoints <-", fit.cmd))
    fitpoints <- eval(parse(text=fit.cmd))

    corr.cmd <- "format(cor(fitpoints[,1], fitpoints[,2])^2, digits=4)"
    appendLog(Rtxt("Obtain the pseudo R2 - a correlation."),
              paste("fitcorr <-", corr.cmd))
    fitcorr <- eval(parse(text=corr.cmd))

    # Plot the points - observed versus predicted.

    # For 2 plots, so as not to overwrite the first plot, if we are
    # about to plot the second plot, initiate a new plot.

    if (numplots == 2 && mtype == model.list[length(model.list)]) newPlot(1)

    par.cmd <- 'par(c(lty="solid", col="blue"))'
    appendLog(Rtxt("Plot settings for the true points and best fit."),
              paste("op <-", par.cmd))
    op <- eval(parse(text=par.cmd))

    # In the plot I originally limited the x and y to (0,1). Not sure
    # why needed. Ed Cox pointed out he was losing values when
    # predicting more than (0,1) (linear regression), so remove the limits
    # for now (080301).

    # 130322 For cforest, the Predicted column gets to be named the
    # same as the target variable! Unlike for other models. Something
    # about the structure of crs$pr So make sure it is called
    # Predicted here.
    
    # vnames <- names(fitpoints)
    vnames <- c(names(fitpoints)[1], "Predicted")
    plot.cmd <-sprintf('plot(%s, fitpoints[[2]], asp=1, xlab="%s", ylab="%s")',
                       ifelse(length(unique(fitpoints[[1]])) < crv$max.categories,
                              "jitter(fitpoints[[1]])",
                              "fitpoints[[1]]"),
                       ifelse(length(unique(fitpoints[[1]])) < crv$max.categories,
                              paste(vnames[1], Rtxt("(Jittered)")),
                              vnames[1]),
                       vnames[2])
    appendLog(Rtxt("Display the observed (X) versus predicted (Y) points."),
              plot.cmd)
    eval(parse(text=plot.cmd))

    # Fit a linear model Predicted ~ Observed.

    lm.cmd <- paste("lm(fitpoints[,2] ~ fitpoints[,1])")
    appendLog(Rtxt("Generate a simple linear fit between predicted and observed."),
              paste("prline <-", lm.cmd))
    prline <- eval(parse(text=lm.cmd))

    ab.cmd <- "abline(prline)"
    appendLog(Rtxt("Add the linear fit to the plot."),
              ab.cmd)
    eval(parse(text=ab.cmd))

    diag.cmd <- paste('par(c(lty="dashed", col="black"))',
                      'abline(0, 1)', sep="\n")
    appendLog(Rtxt("Add a diagonal representing perfect correlation."),
              diag.cmd)
    eval(parse(text=diag.cmd))

    legend("topleft", legend=c(Rtxt("Linear Fit to Points"),
                        Rtxt("Predicted=Observed")),
           lty=c(1, 2),col=c("blue", "black"), bty="n")

    par(op)

    # 100206 Move the pseudo-R squared to the plot itself.

    legend.cmd <- paste('legend("bottomright", ',
                        'sprintf("', Rtxt("Pseudo R-square=%s"), '", fitcorr), ',
                        'bty="n")', sep=" ")
    appendLog(Rtxt("Include a pseudo R-square on the plot"), legend.cmd)
    eval(parse(text=legend.cmd))

    # TODO Add to LOG

    # Add decorations. 100206 Rado suggested not including the
    # filename, but probably still want an indication of whether it is
    # a test or training dataset that is part of the testname - but
    # need to extract that.

    title.cmd <- paste(genPlotTitleCmd(paste(Rtxt("Predicted vs. Observed"), "\n", sep=""),
                                       commonName(mtype),
                                       paste(Rtxt("Model"), "\n", sep=""), testname),
                       '\ngrid()', sep="")
    appendLog(Rtxt("Add a title and grid to the plot."), title.cmd)
    eval(parse(text=title.cmd))

  }
  return(Rtxt("Pr v Ob plot generated."))
}

executeEvaluateHand <- function(probcmd, testset, testname)
{
  # 130119 Replace with the HMeasure package functionality.

  lib.cmd <- "library(hmeasure)"
  if (! packageIsAvailable("hmeasure", "compute Hand's performance measures")) return()
  appendLog(Rtxt("David Hand's Performance Measures"), lib.cmd)
  # 140906 Move to using namespaces within the code, though still
  # expose the interactive commands.
  eval(parse(text=lib.cmd))
  
  model.list <- getEvaluateModels()
  numplots <- length(model.list)

  for (mtype in model.list)
  {
    newPlot(4)
    pr <- try(eval(parse(text=probcmd[[mtype]])), silent=TRUE)
    scoreset <- testset[[mtype]]
    obsset <- sprintf('subset(%s, select=crs$target)', scoreset)

    # 100408 Handle categoric target (assumed binary)
    # Replace
    #   obs <- eval(parse(text=obsset))
    # with the following two lines

    obs <- eval(parse(text=obsset))[[1]]
    if (is.factor(obs)) obs <- as.numeric(obs)-1

    # 130119 Use the now released package rather than David's original
    # code.

    result <- hmeasure::HMeasure(obs, pr)
    hmeasure::plotROC(result, which=1)
    hmeasure::plotROC(result, which=2)
    hmeasure::plotROC(result, which=3)
    hmeasure::plotROC(result, which=4)

    with(result$metric,
         cat(sprintf("%s: H=%f,  Gini=%f, AUC=%f,AUCH=%f, KS=%f\n",
                     mtype, H, Gini, AUC, AUCH, KS)))
  }
}
