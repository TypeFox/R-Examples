# Gnome R Data Miner: GNOME interface to R for Data Mining
#
# Time-stamp: <2015-07-12 16:32:36 gjw>
#
# RANDOM FOREST TAB
#
# Copyright (c) 2006-2011 Graham Williams, Togaware.com, GPL Version 2

#######################################################################
# CALLBACKS

setRFOptions <- function(mtype)
{
  theWidget("model_rf_sample_label")$setSensitive(mtype=="randomForest")
  theWidget("model_rf_sample_entry")$setSensitive(mtype=="randomForest")
  theWidget("model_rf_impute_checkbutton")$setSensitive(mtype=="randomForest")
}

on_model_rf_traditional_radiobutton_toggled <- function(button)
{
  if (button$getActive())
  {
    theWidget("model_rf_builder_label")$setText("randomForrest")
    setRFOptions("randomForest")
  }
}

on_model_rf_conditional_radiobutton_toggled <- function(button)
{
  if (button$getActive())
  {
    theWidget("model_rf_builder_label")$setText("cforest")
    setRFOptions("cforest")
  }
}

on_rf_importance_button_clicked <- function(button)
{
  plotRandomForestImportance()
}

on_rf_errors_button_clicked <- function(button)
{
  plotRandomForestError()
}

on_rf_oob_roc_button_clicked <- function(button)
{
  plotRandomForestOOBROC()
}

on_rf_print_tree_button_clicked <- function(button)
{
  displayRandomForestTree()
}

on_help_randomForest_activate <- function(action, window)
{
  if (showHelpPlus(Rtxt("The randomForest algorithm builds multiple",
                        "decision trees from different samples of the dataset, and while",
                        "building each tree, random subsets of the available variables are",
                        "considered for splitting the data at each node of the tree. A simple",
                        "majority vote is then used for prediction in the case of",
                        "classification (and average for regression).",
                        "RandomForest's are generally robust against over fitting.",
                        "<<>>",
                        "The default is to build 500 trees and to select the square root of the",
                        "number of variables as the subset to choose from at each node. The",
                        "resulting model is generally not very sensitive to the choice of these",
                        "parameters.",
                        "<<>>",
                        "Any observation with missing values will be ignored, which may lead to some",
                        "surprises, like many fewer observations to model when many missing values",
                        "exist. It can also lead to losing all examples of a particular class!",
                        "<<>>",
                        "An estimate of the error rate is provided as the out-of-bag (OOB)",
                        "estimate. This applies each tree to the data that was not used in",
                        "building the tree to give a quite accurate estimate of the error",
                        "rate.",
                        "<<>>",
                        "The Sample Size can be used to down-sample larger classes.",
                        "For a two-class problem with, for example, 5000 in class 0 and 250 in class 1,",
                        "a Sample Size of '250, 250' will usually give a more 'balanced' classifier.",
                        "<<>>",
                        "The R package for building Random Forests is called randomForest.")))
    {
      popupTextviewHelpWindow("randomForest", "randomForest")
    }
}

########################################################################
#
# MODEL -> RF -> Traditional - RANDOM FOREST
#

executeModelRF <- function(traditional=TRUE, conditional=!traditional)
{

  # Decide which random forest to build: randomForest or cforest
  
  FUN <- ifelse(traditional, "randomForest::randomForest", "party::cforest")

  # Make sure the appropriate package is available.
  
  if ((traditional
       && ! packageIsAvailable("randomForest", Rtxt("build a random forest model")))
      ||
      conditional &&
      ! packageIsAvailable("party", Rtxt("build a random forest model")))
    return(FALSE)

  # Identify the appropriate textview widget to display results.
  
  TV <- "rf_textview"

  # Build up the funciton call.
  
  num.classes <- length(levels(as.factor(crs$dataset[[crs$target]])))
  parms <- ""

  # For the traditional random forest make sure there are no included
  # variables which have more than 32 levels, which can not be handled
  # by randomForest.

  if (traditional)
  {
    categoricals <- crs$input[unlist(lapply(crs$input,
                                            function(x) is.factor(crs$dataset[,x])))]

    for (i in categoricals)
      if (length(levels(crs$dataset[,i])) > 32)
      {
        errorDialog(sprintf(Rtxt("This implementation of random forests does not handle",
                                 "categorical variables with more than 32 levels.",
                                 "Having a large number of levels tends to introduce",
                                 "bias into the tree algorithms. The variable %s has",
                                 "%d levels\n\nPlease choose to ignore it in the",
                                 "Data tab if you wish to build a random forest model."),
                            i, length(levels(crs$dataset[,i]))))
        return(FALSE)
      }
  }

  # Retrieve options and set up the function arguments.

  ntree <- theWidget("rf_ntree_spinbutton")$getValue()
  # 100107 Always use the supplied value rather than chekcing if it is
  # changed from the default, which would use:
  #
  # if (ntree != crv$rf.ntree.default)
  parms <- sprintf("\n      ntree=%d", ntree)
  
  mtry <- theWidget("rf_mtry_spinbutton")$getValue()
  # 100107 Always use the supplied value rather than chekcing if it is
  # changed from the default, which would use:
  #
  # if (mtry != crv$rf.mtry.default)
  parms <- sprintf("%s,\n      mtry=%d", parms, mtry)

  if (traditional)
  {
    sampsize <- theWidget("model_rf_sample_entry")$getText()
    if (nchar(sampsize) > 0)
    {
      ss <- as.numeric(unlist(strsplit(sampsize, ",")))
      if (! (length(ss) == 1 || length(ss) == num.classes))
      {
        errorDialog(sprintf(Rtxt("The supplied sample sizes (%s) need to be either",
                                 "a single number or else correspond to the number",
                                 "of classes found in the target variable '%s'.",
                                 "Please supply exactly 1 or %d sample sizes."),
                            sampsize, crs$target, num.classes))
        return(FALSE)
      }
      # TODO Check if sample sizes are larger than the classes!
      parms <- sprintf("%s,\n      sampsize=c(%s)", parms, sampsize)
    }
    
    # By default the MeanDecreaseGini is available for plotting. With
    # importance MeanDecreaseAccuracy is also available, and it seems
    # to also print the relative importance with regard class. So by
    # default, generate them both.
  
    parms <- sprintf("%s,\n      importance=TRUE", parms)
  }
  
  # Build the formula for the model. TODO We assume we will always do
  # classification rather than regression, at least for now.

  frml <- paste(ifelse(is.factor(crs$dataset[[crs$target]]),
                       crs$target,
                       sprintf(ifelse(numericTarget(), "%s", "as.factor(%s)"),
                               crs$target)),
                "~ .")

  # List, as a string of indicies, the variables to be included. 

  # included <- getIncludedVariables()
  included <- "c(crs$input, crs$target)" # 20110102
  
  # Some convenience booleans

  sampling <- not.null(crs$sample)
  including <- not.null(included)
  subsetting <- sampling || including

  # Ensure we have some data - i.e., not all records will be removed
  # because they have missing values.

  dataset <- paste("crs$dataset",
                   if (subsetting) "[",
                   if (sampling) "crs$sample",
                   if (subsetting) ",",
                   if (including) included,
                   if (subsetting) "]",
                   sep="")
  # 130322 Don't na.omit for cforest - not needed?
  #if (! traditional)
  #  dataset <- sprintf("na.omit(%s)", dataset)

  # Replicate rows according to the integer weights variable for
  # randomForest.
  
  if(! is.null(crs$weights))
    if (traditional)
      dataset <- paste(dataset,
                       "[rep(row.names(",
                       dataset,
                       "),\n                                        ",
                       # Use eval since crs$weights could be a formula
                       'as.integer(eval(parse(text = "', crs$weights,
                       '"))',
                       if (sampling) '[crs$sample]',
                       ')),]',
                       sep="")
    else
      dataset <- sprintf("%s,\n      weights=%s%s", dataset,
                       crs$weights,
                       ifelse(sampling, "[crs$sample]", ""))


  # 100107 Deal with missing values. I've not tested whether cforest
  # has issues with missing values.
  
  naimpute <- theWidget("model_rf_impute_checkbutton")$getActive()
  if (traditional)
  {
    missing.cmd <- sprintf('length(attr((na.omit(%s)), "na.action"))', dataset)
    result <- try(missing <- eval(parse(text=missing.cmd)), silent=TRUE)
    if (inherits(result, "try-error")) missing <- 0
    dsrow.cmd <- sprintf("nrow(%s)", dataset)
    result <- try(dsrow <- eval(parse(text=dsrow.cmd)), silent=TRUE)
    if (inherits(result, "try-error")) dsrow <- 0
    if (missing == dsrow && ! naimpute)
    {
      errorDialog(Rtxt("All observations in the dataset have one or more",
                       "missing values for the selected variables.",
                       "The random forest algorithm ignores any observation",
                       "with missing values. No observations remain from which",
                       "to build a model.",
                       "To fix this problem, you might, for example, Ignore any",
                       "variable with many or any missing values.",
                       "Or else enable the Impute check button to impute",
                       "the medium (numeric) or most frequent (categoric) value",
                       "using randomForests' na.roughfix().",
                       "You could also use imputation to fill in default or modelled",
                       "values for the missing values manually."))
      return()
    }
  }
    
  # Start the log.
  
  startLog(commonName(crv$RF))

  # Load the required package into the library.

  lib.cmd <- ifelse(traditional,
                    "library(randomForest, quietly=TRUE)",
                    "library(party, quietly=TRUE)")

  if (traditional)
    appendLog(packageProvides("randomForest", "randomForest"), lib.cmd)
  else
    appendLog(packageProvides("party", "cforest"), lib.cmd)
  
  eval(parse(text=lib.cmd))

  # Build the model.

  rf.cmd <- paste("set.seed(crv$seed)\n",
                  "crs$rf <- ", FUN, "(", frml,
                  ",\n      data=",
                  dataset, ", ",
#                  ifelse(traditional,
#                         ",\n      data=crs$dataset",
#                         ",\n      data=na.omit(crs$dataset"),
#                  if (subsetting) "[",
#                  if (sampling) "crs$sample",
#                  if (subsetting) ",",
#                  if (including) included,
#                  ifelse(subsetting,
#                         ifelse(traditional, "], ", "]), "),
#                         ifelse(traditional, "", ")")),
                  ifelse(traditional, parms,
                         sprintf("controls=cforest_unbiased(%s)", parms)),
                  ifelse(traditional,
                         sprintf(",\n      na.action=%s",
                                 ifelse(naimpute, "randomForest::na.roughfix", "na.omit")),
                         ""),
                  # 100521 Turn subsampling with replacement off since
                  # it is more likely to produce biased imprtance
                  # measures, as explained in
                  # http://www.ncbi.nlm.nih.gov/pmc/articles/PMC1796903/
                  # Note that for cforest, cforest_unbiased uses
                  # replace=FALSE also.
                  if (traditional) ",\n      replace=FALSE",
                  ")", sep="")

  appendLog(sprintf(Rtxt("Build the %s model."), commonName(crv$RF)), rf.cmd)

  start.time <- Sys.time()
  result <- try(eval(parse(text=rf.cmd)), silent=TRUE)

  if (inherits(result, "try-error"))
  {
    if (any(grep("cannot allocate vector", result)))
    {
      errorDialog(sprintf(Rtxt("The call to '%s' appears to have failed.",
                               "This is often due, as in this case,",
                               "to running out of memory.",
                               "A quick solution is to sample the dataset, through the",
                               "Transform tab. On 32 bit machines you may be limited to",
                               "less than 2000 observations."),
                          FUN))
      setTextview(TV)
    }
    else if (any(grep("NA/NaN/Inf", result)))
    {
      # TODO 080520 This error arose when a log transform is done on
      # Deductions where there are many 0's (hence -Inf). To be more
      # helpful, find the variable with the -Inf and suggest ignoring
      # it. We can test this is the error if the following returns
      # non-zero:
      #
      # sum(crs$dataset[crs$sample,c(2:10, 13:14)]==-Inf, na.rm=TRUE)
      
      errorDialog(sprintf(Rtxt("The call to '%s' failed.",
                               "The problem may be with the data",
                               "containing Infinite values.",
                               "A quick solution may be to remove variables",
                               "with any Inf or -Inf values."), FUN))
                          
      setTextview(TV)
    }
    else if (any(grep('inherits', result)))
    {
      # 080825 This is a known error in randomForest 4.5-25 reported
      # to Andy Liau a few months ago, and acknowledge, but has not
      # been fixed. Of course i could try to figure it out myself, but
      # it would probably take some effort!
      
      errorDialog(Rtxt("The call to randomForest failed.",
                       "You probably have version 4.5-25.",
                       "This is a known problem and is fixed in 4.5-26.",
                       "Please install a newer version of randomForest.\n",
                       "\ninstall.packages('randomForest',\n",
                       "    repos='http://rattle.togaware.com')"))
      setTextview(TV)
    }
    else 
      errorDialog(errorMessageFun(FUN, result))
    return(FALSE)
  }

  resetTextview(TV)
  
  # Display the resulting model.

  summary.cmd <- "crs$rf"
  appendLog(sprintf(Rtxt("Generate textual output of '%s' model."), commonName(crv$RF)),
            summary.cmd)

  addTextview(TV, sprintf(Rtxt("Summary of the %s Model"), commonName(crv$RF)),
              "\n==================================",
              "\n\n",
              "Number of observations used to build the model: ",
              ifelse(traditional, length(crs$rf$y), crs$rf@responses@nobs),
              "\n",
              ifelse(naimpute, "Missing value imputation is active.\n", ""),
              sub(") \\n", ")\n\n",
                  sub(", ntree", ",\n              ntree",
                      sub(", data", ",\n              data",
                          gsub(", *", ", ", collectOutput(summary.cmd, TRUE))))))

  # 6 Mar 2012 Add some AUC information as suggested by Akbar Waljee.

  if (traditional && binomialTarget() &&
      packageIsAvailable("pROC", Rtxt("calculate AUC confidence interval")))
  {
    roc.cmd <- "pROC::roc(crs$rf$y, as.numeric(crs$rf$predicted))"
    ci.cmd <- "pROC::ci.auc(crs$rf$y, as.numeric(crs$rf$predicted))"

    appendLog(Rtxt("The `pROC' package implements various AUC functions."))

    appendLog("Calculate the Area Under the Curve (AUC).", roc.cmd)
    appendLog("Calculate the AUC Confidence Interval.", ci.cmd)
    
    addTextview(TV,
                "\n\nAnalysis of the Area Under the Curve (AUC)",
                "  \n==========================================\n",
                collectOutput(roc.cmd),
                "\n\n",
                collectOutput(ci.cmd))
  }
  
  # Display the variable importance.

  # 100107 There is a very good importance measure from cforest that
  # needs to go here.
  
  if (traditional)
  {
    if (numericTarget())
      varimp.cmd <- paste("rn <- round(randomForest::importance(crs$rf), 2)",
                              "rn[order(rn[,1], decreasing=TRUE),]",
                              sep="\n")
    else
      varimp.cmd <- paste("rn <- round(randomForest::importance(crs$rf), 2)",
                          # 100521 Sort on the accuracy measure rather
                          # than the Gini measure, since the Gini is
                          # biased in favour of categoric variables
                          # with many categories.
                          # "rn[order(rn[,3] + rn[,4], decreasing=TRUE),]",
                          "rn[order(rn[,3], decreasing=TRUE),]",
                          sep="\n")
  }
  else
  {
    varimp.cmd <- "data.frame(Importance=sort(varimp(crs$rf), decreasing=TRUE))"
  }

  appendLog(Rtxt("List the importance of the variables."), varimp.cmd)

  result <- try(collectOutput(varimp.cmd), silent=TRUE)
  if (inherits(result, "try-error"))
  {
    msg <- errorMessageFun(ifelse(traditional, "randomForest::importance", "party::varimp"), result)
    errorDialog(msg)
    return(FALSE)
  }
  
  addTextview(TV, sprintf("\n\n%s", Rtxt("Variable Importance")),
              "\n===================\n\n", result, "\n")

  # 100107 What is the purpose of this?

  if (sampling) crs$smodel <- union(crs$smodel, crv$RF)

  # Now that we have a model, make sure the buttons are sensitive.

  showModelRFExists(traditional=traditional, conditional=conditional)

  # Finish up.

  time.taken <- Sys.time()-start.time

  reportTimeTaken(TV, time.taken, commonName(crv$RF))

  return(TRUE)
}

showModelRFExists <- function(traditional=TRUE, conditional=!traditional)
{
  # If an rf model exists then show the various buttons on the Model
  # tab.

  state=!is.null(crs$rf)

  if (state)
  {
    theWidget("rf_importance_button")$show()
    theWidget("rf_importance_button")$setSensitive(TRUE)
    theWidget("rf_errors_button")$show()
    theWidget("rf_errors_button")$setSensitive(TRUE)
    theWidget("rf_oob_roc_button")$show()
    theWidget("rf_oob_roc_button")$setSensitive(TRUE)
    theWidget("rf_print_tree_button")$show()
    theWidget("rf_print_tree_button")$setSensitive(TRUE)
    theWidget("rf_print_tree_spinbutton")$show()
    theWidget("rf_print_tree_spinbutton")$setSensitive(TRUE)
  }
  else
  {
    theWidget("rf_importance_button")$hide()
    theWidget("rf_errors_button")$hide()
    theWidget("rf_oob_roc_button")$hide()
    theWidget("rf_print_tree_button")$hide()
    theWidget("rf_print_tree_spinbutton")$hide()
  }
}

plotRandomForestImportance <- function()
{

  # Make sure there is an rf object first.

  if (is.null(crs$rf))
  {
    errorDialog(Rtxt("E123: This is an unexpected error.",
                     "There is no RF and attempting to plot importance."),
                crv$support.msg)
    return()
  }

  if ( "RandomForest" %in% class(crs$rf))
    if (!packageIsAvailable("ggplot2", Rtxt("display the conditional forest var importance")))
      return()
  
  newPlot()
  if ("RandomForest" %in% class(crs$rf))
    plot.cmd <- paste(#'library(ggplot2)',
                      #'',
                      '# Note that values vary slightly on each call to varimp().',
                      '',
                      'v    <- party::varimp(crs$rf)',
                      'vimp <- data.frame(Variable=as.character(names(v)),',
                      '                   Importance=v,',
                      '                   row.names=NULL, stringsAsFactors=FALSE)',
                      '',
                      'p <- ggplot2::ggplot(vimp, ggplot2::aes(Variable, Importance)) +',
                      '  ggplot2::geom_bar(stat="identity", position="identity") +',
                      '  ggplot2::scale_x_discrete(limits=with(vimp, Variable[order(Importance)])) +',
                      '  ggplot2::theme(legend.position="none",',
                      '                 axis.title.x = ggplot2::element_blank(),',
                      '                 axis.title.y = ggplot2::element_blank()) +',
                      '  ggplot2::coord_flip()',
                      'print(p)',
                      sep="\n")
  else
    plot.cmd <- paste('randomForest::varImpPlot(crs$rf, main="")\n',
                      genPlotTitleCmd(Rtxt("Variable Importance"),
                                      commonName(crv$RF), crs$dataname),
                      sep="")

  startLog()
  appendLog(Rtxt("Plot the relative importance of the variables."),
            sub("print\\(p\\)", "p", plot.cmd))

  set.cursor("watch")
  on.exit(set.cursor())

  eval(parse(text=plot.cmd))

  setStatusBar(Rtxt("Variable Importance has been plotted."))
}
  
plotRandomForestError <- function()
{
  # Make sure there is an rf object first.

  if (is.null(crs$rf))
  {
    errorDialog(Rtxt("E129: This is an unexpected error.",
                     "There is no RF and attempting to plot errors."),
                crv$support.msg)
    return()
  }
  
  newPlot()
  plot.cmd <- paste('plot(crs$rf, main="")',
                    sprintf(paste('legend("topright", %s, text.col=1:6,',
                                  'lty=1:3, col=1:3)'),
                            sprintf("c(%s)", paste('"', colnames(crs$rf$err.rate),
                                                   '"', sep="", collapse=", "))),
                    genPlotTitleCmd(Rtxt("Error Rates"), commonName(crv$RF), crs$dataname),
                    sep="\n")

  appendLog(Rtxt("Plot the error rate against the number of trees."), plot.cmd)
  eval(parse(text=plot.cmd))
  
  setStatusBar(Rtxt("The error rates plot has been generated."))
}

plotRandomForestOOBROC <- function()
{
  # Make sure there is an rf object first.

  if (is.null(crs$rf))
  {
    errorDialog(Rtxt("There is no RF and attempting to plot OOB ROC."),
                crv$support.msg)
    return()
  }

  
  if (! packageIsAvailable("verification", Rtxt("plot ROC curve")))
    return()

  newPlot()

  # 111115 Akbar Waljee suggested using roc.plot from
  # verification. Would also like a confidence interval but not sure
  # how to. 130129 Akbar noted that when Impute button is disabled
  # then the OOB plot fails. It fails because ommitted observations
  # for the model are not being ommitted from the attempted plot.
  # Need to add na.omit around the dataset conditional on naimpute.
  
  naimpute <- theWidget("model_rf_impute_checkbutton")$getActive()
  plot.cmd <- paste('library(verification)',
                    # 111116 the as.integer as.factor to ensure
                    # 0/1. Need to fix to ensure button is only
                    # avilable for binary classification.
                    '\naucc <- verification::roc.area(as.integer(as.factor(',
                    if (!naimpute) 'na.omit(',
                    'crs$dataset[',
                    if (not.null(crs$sample)) 'crs$sample',
                    if (!naimpute) ',])[',
                    ', crs$target]))-1,',
                    '\n                 crs$rf$votes[,2])$A',
                    '\nverification::roc.plot(as.integer(as.factor(',
                    if (!naimpute) 'na.omit(',
                    'crs$dataset[',
                    if (not.null(crs$sample)) 'crs$sample',
                    if (!naimpute) ',])[',
                    ', crs$target]))-1,',
                    '\n         crs$rf$votes[,2], main="")',
                    '\nlegend("bottomright", bty="n",',
                    '\n       sprintf("Area Under the Curve (AUC) = %1.3f", aucc))\n',
                    genPlotTitleCmd(Rtxt("OOB ROC Curve"),
                                    commonName(crv$RF), crs$dataname),
                    sep="")

  appendLog(Rtxt("Plot the OOB ROC curve."), plot.cmd)
  eval(parse(text=plot.cmd))
  
  setStatusBar(Rtxt("The OOB ROC curve has been generated."))
}

displayRandomForestTree <- function()
{
  # Initial setup.
  
  TV <- "rf_textview"

  # Obtain which tree to display.
  
  tree.num <- theWidget("rf_print_tree_spinbutton")$getValue()

  # If tree.num is zero and there are very many trees, first warn.

  if (tree.num == 0 && crs$rf$ntree > 5)
    if (! questionDialog(Rtxt("Displaying all rules can take quite a while.",
                              "\n\nDo you wish to continue?")))
      return()
  
  # Command to run.

  display.cmd <- ifelse(class(crs$rf) == "RandomForest",
                        sprintf(paste('party:::prettytree(crs$rf@ensemble[[%d]],',
                                      'names(crs$rf@data@get("input")))'), tree.num),
                        sprintf("printRandomForests(crs$rf, %d)", tree.num))

  # Perform the action.

  appendLog(sprintf(Rtxt("Display tree number %d."), tree.num), display.cmd)
  set.cursor("watch")
  setStatusBar(Rtxt("The rules are being generated ..."))
  addTextview(TV, collectOutput(display.cmd, TRUE), textviewSeparator())
  set.cursor()
  setStatusBar(paste(ifelse(tree.num == 0,
                            Rtxt("Rules from all trees have been added to the textview."),
                            sprintf(Rtxt("Rules from tree %d have been added to the textview."),
                                    tree.num)),
                     Rtxt("You may need to scroll the textview to view them.")))
}

printRandomForests <- function(model, models=NULL, include.class=NULL,
                               format="")
{
  # format=
  #    "VB"	Generate code that looks like VisualBasic
  
  if (! packageIsAvailable("randomForest", Rtxt("print the rule sets")))
    return()

  if (is.null(models) || models == 0) models <- 1:model$ntree

  for (i in models)
    printRandomForest(model, i, include.class, format)
}

## Move to using the more generic functions

treeset.randomForest <- function(model, n=1, root=1, format="R")
{
  ## Return a string listing the decision tree form of the chosen tree
  ## from the random forest.
  
  tree <- randomForest::getTree(model, n)
  if (format == "R")
  {
    cassign <- "<-"
    cif <- "if"
    cthen <- ""
    celse <- "else"
    cendif <- ""
    cin <- "%in%"
  }
  else if (format == "VB")
  {
    cassign <- "="
    cif <- "If"
    cthen <- "Then"
    celse <- "Else"
    cendif <- "End If"
    cin <- "In"
  }

  ## Traverse the tree

  tr.vars <- attr(model$terms, "dataClasses")[-1]
  var.names <- names(tr.vars)
  
  result <- ""
  if (tree[root, 'status'] == -1) # Terminal node
  {
    result <- sprintf("Result %s %s", cassign,
                      levels(model$y)[tree[root,'prediction']])
  }
  else
  {
    var.class <- tr.vars[tree[root, 'split var']]
    node.var <- var.names[tree[root,'split var']]
    if("character" %in% var.class | "factor" %in% var.class | "ordered" %in% var.class)
    {
      ## Convert the binary split point to a 0/1 list for the levels.
      
      var.levels <- levels(eval(model$call$data)[[tree[root,'split var']]])
      bins <- sdecimal2binary(tree[root, 'split point'])
      bins <- c(bins, rep(0, length(var.levels)-length(bins)))
      node.value <- var.levels[bins==1]
      node.value <- sprintf('("%s")', paste(node.value, collapse='", "'))
      condition <- sprintf("%s %s %s%s", node.var, cin,
                           ifelse(format=="R", "c", ""), node.value)
    }
    else if ("integer" %in% var.class | "numeric" %in% var.class)
    {
      ## Assume spliting to the left means "<=", and right ">",
      ## which is not what the man page for getTree claims!

      node.value <- tree[root, 'split point']
      condition <- sprintf("%s <= %s", node.var, node.value)

    }
    else
    {
      stop(sprintf("Rattle E222: getRFRuleSet: class %s not supported.",
                   var.class))
    }
    

    condition <- sprintf("%s (%s)", cif, condition)
    
    lresult <- treeset.randomForest(model, n, tree[root,'left daughter'],
                                    format=format)
    if (cthen == "")
      lresult <- c(condition, lresult)
    else
      lresult <- c(condition, cthen, lresult)
    rresult <- treeset.randomForest(model, n, tree[root,'right daughter'],
                                    format=format)
    rresult <- c(celse, rresult)
    result <- c(lresult, rresult)
    if (cendif != "") result <- c(result, cendif)
  }
  return(result)
}

ruleset.randomForest <- function(model, n=1, include.class=NULL)
{

  ## NOT YET WORKING
  
  ## Same as printRandomForest, for now, but returns the string rather
  ## than printing it. Perhaps it is a list of vectors of strings,
  ## each in the list is one rule, and each string in the vector
  ## represents a single test. I HAVE SINCE FINISHED
  ## treeset.randomForest AND PERHAPS CAN USE THAT AS A TEMPLATE.
  
  # include.class	Vector of predictions to include
  
  if (! packageIsAvailable("randomForest", Rtxt("generate a rule set")))
    stop(Rtxt("the 'randomForest' package is required to generate rule sets"))

  if (!inherits(model, "randomForest"))
    stop(Rtxt("the model is not of the 'randomForest' class"))

  tr <- randomForest::getTree(model, n)
  tr.paths <- getRFPathNodesTraverse(tr)
  tr.vars <- attr(model$terms, "dataClasses")[-1]

  ruleset <- list()
  nameset <- c()
  
  ## Generate a simple form for each rule.

  for (i in seq_along(tr.paths))
  {
    tr.path <- tr.paths[[i]]
    nodenum <- as.integer(names(tr.paths[i]))
    target <- levels(model$y)[tr[nodenum,'prediction']]

    if (! is.null(include.class) && target %notin% include.class) next()

    rule <- c()
    
    #cat(sprintf("%sTree %d Rule %d Node %d Decision %s\n\n",
    #            comment, n, i, nodenum, target))

    #nrules <- nrules + 1
    
    ## Indicies of variables in the path
    
    var.index <- tr[,3][abs(tr.path)] # 3rd col is "split var"
    var.names <- names(tr.vars)[var.index]
    var.values <- tr[,4][abs(tr.path)] # 4th col is "split point"
    
    for (j in 1:(length(tr.path)-1))
    {
      var.class <- tr.vars[var.index[j]]
      if (var.class == "character" | var.class == "factor" | var.class == "ordered")
      {
        node.op <- "IN"

        ## Convert the binary to a 0/1 list for the levels.
        
        var.levels <- levels(eval(model$call$data)[[var.names[j]]])
        bins <- sdecimal2binary(var.values[j])
        bins <- c(bins, rep(0, length(var.levels)-length(bins)))
        if (tr.path[j] > 0)
          node.value <- var.levels[bins==1]
        else
          node.value <- var.levels[bins==0]
        node.value <- sprintf('%s', paste(node.value, collapse=', '))
      }
      else if (var.class == "integer" | var.class == "numeric")
      {
        ## Assume spliting to the left means "<=", and right ">",
        ## which is not what the man page for getTree claims!
        if (tr.path[j]>0)
          node.op <- "<="
        else
          node.op <- ">"
        node.value <- var.values[j]
      }
      else
        stop(sprintf("Rattle E223: getRFRuleSet: class %s not supported.",
                     var.class))

      rule <- c(rule, sprintf("%s %s %s\n", var.names[j], node.op, node.value))
    }
    ruleset <- c(ruleset, rule)
    nameset <- c(nameset, nodenum)
    break()
  }
  names(ruleset) <- nameset
  return(ruleset)
}


printRandomForest <- function(model, n=1, include.class=NULL,
                              format="", comment="")
{
  # include.class	Vector of predictions to include
  
  if (! packageIsAvailable("randomForest", Rtxt("generate the rule sets")))
    return()

  if (!inherits(model, "randomForest"))
    stop(Rtxt("the model is not of the 'randomForest' class"))

  if (format=="VB") comment="'"
  
  tr <- randomForest::getTree(model, n)
  tr.paths <- getRFPathNodesTraverse(tr)
  tr.vars <- attr(model$terms, "dataClasses")[-1]
  
  ## Initialise the output

  cat(sprintf("%sRandom Forest Model %d", comment, n), "\n\n")

  ## Generate a simple form for each rule.

  cat(paste(comment,
            "-------------------------------------------------------------\n",
            sep=""))

  if (format=="VB")
    cat("IF FALSE THEN\n' This is a No Op to simplify the code\n\n")
  
  ## Number of rules generated

  nrules <- 0
  
  for (i in seq_along(tr.paths))
  {
    tr.path <- tr.paths[[i]]
    nodenum <- as.integer(names(tr.paths[i]))
    # 090925 This needs work to make it apply in the case of a
    # regression model. For now simply note this in the output.
    target <- levels(model$y)[tr[nodenum,'prediction']]

    if (! is.null(include.class) && target %notin% include.class) next()
    
    cat(sprintf("%sTree %d Rule %d Node %d %s\n \n",
                comment, n, i, nodenum,
                ifelse(is.null(target), "Regression (to do - extract predicted value)",
                       paste("Decision", target))))

    if (format=="VB") cat("ELSE IF TRUE\n")

    nrules <- nrules + 1
    
    ## Indicies of variables in the path
    
    var.index <- tr[,3][abs(tr.path)] # 3rd col is "split var"
    var.names <- names(tr.vars)[var.index]
    var.values <- tr[,4][abs(tr.path)] # 4th col is "split point"
    
    for (j in 1:(length(tr.path)-1))
    {
      var.class <- tr.vars[var.index[j]]
      if (var.class == "character" | var.class == "factor" | var.class == "ordered")
      {
        node.op <- "IN"

        ## Convert the binary to a 0/1 list for the levels.
        
        var.levels <- levels(eval(model$call$data)[[var.names[j]]])
        bins <- sdecimal2binary(var.values[j])
        bins <- c(bins, rep(0, length(var.levels)-length(bins)))
        if (tr.path[j] > 0)
          node.value <- var.levels[bins==1]
        else
          node.value <- var.levels[bins==0]
        node.value <- sprintf('("%s")', paste(node.value, collapse='", "'))
      }
      else if (var.class == "integer" | var.class == "numeric")
      {
        ## Assume spliting to the left means "<", and right ">=",
        ## which is not what the man page for getTree claims!
        if (tr.path[j]>0)
          node.op <- "<="
        else
          node.op <- ">"
        node.value <- var.values[j]
      }
      else
        stop(sprintf("Rattle E234: getRFRuleSet: class %s not supported.",
                     var.class))

      if (format=="VB")
        cat(sprintf("AND\n%s %s %s\n", var.names[j], node.op, node.value))
      else
        cat(sprintf("%d: %s %s %s\n", j, var.names[j], node.op, node.value))
    }
    if (format=="VB") cat("THEN Count = Count + 1\n")
    cat("-----------------------------------------------------------------\n")
  }
  if (format=="VB") cat("END IF\n\n")
  cat(sprintf("%sNumber of rules in Tree %d: %d\n\n", comment, n, nrules))
}

randomForest2Rules <- function(model, models=NULL)
{
  if (! packageIsAvailable("randomForest", Rtxt("generate the rule sets")))
    return()

  if (is.null(models)) models <- 1:model$ntree

  ## Obtain information we need about the data
  
  vars <- attr(model$terms, "dataClasses")[-1]

  ruleset <- list()
  
  for (i in models)
  {
    ruleset[[i]] <- list(ruleset=getRFRuleSet(model, i))
  }
  return(ruleset)
}

## Generate a list of rules (conditions and outcomes) for RF MODEL
## number N.

getRFRuleSet <- function(model, n)
{
  if (! packageIsAvailable("randomForest", Rtxt("generate the rule sets")))
    return()

  tr <- randomForest::getTree(model, n)
  tr.paths <- getRFPathNodes(tr)
  tr.vars <- attr(model$terms, "dataClasses")[-1]
  
  ## Initialise the output
  
  rules <- list()

  ## Generate rpart form for each rule.

  for (i in seq_along(tr.paths))
  {
    tr.path <- tr.paths[[i]]

    ## Indicies of variables in the path
    
    var.index <- tr[,3][abs(tr.path)] # 3rd col is "split var"
    var.names <- names(tr.vars)[var.index]
    var.values <- tr[,4][abs(tr.path)] # 4th col is "split point"
    
    tr.rule <- c("root")

    for (j in seq_along(tr.path))
    {
#print(j)
#print(tr.vars)
#print(var.index)
      var.class <- tr.vars[var.index[j]]
#print(var.class)
      # 6 Mar 2012 This is currently failing. Needs debugging.
      if (var.class == "character" | var.class == "factor" | var.class == "ordered")
      {
        node.op <- "="

        ## Convert the binary to a 0/1 list for the levels.
        
        var.levels <- levels(eval(model$call$data)[[var.names[j]]])
        bins <- sdecimal2binary(var.values[j])
        bins <- c(bins, rep(0, length(var.levels)-length(bins)))
        if (tr.path[j] > 0)
          node.value <- var.levels[bins==1]
        else
          node.value <- var.levels[bins==0]
        node.value <- paste(node.value, collapse=",")
      }
      else if (var.class == "integer" | var.class == "numeric")
      {
        ## Assume spliting to the left means "<", and right ">=",
        ## which is not what the man page for getTree claims!
        if (tr.path[j]>0)
          node.op <- "<="
        else
          node.op <- ">"
        node.value <- var.values[j]
      }
      else
        stop(sprintf("Rattle E235: getRFRuleSet: class %s not supported.",
                     var.class))

      tr.rule <- c(tr.rule, paste(var.names[j], node.op, node.value))
    }
    
    ## TODO Walk through tr.rule and remove all but the last "VAR<="
    ## and "VAR>" conditions.
    
    rules[[i]] <- list(rule=tr.rule)
  }
  return(rules)
}

getRFPathNodesTraverse <- function(tree, root=1)
{
  # Traverse the paths through the binary decision tree represented as
  # a matrix.
  #
  # The columns in the RF tree matrix are:
  #   1. left daughter;
  #   2. right daughter;
  #   3. split var;
  #   4. split point;
  #   5. status;
  #   6. prediction
  
  paths <- list()
  if (tree[root,'status'] == -1) # Terminal node
  {
    paths <- list(root)
    names(paths) <- root
  }
  else
  {
    lpaths <- getRFPathNodesTraverse(tree, tree[root,'left daughter'])
    lpaths <- lapply(lpaths, append, root, 0)
    rpaths <- getRFPathNodesTraverse(tree, tree[root,'right daughter'])
    rpaths <- lapply(rpaths, append, -root, 0)
    paths <- c(lpaths, rpaths)
  }
  return(paths)
}

getRFPathNodes <- function(tree.matrix)
{

  ## I am now using Traverse as above, as it is a more logical
  ## ordering of the rules. This function can disappear some day.
  
  # The columns in the RF tree matrix are:
  #   1. left daughter;
  #   2. right daughter;
  #   3. split var;
  #   4. split point;
  #   5. status;
  #   6. prediction
  #
  # The traversal is essentially from the shortest paths to the longer
  # paths, since we simply list the leaf nodes in the order they
  # appear in the matrix, and then list the rules in this same
  # order. For each rule we work backwards to build the path. An
  # alternative might be to order the rules in a recursive left
  # traversal manner, in which case we get a more logical ordering to
  # the rules.

  # Number of nodes in the tree
  
  nnodes <- dim(tree.matrix)[1] 

  # Leaf node indices: leaf (or terminal) nodes have a value of -1 for
  # the fifth column (status) of the tree matrix. 

  lnodes <- which(tree.matrix[,5] == -1)
  
  # Initial the paths, being a list, each element is a vector of the
  # indices of a path from the root to the leaf.
  
  paths <- list() 

  # Process each leaf's path back to the root node.
  
  for (i in seq_along(lnodes))
  {
    # Initialise the node to the leaf index
    
    node <- lnodes[[i]]
    pathI <- node
    repeat
    {
      leftV <- 1:nnodes * as.integer(tree.matrix[,1]==abs(node))
      leftNode <- leftV[leftV!=0]
      if (length(leftNode)!= 0)
      {
        node <- leftNode
      }
      else # else must not be in the next line
      {
        rightV <- 1:nnodes * as.integer(tree.matrix[,2]==abs(node))
        node <- -rightV[rightV!=0] # rhs node identified with negative index
      }

      pathI <-c(node, pathI)

      ## If the node is the root node (first row in the matrix), then
      ## the path is complete.
      
      if (abs(node)==1) break
    }
    paths[[i]] <- pathI
  }

  ## Each path is named after its leaf index number
  
  names(paths) <- as.character(lnodes)

  return(paths)
}

sdecimal2binary <- function(x)
{
  return(rev(sdecimal2binary.smallEndian(x)))
}
	
sdecimal2binary.smallEndian <- function(x)
{
  if (x==0) return(0)
  if (x<0) stop(Rtxt("the input must be positive"))
  dec <- x
	 
  n <- floor(log(x)/log(2))
  bin <- c(1)
  dec <- dec - 2 ^ n
	  
  while(n > 0)
  {
    if (dec >= 2 ^ (n-1)) {bin <- c(bin,1); dec <- dec - 2 ^ (n-1)}
    else bin <- c(bin,0)
    n <- n - 1
  }
  return(bin)
}
	
exportRandomForestModel <- function()
{
  # Make sure we have a model first!

  if (noModelAvailable(crs$rf, crv$RF)) return(FALSE)

  startLog(paste(Rtxt("Export"), commonName(crv$RF)))

  save.name <- getExportSaveName(crv$RF)
  if (is.null(save.name)) return(FALSE)
  ext <- tolower(get.extension(save.name))

  # Construct the command to produce PMML. Currently
  # dataset=/transforms= are not needed/supported. TODO.

  pmml.cmd <- sprintf("pmml(crs$rf%s, dataset=crs$dataset)",
                      ifelse(length(crs$transforms),
                             ", transforms=crs$transforms", ""))

  # We can't pass "\" in a filename to the parse command in MS/Windows
  # so we have to run the save/write command separately, i.e., not
  # inside the string that is being parsed.

  if (ext == "xml")
  {
    appendLog(sprintf(Rtxt("Export %s as PMML."), commonName(crv$RF)),
              sprintf('saveXML(%s, "%s")', pmml.cmd, save.name))
    XML::saveXML(eval(parse(text=pmml.cmd)), save.name)
  }
  else if (ext == "c")
  {
    # 090103 gjw Move to a function: saveC(pmml.cmd, save.name,
    # "random forest")

    # 090223 Why is this tolower being used? Under GNU/Linux it is
    # blatantly wrong. Maybe only needed for MS/Widnows

    if (isWindows()) save.name <- tolower(save.name)

    model.name <- sub("\\.c", "", basename(save.name))

    export.cmd <- generateExportPMMLtoC(model.name, save.name, "rf_textview")
    
    appendLog(sprintf(Rtxt("Export %s as a C routine."), commonName(crv$RPART)),
              sprintf('pmml.cmd <- "%s"\n\n', pmml.cmd),
              export.cmd)

    eval(parse(text=export.cmd))
  }      
  setStatusBar(sprintf(Rtxt("The model has been exported to '%s'."), save.name))
}
