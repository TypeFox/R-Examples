# Gnome R Data Miner: GNOME interface to R for Data Mining
#
# Time-stamp: <2015-05-17 08:57:27 gjw>
#
# MODEL TAB
#
# Copyright (c) 2009 Togaware Pty Ltd
#
# This file is part of Rattle.
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

# When radio button is selected, display appropriate tab page

on_model_linear_radiobutton_toggled <- function(button)
{
  if (button$getActive())
  {
    crv$MODEL$setCurrentPage(crv$MODEL.GLM.TAB)
    # 090222 Do we really want to reset the textview? Obviously I
    # decided to do so some time ago, so continut to do so, but call
    # resetTextviews for this textview).
    ## setTextview("confusion_textview")
    resetTextviews("confusion_textview")
  }
  setStatusBar()
}

on_dtree_radiobutton_toggled <- function(button)
{
  if (button$getActive())
  {
    crv$MODEL$setCurrentPage(crv$MODEL.RPART.TAB)
    ## setTextview("confusion_textview")
    resetTextviews("confusion_textview")
  }
  setStatusBar()
}

on_boost_radiobutton_toggled <- function(button)
{
  if (button$getActive())
  {
    ## crv$MODEL$setCurrentPage(crv$MODEL.GBM.TAB)
    crv$MODEL$setCurrentPage(crv$MODEL.ADA.TAB)
    ## setTextview("confusion_textview")
    resetTextviews("confusion_textview")
  }
  setStatusBar()
}

on_nnet_radiobutton_toggled <- function(button)
{
  if (button$getActive())
  {
    crv$MODEL$setCurrentPage(crv$MODEL.NNET.TAB)
    ## setTextview("confusion_textview")
    resetTextviews("confusion_textview")
  }
  setStatusBar()
}

on_rf_radiobutton_toggled <- function(button)
{
  if (button$getActive())
  {
    crv$MODEL$setCurrentPage(crv$MODEL.RF.TAB)
    ## setTextview("confusion_textview")
    resetTextviews("confusion_textview")
  }
  setStatusBar()
}

on_svm_radiobutton_toggled <- function(button)
{
  if (button$getActive())
  {
    crv$MODEL$setCurrentPage(crv$MODEL.SVM.TAB)
    ## setTextview("confusion_textview")
    resetTextviews("confusion_textview")
  }
  setStatusBar()
}

on_e1071_radiobutton_toggled <- function(button)
{
  if (button$getActive())
  {
    crv$SVMNB$setCurrentPage(crv$SVMNB.ESVM.TAB)
  }
  setStatusBar()
}

on_kernlab_radiobutton_toggled <- function(button)
{
  if (button$getActive())
  {
    crv$SVMNB$setCurrentPage(crv$SVMNB.KSVM.TAB)
  }
  setStatusBar()
}

on_model_survival_radiobutton_toggled <- function(button)
{
  if (button$getActive())
  {
    crv$MODEL$setCurrentPage(crv$MODEL.SURVIVAL.TAB)
    # 090222 Do we really want to reset the textview? Obviously I
    # decided to do so some time ago, so continut to do so, but call
    # resetTextviews for this textview).
    ## setTextview("confusion_textview")
    resetTextviews("confusion_textview")
  }
  setStatusBar()
}

on_model_linear_plot_button_clicked <- function(button)
{
  # Make sure there is an appropriate model first.

  if (is.null(crs$glm))
  {
    errorDialog("There is no GLM and attempting to plot it.", crv$support.msg)
    return()
  }
  newPlot(4)
  plot.cmd <- paste('ttl <- genPlotTitleCmd("Linear Model",crs$dataname,vector=TRUE)\n',
                    'plot(crs$glm, main=ttl[1])',
                    sep="")
  appendLog("Plot the model evaluation.", plot.cmd)
  eval(parse(text=plot.cmd))

  setStatusBar("Linear model evaluation has been plotted.")
}

#-----------------------------------------------------------------------
# Model -> Tree
#
# Set the model builder label appropriately.

setTreeOptions <- function(mtype)
{
  theWidget("model_tree_priors_label")$setSensitive(mtype=="rpart")
  theWidget("model_tree_priors_entry")$setSensitive(mtype=="rpart")
  theWidget("model_tree_loss_label")$setSensitive(mtype=="rpart")
  theWidget("model_tree_loss_entry")$setSensitive(mtype=="rpart")
  theWidget("model_tree_cp_label")$setSensitive(mtype=="rpart")
  theWidget("model_tree_cp_spinbutton")$setSensitive(mtype=="rpart")
  theWidget("model_tree_include_missing_checkbutton")$
  setSensitive(mtype=="rpart")
}
  

on_model_tree_rpart_radiobutton_toggled <- function(button)
{
  if (button$getActive())
  {
    theWidget("model_tree_builder_label")$setText("rpart")
    setTreeOptions("rpart")
  }
}

on_model_tree_ctree_radiobutton_toggled <- function(button)
{
  if (button$getActive())
  {
    theWidget("model_tree_builder_label")$setText("ctree")
    setTreeOptions("ctree")
  }
}

# Model -> Linear
#
# When any of the regression radion buttons change then ensure the
# Model Builder label is updated to indicate the right model bulder.

on_glm_linear_radiobutton_toggled <- function(button)
{
  if (button$getActive())
    theWidget("model_linear_builder_label")$setText("lm")
}

on_glm_gaussian_radiobutton_toggled <- function(button)
{
  if (button$getActive())
    theWidget("model_linear_builder_label")$setText("glm (gaussian)")
}

on_model_linear_poisson_radiobutton_toggled <- function(button)
{
  if (button$getActive())
    theWidget("model_linear_builder_label")$setText("glm (poisson)")
}

on_glm_logistic_radiobutton_toggled <- function(button)
{
  if (button$getActive())
    theWidget("model_linear_builder_label")$setText("glm (logit)")
}

on_model_linear_probit_radiobutton_toggled <- function(button)
{
  if (button$getActive())
    theWidget("model_linear_builder_label")$setText("glm (probit)")
}

on_glm_multinomial_radiobutton_toggled <- function(button)
{
  if (button$getActive())
    theWidget("model_linear_builder_label")$setText("multinom")
}

on_evaluate_model_checkbutton_toggled <- function(button)
{
  configureEvaluateTab()
}

########################################################################
# UTILITIES

commonName <- function(mtype)
{
  name.map <- data.frame(ada=Rtxt("Ada Boost"),
                         arules=Rtxt("Association Rules"),
                         biclust=Rtxt("BiCluster"),
                         cforest=Rtxt("Random Forest"),
                         ctree=Rtxt("Conditional Tree"),
                         ewkm=Rtxt("Entropy Weighted KMeans"),
                         hclust=Rtxt("Hierarchical"),
                         kmeans=Rtxt("KMeans"),
                         rf=Rtxt("Random Forest"),
                         rpart=Rtxt("Decision Tree"),
                         svm=Rtxt("SVM"),
                         ksvm=Rtxt("SVM"),
                         glm=Rtxt("Linear"),
                         linear=Rtxt("Linear"),
                         multinom=Rtxt("Neural Net"),
                         nnet=Rtxt("Neural Net"),
                         survival=Rtxt("Survival"))
  return(as.character(name.map[[mtype]]))
}

numericTarget <- function()
{
  if (length(getSelectedVariables("target")) == 0)
    return(FALSE)

  # 091206 Move to not using the auto radio button.091223 Gerry and
  # Rado want it back
  
  else if (theWidget("data_target_auto_radiobutton")$getActive())

    # 080505 TODO we should put 10 as a global CONST
    
    return(is.numeric(crs$dataset[[crs$target]]) &&
           length(levels(as.factor(crs$dataset[[crs$target]]))) > 10)

  else if (theWidget("data_target_categoric_radiobutton")$getActive())
    return(FALSE)
  else if (theWidget("data_target_numeric_radiobutton")$getActive())
    return(TRUE)
  else
    return(FALSE)
  
}

countTarget <- function()
{
  # 080913 Return TRUE if the target is numeric, and is or looks like
  # an integer (i.e., no decimal points), and is all positive.
  
  return(numericTarget()
         && (is.integer(crs$dataset[[crs$target]])
             || ! any(grep("\\.", crs$dataset[[crs$target]])))
         && ! any(crs$dataset[[crs$target]] < 0, na.rm=TRUE))
}

categoricTarget <- function()
{
  if (length(getSelectedVariables("target")) == 0)
    return(FALSE)

  # 091206 Move to not using the auto radio button. 091223 Gerry and
  # Rado want it back
  
  else if (theWidget("data_target_auto_radiobutton")$getActive())
    
    # 080505 TODO we should put 10 as a global CONST

    return(is.factor(crs$dataset[[crs$target]]) ||
           (is.numeric(crs$dataset[[crs$target]]) &&
            length(levels(as.factor(crs$dataset[[crs$target]]))) <= 10))

  else if (theWidget("data_target_categoric_radiobutton")$getActive())
    return(TRUE)
  else if (theWidget("data_target_numeric_radiobutton")$getActive())
    return(FALSE)
  else
    return(FALSE)
}

binomialTarget <- function()
{
  return(categoricTarget() &&
         length(levels(as.factor(crs$dataset[[crs$target]]))) == 2)
}

multinomialTarget <- function()
{
  return(categoricTarget() &&
         length(levels(as.factor(crs$dataset[[crs$target]]))) > 2)
}

survivalTarget <- function()
{
  return(theWidget("data_target_survival_radiobutton")$getActive())
}

currentModelTab <- function()
{
  lb <- getCurrentPageLabel(crv$MODEL)
  if (lb == crv$SVM && theWidget("kernlab_radiobutton")$getActive()) lb <- crv$KSVM
  return(lb)
}

existsCategoricModel <- function()
{
  # TRUE if there is a classification model - i.e., a model to predict
  # a class.
  return(categoricTarget() && existsPredictiveModel())
}

existsPredictiveModel <- function()
{
  # TRUE if there is a predictive model as distinct from a descriptive
  # model.
  return(! is.null(listBuiltModels(c(crv$KMEANS, crv$HCLUST, crv$APRIORI))))
}

noModelAvailable <- function(model, model.class)
{
  # 090812 A general test for a model existing. Designed for use for
  # each of the export<model>Tab functions. We need to also know the
  # model.class since if model is NULL we won't be able to determine
  # it.

  if (is.null(model))
    errorDialog("No", commonName(model.class),
                "model is available. Be sure to build",
                "the model before trying to export it! You will need",
                "to press the Execute button (F2) in order to build the",
                "model.")
  return(is.null(model))
}

########################################################################
# EXECUTE MODEL TAB

executeModelTab <- function()
{
  # Perform the actions requested from the Model tab.
  
  # Check for prerequisites: Can not build a model without a dataset.

  if (noDatasetLoaded()) return()

  # If VARIABLES has some ignores but crs$ignore is NULL, complain.

  if (variablesHaveChanged(Rtxt("building a model"))) return()

  # If the WeightCalculator has changed but it is not the same as
  # crs$weight, complain. This doesn't work any more since we add
  # crs$dataset to the variable names in the Weights Calculator, so
  # they are different! But, let's remove the crs$dataset and compare.

  weights.display <- gsub('crs\\$dataset\\$', '', crs$weights)

  if (theWidget("weight_entry")$isSensitive()
      && not.null(crs$weights)
      && weights.display != theWidget("weight_entry")$getText())
  {
    errorDialog(sprintf(Rtxt("You appear to have changed the formula for calculating",
                             "the weights on the Data tab without executing the tab.",
                             "The previous formula was '%s' and it is now '%s'.",
                             "Please be sure to execute the Data tab",
                             "before continuing."),
                        crs$weights, theWidget("weight_entry")$getText()))
    return()
  }
    
  # Retrieve the target and make sure there is one.

  if (length(crs$target) == 0)
  {
    errorDialog(Rtxt("No target has been specified.",
                     "Please identify the target using the Data tab.",
                     "Be sure to Execute the tab once the target has",
                     "been identified."))
    return()
  }

  # Check if sampling needs executing.

  if (sampleNeedsExecute()) return()
    
  # If the target is a categorical and has more than 2 levels then
  # disable the ROCR and Risk plots, and place a message on the first
  # textview of the Evaluate tab. We make this word wrap here and then
  # turn that off once the tab is Executed.

  #091113 configureEvaluateTab()

  if (multinomialTarget())
  {
    theWidget("confusion_textview")$setWrapMode("word")
    resetTextview("confusion_textview")
    appendTextview("confusion_textview",
                   Rtxt("Note that the target you have chosen has more than",
                        "2 classes. Some functionality on the Evaluate tab",
                        "will not be available. In particular, the ROCR",
                        "package (Lift, ROC, Precision, and Sensitivity",
                        "charts) and the Risk Chart only handle binary",
                        "classification."))
  }
  else if (numericTarget())
  {
    resetTextviews("confusion_textview")
  }
  else
  {
    resetTextviews("confusion_textview")
  }

  # Dispatch

  build.all <- theWidget("all_models_radiobutton")$getActive()

  # Reset all Evaluate options to unchecked.

  #091112 resetEvaluateTab("all_inactive")
  #resetEvaluateTab(desensitise=FALSE)
  
  # The following work for ada, do they work for the rest?
  
  formula <- paste(crs$target, "~ .")
  included <- "c(crs$input, crs$target)" # 20110102 getIncludedVariables()
  sampling <- not.null(crs$train)
  including <- not.null(included)
  subsetting <- sampling || including
  dataset <- paste("crs$dataset",
                   if (subsetting) "[",
                   if (sampling) "crs$train",
                   if (subsetting) ",",
                   if (including) included,
                   if (subsetting) "]",
                   sep="")
  
  # This order of execution should correspond to the order in the
  # GUI as this makes most logical sense to the user.

  start.time <- Sys.time()

  if (build.all || currentModelTab() == crv$RPART)
  {
    if (theWidget("rpart_build_radiobutton")$getActive())
    {
      setStatusBar(sprintf(Rtxt("Building %s model ..."), commonName(crv$RPART)))

      if (theWidget("model_tree_ctree_radiobutton")$getActive())
      {
        if (executeModelCTree())
          theWidget("evaluate_rpart_checkbutton")$setActive(TRUE)
        else
          setStatusBar(sprintf(Rtxt("Building %s model ... failed."),
                               commonName("ctree")))
      }
      else
      {
        if (executeModelRPart())
            theWidget("evaluate_rpart_checkbutton")$setActive(TRUE)
        else
          setStatusBar(sprintf(Rtxt("Building %s model ... failed."),
                               commonName(crv$RPART)))
      }
    }
    else if (theWidget("rpart_tune_radiobutton")$getActive())
    {
      setStatusBar(sprintf(Rtxt("Tuning %s model ..."), commonName(crv$RPART)))
      if (! executeModelRPart("tune"))
        setStatusBar(sprintf(Rtxt("Tuning %s model ... failed"), commonName(crv$RPART)))
    }
    else if (theWidget("rpart_best_radiobutton")$getActive())
    {
      setStatusBar(sprintf(Rtxt("Building best %s model ..."), commonName(crv$RPART)))
      if (! executeModelRPart("best"))
        setStatusBar(sprintf(Rtxt("Building best %s model ... failed"),
                             commonName(crv$RPART)))
    }
    else # That's all the radio buttons - we should not be here.
    {
      errorDialog(Rtxt("Tried building an rpart model with option not",
                       "one of build/tune/best. This should not be possible."),
                  crv$support.msg)
      return(FALSE)
      
    }
  }
  if ((binomialTarget() && build.all)
      || currentModelTab() == crv$ADA)
  {
    setStatusBar(sprintf(Rtxt("Building %s model ..."), commonName(crv$ADA)))
    crs$ada <-
      buildModelAda(formula,
                    dataset,
                    tv=theWidget("ada_textview"),
                    maxdepth=theWidget("ada_maxdepth_spinbutton")$getValue(),
                    minsplit=theWidget("ada_minsplit_spinbutton")$getValue(),
                    cp=theWidget("ada_cp_spinbutton")$getValue(),
                    xval=theWidget("ada_xval_spinbutton")$getValue(),
                    ntree=theWidget("ada_ntree_spinbutton")$getValue())
    if (not.null(crs$ada))
    {
      showModelAdaExists()
      theWidget("evaluate_ada_checkbutton")$setActive(TRUE)
    }
    else
      setStatusBar(sprintf(Rtxt("Building %s model ... failed."), commonName(crv$ADA)))
  }

  if (build.all || currentModelTab() == crv$RF)
  {
    setStatusBar(sprintf(Rtxt("Building %s model ..."), commonName(crv$RF)))
    if (executeModelRF(traditional=theWidget("model_rf_traditional_radiobutton")$
                       getActive(),
                       conditional=theWidget("model_rf_conditional_radiobutton")$
                       getActive()))
      theWidget("evaluate_rf_checkbutton")$setActive(TRUE)
    else
      setStatusBar(sprintf(Rtxt("Building %s model ... failed."), commonName(crv$RF)))
  }
  
  if ((categoricTarget() && build.all)
      || currentModelTab() %in% c(crv$SVM, crv$KSVM))
  {
    setStatusBar(sprintf(Rtxt("Building %s model ..."), commonName(crv$KSVM)))
    if (executeModelSVM())
      theWidget("evaluate_ksvm_checkbutton")$setActive(TRUE)
    else
      setStatusBar(sprintf(Rtxt("Building %s model ... failed."), commonName(crv$KSVM)))

  }
  if (build.all || currentModelTab() == crv$GLM)
  {
    setStatusBar(sprintf(Rtxt("Building %s model ..."), commonName(crv$GLM)))
    if (executeModelGLM())
      theWidget("evaluate_glm_checkbutton")$setActive(TRUE)
    else
      setStatusBar(sprintf(Rtxt("Building %s model ... failed."), commonName(crv$GLM)))
  }
  if ((theWidget("nnet_radiobutton")$isSensitive() && build.all)
      || currentModelTab() == crv$NNET)
  {
    setStatusBar(sprintf(Rtxt("Building %s model ..."), commonName(crv$NNET)))
    if (executeModelNNet())
      theWidget("evaluate_nnet_checkbutton")$setActive(TRUE)
    else
      setStatusBar(sprintf(Rtxt("Building %s model ... failed."), commonName(crv$NNET)))
  }
  # 100416 The currentModelTab returns "Survival" for the
  # translations, but "survival" for English? Yet the above work? Need
  # to know why, but as a temporary fix, convert to lower here.
  if (tolower(currentModelTab()) == crv$SURVIVAL)
  {
    # 20110102 included <- getIncludedVariables(risk=TRUE)
    included <- "c(crs$input, crs$risk, crs$target)"
    including <- not.null(included)
    dataset <- paste("crs$dataset",
                     if (subsetting) "[",
                     if (sampling) "crs$train",
                     if (subsetting) ",",
                     if (including) included,
                     if (subsetting) "]",
                     sep="")
    setStatusBar(sprintf(Rtxt("Building %s model ..."), commonName(crv$SURVIVAL)))
      # ?? survfit(Surv(RISK_MM, as.numeric(RainTomorrow)) ~ MinTemp, data=crs$dataset)
    formula <- sprintf("Surv(%s, %s) ~ %s",
                       crs$target, crs$risk, paste(crs$input, collapse=" + "))
    crs$survival <-
      buildModelSurvival(formula,
                         dataset,
                         tv=theWidget("model_survival_textview"),
                         method=ifelse(theWidget("model_survival_para_radiobutton")$
                           getActive(), "para", "coxph"))
    if (not.null(crs$survival))
    {
      showModelSurvivalExists()
      theWidget("evaluate_survival_checkbutton")$setActive(TRUE)
    }
    else
      setStatusBar(sprintf(Rtxt("Building %s model ... failed."),
                           commonName(crv$SURVIVAL)))
  }
  
  if (build.all)
  {
    time.taken <- Sys.time()-start.time
    time.msg <- sprintf(Rtxt("Time taken: %0.2f %s"), time.taken,
                        attr(time.taken, "units"))
    setStatusBar(Rtxt("All models have been generated."), time.msg)
  }
}

#----------------------------------------------------------------------
#
# MODEL GLM
#

executeModelGLM <- function()
{
  # Initial setup. 
  
  TV <- "glm_textview"
  mtype <- "linear"

  # Obtain the family

  if (theWidget("glm_linear_radiobutton")$getActive())
    family <- "Linear"
  else if (theWidget("glm_gaussian_radiobutton")$getActive())
    family <- "Gaussian"
  else if (theWidget("model_linear_poisson_radiobutton")$getActive())
    family <- "Poisson"
  else if (theWidget("glm_logistic_radiobutton")$getActive())
    family <- "Logistic"
  else if (theWidget("model_linear_probit_radiobutton")$getActive())
    family <- "Probit"
  else if (theWidget("glm_multinomial_radiobutton")$getActive())
    family <- "Multinomial"
  
  # Build the formula for the model. 080719 If the user has requested
  # a numeric target and the target is actually a factor, then convert
  # to a numeric, else the algorithms complain.

  if (family %in% c("Linear", "Gaussian", "Poisson")
      && "factor" %in% class(crs$dataset[[crs$target]]))
    frml <- sprintf("as.numeric(%s) ~ .", crs$target)
  else
    frml <- paste(crs$target, "~ .")

  # List, as a string, the variables to be included. 
  
  included <- "c(crs$input, crs$target)" # 20110102 getIncludedVariables()
  
  # Some convenience booleans.

  sampling  <- not.null(crs$train)
  including <- not.null(included)
  subsetting <- sampling || including
  
  startLog(Rtxt("Regression model"))

  if (family == "Logistic" || family == "Probit")
  {
    # For a categoric variable we usually default to assuming
    # proprtions data, and so we perform logistic regression, which
    # uses a binomial distribution and a logit link function. However,
    # the user could eventually choose a different distriubtion/link
    # pair.
    #
    # If we have a binary response it may be that we might consider
    # using a loglog link rather than a logit link.

    model.cmd <- paste("crs$glm <- glm(", frml, ",\n    data=crs$dataset",
                       if (subsetting) "[",
                       if (sampling) "crs$train",
                       if (subsetting) ", ",
                       if (including) included,
                       if (subsetting) "]",
                       sprintf(',\n    family=binomial(link="%s")',
                               ifelse(family=="Probit", "probit", "logit")),
                       #", na.action=na.pass",
                       if (! is.null(crs$weights))
                              sprintf(",\n    weights=(%s)%s",
                                      crs$weights,
                                      ifelse(sampling, "[crs$train]", "")),
                       ")", sep="")

    # In addition to the default summary, add the chi-square test of
    # the difference between the null model and the current model as
    # presented in http://www.ats.ucla.edu/stat/R/dae/probit.htm.

    
    
    summary.cmd <- paste("print(summary(crs$glm))",
                         paste('cat(sprintf("Log likelihood: %.3f (%d df)\\n",',
                               '            logLik(crs$glm)[1],',
                               '            attr(logLik(crs$glm), "df")))', sep="\n"),
                         paste('cat(sprintf("Null/Residual deviance difference:',
                               '%.3f (%d df)\\n",'),
                         '            crs$glm$null.deviance-crs$glm$deviance,',
                         '            crs$glm$df.null-crs$glm$df.residual))',
                         'cat(sprintf("Chi-square p-value: %.8f\\n",',
                         '            dchisq(crs$glm$null.deviance-crs$glm$deviance,',
                         '                   crs$glm$df.null-crs$glm$df.residual)))',
                         'cat(sprintf("Pseudo R-Square (optimistic): %.8f\\n",',
                         '             cor(crs$glm$y, crs$glm$fitted.values)))',
                         "cat('\\n==== ANOVA ====\\n\\n')",
                         'print(anova(crs$glm, test="Chisq"))',
                         'cat("\\n")',
                         sep="\n")
  }
  else if (family == "Linear")
  {

    # For a numeric target we expect to produce the usual linear
    # model. We could use glm to generate the model using the gaussian
    # distribution and the identity link function. This will produce
    # the same model as lm. But lm is faster (glm is an iterative
    # algorithm) and it also produces the R squared stats, so we use
    # lm.
    
    model.cmd <- paste("crs$glm <- lm(", frml, ", data=crs$dataset",
                       if (subsetting) "[",
                       if (sampling) "crs$train",
                       if (subsetting) ",",
                       if (including) included,
                       if (subsetting) "]",
                       if (! is.null(crs$weights))
                       sprintf(",\n    weights=(%s)%s",
                               crs$weights,
                               ifelse(sampling, "[crs$train]", "")),
                      ")", sep="")

    summary.cmd <- paste("print(summary(crs$glm))",
                         "cat('==== ANOVA ====\n\n')",
                         "print(anova(crs$glm))",
                         'print("\n")',
                         sep="\n")
  }
  else if (family == "Gaussian")
  {
    # Whilst this is a less efficient equivalent of the Linear model
    # using lm, it is identified that some users perceive value in
    # having both lm and glm options for numeric regression. This uses
    # a gaussian distribution and an identity link function.

    model.cmd <- paste("crs$glm <- glm(", frml, ", data=crs$dataset",
                       if (subsetting) "[",
                       if (sampling) "crs$train",
                       if (subsetting) ",",
                       if (including) included,
                       if (subsetting) "]",
                       ", family=gaussian(identity)",
                       ")", sep="")

    summary.cmd <- paste("print(summary(crs$glm))",
                         "cat('==== ANOVA ====\n\n')",
                         "print(anova(crs$glm))",
                         'print("\n")',
                         sep="\n")
  }
  else if (family == "Poisson")
  {
    # 080912 Added

    model.cmd <- paste("crs$glm <- glm(", frml, ", data=crs$dataset",
                       if (subsetting) "[",
                       if (sampling) "crs$train",
                       if (subsetting) ",",
                       if (including) included,
                       if (subsetting) "]",
                       ", family=poisson(log)",
                       ")", sep="")

    summary.cmd <- paste("print(summary(crs$glm))",
                         "cat('==== ANOVA ====\n\n')",
                         "print(anova(crs$glm))",
                         'print("\n")', sep="\n")
  }
  else if (family == "Multinomial")
  {
    lib.cmd <-  "library(nnet, quietly=TRUE)"
    if (! packageIsAvailable("nnet", Rtxt("build a mulitnomial model"))) return(FALSE)
    appendLog(Rtxt("Build a multinomial model using the nnet package."), lib.cmd)
    eval(parse(text=lib.cmd))

    car.available <- TRUE
    lib.cmd <- "library(car, quietly=TRUE)"
    if (! packageIsAvailable("car", Rtxt("use Anova to evaluate a mulitnomial model")))
      car.available <- FALSE
    else
    {
      appendLog(Rtxt("Summarise multinomial model using Anova from the car package."), lib.cmd)
      eval(parse(text=lib.cmd))
    }
    
    model.cmd <- paste("crs$glm <- ",
                       "multinom",
                       "(", frml, ", data=crs$dataset",
                       if (subsetting) "[",
                       if (sampling) "crs$train",
                       if (subsetting) ",",
                       if (including) included,
                       if (subsetting) "]",
                       ", trace=FALSE, maxit=1000",
                       ")", sep="")

    summary.cmd <- paste("rattle.print.summary.multinom(summary(crs$glm,",
                         "                              Wald.ratios=TRUE))",
                         paste('cat(sprintf("Log likelihood: %.3f (%d df)\n",',
                               'logLik(crs$glm)[1], attr(logLik(crs$glm), "df")))'),
                         paste('if (is.null(crs$glm$na.action)) omitted <- TRUE',
                               'else omitted <- -crs$glm$na.action'),
                         paste('cat(sprintf("Pseudo R-Square: %.8f\n\n",',
                               'cor(apply(crs$glm$fitted.values, 1, ',
                               'function(x) which(x == max(x))),\n',
                               'as.integer(crs$dataset',
                               ifelse(sampling, '[crs$train,]', ''),
                               '[omitted,]$',
                               crs$target, '))))\n', sep=""),
                         "cat('==== ANOVA ====\n')",
                         "print(Anova(crs$glm))",
                         'print("\n")',
                         sep="\n")

  }
  
  # Build the model.

  appendLog(Rtxt("Build a Regression model."),
            model.cmd, sep="")
  start.time <- Sys.time()
  result <- try(eval(parse(text=model.cmd)), silent=TRUE)
  if (inherits(result, "try-error"))
  {
    if (any(grep("too many (.*) weights", result)))
    {
      find.num.weights <- regexpr('\\([0-9]*\\)', result)
      num.weights <- substr(result, find.num.weights,
                            find.num.weights + attr(find.num.weights, "match.length")-1)
      errorDialog(sprintf(Rtxt("The Multinomial model build has failed,",
                               "with too many weights (%d)",
                               "needing to be calculated.",
                               "Perhaps consider reducing the",
                               "number of categoric variables with",
                               "unique values (if you have",
                               "such variables in your input data)",
                               "or perhaps treating the target",
                               "variable as numeric and perform a",
                               "numeric linear regression."),
                          num.weights))
      setTextview(TV)
    }        
    else if (any(grep("contrasts can be applied only to factors with 2", result)))
    {
      factors <- crs$input[sapply(crs$input, function(x)
                                  is.factor(crs$dataset[[x]]))]
      single <- factors[sapply(factors, function(x)
                               length(levels(crs$dataset[[x]]))==1)]
      one <- length(single)==1
      errorDialog("It appears that", ifelse(one, "a", "some"),
                  "categoric input",
                  ifelse(one, "variable is", "variables are"), "constant.",
                  "The regression model algorithm can not handle such",
                  ifelse(one, "a variable.", "variables."),
                  "You may like to Ignore the",
                  ifelse(one, "variable", "variables"),
                  "through the Data tab:\n\n",
                  paste(single, collapse=", "))
      setTextview(TV)
    }
    else
      errorDialog(Rtxt("The regression model appears to have failed.",
                       "The error message was:"),
                  result, crv$support.msg)
    return(FALSE)
  }
  
  # Summarise the model.

  appendLog(sprintf(Rtxt("Generate a textual view of the %s model."),
                    commonName(crv$GLM)), summary.cmd)
  
  resetTextview(TV)
  setTextview(TV, sprintf(Rtxt("Summary of the %s %s model",
                               "(built using %s):\n"),
                          family, "Regression",
                          ifelse(family == "Linear", "lm",
                                 ifelse(family == "Multinomial", "multinom", "glm"))),
              ifelse(any(is.na(coef(crs$glm))),
                     paste("\n***Note*** Singularities were found in the modeling
and are indicated by an NA in the following table.
This is often the case when variables are linear
combinations of other variables, or the variable
has a constant value.  These variables will be ignored
when using the model to score new data and will not be
included as parameters in the exported scoring routine.\n",
                           sep="\n"), ""),
              collectOutput(summary.cmd))

# 090223 Got it working above - just wrap a print around it?
#
#  if (family == "Multinomial" && car.available)
#  {
#    # Couldn't get this working within the summary.cmd
#    appendTextview(TV, paste("\n\n",
#                             collectOutput("Anova(crs$glm)", use.print=TRUE)),
#                   tvsep=FALSE)
#  }
  
  
  if (sampling) crs$smodel <- union(crs$smodel, crv$GLM)

  # Enable the plot button

  if (family == "Multinomial")
    theWidget("model_linear_plot_button")$setSensitive(FALSE)
  else
    theWidget("model_linear_plot_button")$setSensitive(TRUE)
  
  # Finish up.
  
  time.taken <- Sys.time()-start.time

  reportTimeTaken(TV, time.taken,
                  msg=paste(sprintf(Rtxt("The %s model has been built."),
                    commonName(mtype)),
                    ifelse(any(is.na(coef(crs$glm))),
                           Rtxt("Singularities exist."), "")))
  return(TRUE)
}

rattle.print.summary.multinom <- function (x, digits = x$digits, ...) 
{
  # All I want is to add "n=XXX" here!!!
  if (!is.null(cl <- x$call)) {
        cat("Call:\n")
        dput(cl, control = NULL)
    }
    cat(sprintf("\nn=%d\n", nrow(x$weights)))
    cat("\nCoefficients:\n")
    if (x$is.binomial) {
        print(cbind(Values = x$coefficients, "Std. Err." = x$standard.errors, 
            "Value/SE" = x$Wald.ratios), digits = digits)
    }
    else {
        print(x$coefficients, digits = digits)
        cat("\nStd. Errors:\n")
        print(x$standard.errors, digits = digits)
        if (!is.null(x$Wald.ratios)) {
            cat("\nValue/SE (Wald statistics):\n")
            print(x$coefficients/x$standard.errors, digits = digits)
        }
    }
    cat("\nResidual Deviance:", format(x$deviance), "\n")
    cat("AIC:", format(x$AIC), "\n")
    if (!is.null(correl <- x$correlation)) {
        p <- dim(correl)[2]
        if (p > 1) {
            cat("\nCorrelation of Coefficients:\n")
            ll <- lower.tri(correl)
            correl[ll] <- format(round(correl[ll], digits))
            correl[!ll] <- ""
            print(correl[-1, -p], quote = FALSE, ...)
        }
    }
    invisible(x)
}

exportRegressionModel <- function()
{
  # Make sure we have a model first!

  if (noModelAvailable(crs$glm, crv$GLM)) return(FALSE)

  startLog(Rtxt("Export regression model."))

  save.name <- getExportSaveName(crv$GLM)
  if (is.null(save.name)) return(FALSE)
  ext <- tolower(get.extension(save.name))

  # Generate appropriate code.
  
  if (! is.null(crs$weights))
    wt <- gsub("^\\(|\\)$", "",
               gsub("crs\\$dataset\\$|\\[.*\\]", "",
                    capture.output(print(crs$weights))))
  
  pmml.cmd <- sprintf("pmml(crs$glm%s%s)",
                      ifelse(length(crs$transforms) > 0,
                             ", transforms=crs$transforms", ""),
                      ifelse(is.null(crs$weights), "",
                             paste(", weights=", wt, sep="")))

  if (ext == "xml")
  {
    appendLog(Rtxt("Export regression as PMML."),
              sprintf('saveXML(%s, "%s")', pmml.cmd, save.name))
    XML::saveXML(eval(parse(text=pmml.cmd)), save.name)
  }
  else if (ext == "c")
  {
    # 090103 gjw Move to a function: saveC(pmml.cmd, save.name, "regression")

    # 090223 Why is this tolower being used? Under GNU/Linux it is
    # blatantly wrong. Maybe only needed for MS/Widnows
    
    if (isWindows()) save.name <- tolower(save.name)
    
    model.name <- sub("\\.c", "", basename(save.name))
    
    export.cmd <- generateExportPMMLtoC(model.name, save.name, "glm_textview")
    
    appendLog(sprintf(Rtxt("Export %s as a C routine."), commonName(crv$GLM)),
              sprintf('pmml.cmd <- "%s"\n\n', pmml.cmd),
              export.cmd)

    eval(parse(text=export.cmd))
  }
  
  setStatusBar(sprintf(Rtxt("The %s file '%s' has been written."),
                       toupper(ext), save.name))
}

#------------------------------------------------------------------------
# MODEL SVM - SUPPORT VECTOR MACHINE
#
# Eventually the following will go into svm_gui.R, using the same
# conventions as in ada.R.

on_svm_kernel_combobox_changed <- function(action, window)
{
  krnl <- theWidget("svm_kernel_combobox")$getActiveText()
  krnl <- gsub(").*$", "", gsub("^.*\\(", "",  krnl))
  if (length(krnl)) setGuiDefaultsSVM(krnl)
}

setGuiDefaultsSVM <- function(kernel=NULL)
{
  if (is.null(kernel))
  {
    theWidget("svm_kernel_combobox")$setActive(0)
    theWidget("svm_classweights_entry")$setText("")
    theWidget("svm_poly_degree_spinbutton")$setValue(1)
  }
  else
  {
    if (kernel=="polydot")
    {
      theWidget("svm_poly_degree_label")$show()
      theWidget("svm_poly_degree_spinbutton")$show()
    }
    else
    {
      theWidget("svm_poly_degree_label")$hide()
      theWidget("svm_poly_degree_spinbutton")$hide()
    }
  }
}

## Eventually the following will go into svm.R, using the same
## conventions as in ada.R.

executeModelSVM <- function()
{
  # DESCRIPTION
  # Build a support vector machine predictor.
  #
  # RETURNS
  # Ignored.
  #
  # DETAILS There are two model builders for SVMs: The e1071 version
  # is older and is supported by tune, and the kernlab version is
  # much more extensive. I did move back to e1071 because I thought
  # issues around the handling of NAs in kernlab a problem, but
  # essentially I think it is an issue with svm using all variables,
  # so I had to clean up my handling of NAs.
  
  useKernlab <- theWidget("kernlab_radiobutton")$getActive()

  TV <- ifelse(useKernlab, "ksvm_textview", "esvm_textview")
  
  startLog(Rtxt("Support vector machine."))

  ## Library.

  if (useKernlab)
  {
    if (packageIsAvailable("kernlab", Rtxt("build an SVM model using ksvm")))
    {
      libCmd <- "library(kernlab, quietly=TRUE)"
      appendLog(packageProvides('kernlab', 'ksvm'), libCmd)
    }
    else
      return(FALSE)
  }
  else
  {
    if (packageIsAvailable("e1071", Rtxt("build an SVM model using svm")))
    {
      libCmd <- "library(e1071, quietly=TRUE)"
      appendLog(packageProvides('e1071', 'svm'), libCmd)
    }
    else
      return(FALSE)
   }
  eval(parse(text=libCmd))

  ## Formula. TODO For kernlab we assume we will always do
  ## classification rather than regression, at least for now.

  if (useKernlab)
    frml <- paste(ifelse(numericTarget(),
                         crs$target,
                         sprintf("as.factor(%s)", crs$target)),
                  "~ .")
  else
    frml <- paste(crs$target, "~ .")

  ## Interface options.

  krnl <- theWidget("svm_kernel_combobox")$getActiveText()
  krnl <- gsub(").*$", "", gsub("^.*\\(", "",  krnl))

  opts <- theWidget("model_svm_options_entry")$getText()
  
  if (krnl == "polydot")
    degree <- theWidget("svm_poly_degree_spinbutton")$getValue()
  
  cweights <- theWidget("svm_classweights_entry")$getText()
  
  ## Included variables.

  included <- "c(crs$input, crs$target)" # 20110102 getIncludedVariables()
  
  ## Convenience booleans.

  sampling   <- not.null(crs$train)
  including  <- not.null(included)
  subsetting <- sampling || including

  # Parameters.

  parms <- ""
  if (krnl != "")
    parms <- sprintf('%s,\n      kernel="%s"', parms, krnl)
  if (cweights != "")
    parms <- sprintf('%s,\n      class.weights=%s', parms, cweights)
  if (krnl == "polydot")
    parms <- sprintf('%s,\n      kpar=list("degree"=%s)', parms, degree)
  if (nchar(opts) > 0)
    parms <- sprintf('%s,\n      %s', parms, opts)

  # Dataset

  dataset <- paste("crs$dataset",
                   if (subsetting) "[",
                   if (sampling) "crs$train",
                   if (subsetting) ",",
                   if (including) included,
                   if (subsetting) "]",
                   sep="")

  # Replicate rows according to the integer weights variable.
  
  if(! is.null(crs$weights))
    dataset <- paste(dataset,
                     "[rep(row.names(",
                     dataset,
                     "),\n                                        ",
                     # Use eval since crs$weights could be a formula
                     'as.integer(eval(parse(text = "', crs$weights,
                     '"))[crs$sample])),]',
                     sep="")
  
  
  # Build the model.

  if (useKernlab)
    svmCmd <- paste("crs$ksvm <- ksvm(", frml, ",\n      data=", dataset, sep="")
  else
    svmCmd <- paste("crs$svm <- svm(", frml, ",\n      data=", dataset, sep="")
  svmCmd <- paste(svmCmd, parms, sep="")

  # We need to add an option to ensure the probabilities are also
  # recorded.

  if (useKernlab)
    svmCmd <- paste(svmCmd, ",\n      prob.model=TRUE", sep="")  # Probabilities
  else
    svmCmd <- paste(svmCmd, ",\n      probability=TRUE", sep="")  # Probabilities
  svmCmd <- paste(svmCmd, ")", sep="")

  # 100614 Ensure each run generates the same model.

  seed.cmd <- "set.seed(crv$seed)"
  svmCmd <- paste(seed.cmd, svmCmd, sep="\n")
  
  start.time <- Sys.time()
  appendLog(Rtxt("Build a Support Vector Machine model."), svmCmd)
  result <- try(eval(parse(text=svmCmd)), silent=TRUE)
  if (inherits(result, "try-error"))
  {
    if (any(grep("cannot allocate vector", result)))
    {
      errorDialog(Rtxt("The call to svm has failed, running out of memory.",
                       "The support vector machine algorithm is rather memory hungry.",
                       "A quick solution is to down sample the dataset, through the",
                       "Data tab. On 32 bit machines you may be limited to",
                       "less than 10,000 entities."))
      setTextview(TV)
    }
    else
      errorDialog(errorMessageFun("svm", result))
    return(FALSE)
  }

  # Display the resulting model.

  if (useKernlab)
    summaryCmd <- "crs$ksvm"
  else
    summaryCmd <- "crs$svm"
  appendLog(sprintf(Rtxt("Generate a textual view of the %s model."),
                    commonName(crv$KSVM)), summaryCmd)
  resetTextview(TV)
  setTextview(TV, sprintf(Rtxt("Summary of the %s model",
                               "(built using ksvm):"),
                          commonName("ksvm")),
              "\n\n",
              collectOutput(summaryCmd, TRUE), "\n")

  if (sampling)
    if (useKernlab)
      crs$smodel <- union(crs$smodel, crv$KSVM)
    else
      crs$smodel <- union(crs$smodel, crv$SVM)

  # Finish up.

  time.taken <- Sys.time()-start.time

  reportTimeTaken(TV, time.taken, ifelse(useKernlab, crv$KSVM, crv$SVM))

  return(TRUE)
}

exportSVMModel <- function()
{
  # Make sure we have a model first!
  
  if (is.null(crs$ksvm))
  {
    errorDialog(Rtxt("No SVM model is available. Be sure to build",
                     "the model before trying to export it! You will need",
                     "to press the Execute button (F2) in order to build the",
                     "model."))
    return()
  }

  # Require the pmml package
  
  lib.cmd <- "library(pmml, quietly=TRUE)"
  if (! packageIsAvailable("pmml", Rtxt("export SVM model"))) return(FALSE)
  appendLog(Rtxt("Load the PMML package to export a SVM model."), lib.cmd)
  # Load the package unless we already have a pmml defined (through source).
  if (! exists("pmml")) eval(parse(text=lib.cmd))
  
  # Obtain filename to write the PMML to.
  
  dialog <- RGtk2::gtkFileChooserDialog(Rtxt("Export PMML"), NULL, "save",
                                 "gtk-cancel", RGtk2::GtkResponseType["cancel"],
                                 "gtk-save", RGtk2::GtkResponseType["accept"])
  dialog$setDoOverwriteConfirmation(TRUE)

  if(not.null(crs$dataname))
    dialog$setCurrentName(paste(get.stem(crs$dataname), "_ksvm.xml", sep=""))

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

#  if (get.extension(save.name) == "") save.name <- sprintf("%s.xml", save.name)

  if (! is.null(crs$weights))
    wt <- gsub("^\\(|\\)$", "",
               gsub("crs\\$dataset\\$|\\[.*\\]", "",
                    capture.output(print(crs$weights))))
  
  pmml.cmd <- paste('pmml(crs$ksvm, dataset=crs$dataset',
                    if (! is.null(crs$weights))
                    paste(', weights=', wt, sep=""),
                    ')', sep="")
  appendLog(Rtxt("Export a SVM model as PMML."),
            sprintf('saveXML(%s, "%s")', pmml.cmd, save.name))
  XML::saveXML(eval(parse(text=pmml.cmd)), save.name)

  setStatusBar(sprintf(Rtxt("The PMML file '%s' has been written."), save.name))
}

##----------------------------------------------------------------------
##
## MARS
##
## y <- pchresp[, c(1)]
## x <- pchresp[, -c(1)]
##
## m1 <- mars(x, y)
##
## showcuts <- function(obj)
## {
##   tmp <- obj$cuts[obj$sel, ]
##   dimnames(tmp) <- list(NULL, dimnames(x)[[2]])
##   tmp
## }
##
## m2 <- mars(x, y, degree=2)

##----------------------------------------------------------------------
##
## SVM
##
## > >
## > > I am using the "svm" command in the e1071 package.
## > >
## > > Does it have an automatic way of setting the "cost" parameter?
## >
## > See ?best.svm in that package.
## >
## > > I changed a few values for the "cost" parameter but I hope there is a
## > > systematic way of obtaining the best "cost" value.
## > >
## > > I noticed that there is a "cross" (Cross validation)
## > > parameter in the "svm"
## > > function.
## > >
## > > But I did not see how it can be used to optimize the "cost" parameter.
## > >
## > > By the way, what does a 0 training error and a high testing
## > > error mean?
## > > Varying "cross=5", or "cross=10", etc. does not change the
## > > training error
## > > and testing error at all. How to improve?
## >
## > Overfitting, which varying different validation method will not solve.


## You might find http://www.csie.ntu.edu.tw/~cjlin/papers/guide/guide.pdf
## <http://www.csie.ntu.edu.tw/~cjlin/papers/guide/guide.pdf>  helpful.

## Parameter tuning is essential for avoiding overfitting.


########################################################################
# Export

exportModelTab <- function()
{
  # 090812 The Export button has been clicked whilst the Model tab is
  # active. Export the active Model as appropriate (either PMML or C
  # code, if C code export is available.)
  
  if (noDatasetLoaded()) return()

  # 090518 Test if each of the variables is exportable. If not (i.e.,
  # one or more are transforms that are not supported) then put up a
  # warning but continue. RStat may decide to not continue, but Rattle
  # should since perhaps the PMML is provided as a template to then
  # allow a user to edit.

  ## 131114 pmmlCanExport no longer available in pmml package.
  ##
  ## if (any(!sapply(crs$input, pmmlCanExport)))
  ## {
  ##   if (!questionDialog(Rtxt("In exporting the model the following variables appear",
  ##                            "to be transformations that are not currently supported",
  ##                            "for export. Be aware that the results will not",
  ##                            "perform the transformations."),
  ##                       "\n\n",
  ##                       paste(names(which(!sapply(crs$input, pmmlCanExport))),
  ##                             collapse=", "),
  ##                       "\n\n",
  ##                       Rtxt("Do you wish to continue?")))
  ##     return()
  ## }
  
  if (theWidget("rpart_radiobutton")$getActive())
  {
    exportRpartModel()
  }
  else if (theWidget("rf_radiobutton")$getActive())
  {
    exportRandomForestModel()
  }
  else if (theWidget("model_linear_radiobutton")$getActive())
  {
    exportRegressionModel()
  }
  else if (theWidget("svm_radiobutton")$getActive())
  {
    exportSVMModel()
  }
  else if (theWidget("nnet_radiobutton")$getActive())
  {
    exportNNetModel()
  }
  else if (theWidget("model_survival_radiobutton")$getActive())
  {
    exportSurvivalModel()
  }
  else
  {
    errorDialog(Rtxt("PMML export for this model is not yet implemented."))
    return()
  }
}
