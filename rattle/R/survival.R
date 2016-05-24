# Rattle Survival
#
# Time-stamp: <2015-05-17 08:58:46 gjw>
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
# GUI

setGuiDefaultsSurvival <- function()
{
  theWidget("model_survival_time_var_label")$setText(Rtxt("No time variable selected"))
  theWidget("model_survival_status_var_label")$setText(Rtxt("No status variable selected"))
  theWidget("model_survival_coxph_radiobutton")$setActive(TRUE)
  theWidget("model_survival_plots_label")$setSensitive(FALSE)
  theWidget("model_survival_plot_survival_button")$setSensitive(FALSE)
  theWidget("model_survival_plot_residual_button")$setSensitive(FALSE)
}

on_model_survival_coxph_radiobutton_toggled <- function(button)
{
  if (button$getActive())
    theWidget("model_survival_function_label")$setText("coxph")
  else
    # 091114 This is the only other alternative for now, so only use
    # the one callback unless we add alternative buidlers.
    theWidget("model_survival_function_label")$setText("survreg")

  activate.buttons <- button$getActive() && ! is.null(crs$survival)

  theWidget("model_survival_plots_label")$setSensitive(activate.buttons)
  theWidget("model_survival_plot_survival_button")$setSensitive(activate.buttons)
  theWidget("model_survival_plot_residual_button")$setSensitive(activate.buttons)
}

on_model_survival_plot_survival_button_clicked <- function(button)
{
  plotSurvivalModel()
}


on_model_survival_plot_residual_button_clicked <- function(button)
{
  plotResidualModels()
}


########################################################################
# Model Tab

buildModelSurvival <- function(formula, dataset, tv=NULL, method=c("para", "coxph"))
{
  # If tv is not NULL, then we will be updating the textview object as
  # we proceed, as well as sending information to the log. The aim is
  # for this function to run totally independent of the GUI, but to
  # also support it. A developer can use this function, supply their
  # own textview object and their own implementations of resetTextview
  # and appendTextview for the modelling output, and startLog and
  # appendLog for a log of the commands, and setStatusBar for a
  # summary of what has been done.

  gui <- not.null(tv)
  if (gui) startLog(Rtxt("Survival Model"))

  sampling <- not.null(crs$sample)

  # Load the required package into the library.

  lib.cmd <-  "library(survival, quietly=TRUE)"
  if (! packageIsAvailable("survival", Rtxt("build a Survival model"))) return(NULL)
  if (gui) appendLog(Rtxt("Require the survival package."), lib.cmd)
  eval(parse(text=lib.cmd))

  # Build a model. 

  method <- ifelse(method=="para", "survreg", "coxph")
  model.cmd <- paste(method, "(", formula,
                     ",\n      data=", dataset,
                     if (! is.null(crs$weights))
                        sprintf(",\n    weights=(%s)%s",
                                crs$weights,
                                ifelse(sampling, "[crs$train]", "")),
                     ")", sep="")

  if (gui) appendLog(Rtxt("Build the Survival model."),
                     sprintf('crs$survival <- %s', model.cmd))

  # Note that this crs$survival is not the global crs$survival! We use
  # it here to be consistent in terms of the commands that are
  # reported to the log, but we return this value and in the outer
  # call we globally assign to crs$survival, at least in the context
  # of the Rattle GUI.
  
  start.time <- Sys.time()
  crs$survival <- try(eval(parse(text=model.cmd)), silent=TRUE)
  time.taken <- Sys.time()-start.time

  if (inherits(crs$survival, "try-error"))
  {
    msg <- errorMessageFun(method, crs$survival)
    
    if (any(grep(Rtxt("Invalid survival times for this distribution"), crs$survival)))
    {
      errorDialog(Rtxt("The building of the survival model failed.",
                       "The error indicates an invalid Time variable.",
                       "This can be the case when using survreg and there is",
                       "a zero time value (as might result from an imputation).",
                       "Please review the source data and ensure the Time values",
                       "are correct."))
      setTextview(tv)
    }
    else if (any(grep(Rtxt("NA/NaN/Inf in foreign function call"), crs$survival)))
    {
      errorDialog(Rtxt("Your data contains variables with too many categoric values.",
                       "Please reduce the number of categoric values or remove any",
                       "identifier variables from the input variables in order to",
                       "generate a survival model.",
                       "\n\nThe actual error message was:"),
                  "\n\n", paste(crs$survival, "\n"))
      setTextview(tv)
    }
    else
    {
      if (gui)
      {
        errorDialog(msg)
        return(NULL)
      }
    }
    return(FALSE)
  }

  # Print the results of the modelling.

  if (gui)
  {
    print.cmd <- "summary(crs$survival)"
    appendLog(Rtxt("Print the results of the modelling."), print.cmd)
    resetTextview(tv, tvsep=FALSE,
                  sprintf(Rtxt("Summary of the Survival model (built using %s):"),
                          method),
                  "\n\n",
                  collectOutput(print.cmd))
    if (method=="coxph")
    {
      print.cmd <- paste(print.cmd, "cox.zph(crs$survival)", sep="; ")
      appendTextview(tv, tvsep=FALSE, "\n\n",
                     Rtxt("Test the proportional hazards ",
                          "assumption for a Cox regression model:"),
                     "\n\n",
                     collectOutput(print.cmd))
    }
  }

  # Finish up.
  
  if (gui)
  {
    time.msg <- sprintf("\nTime taken: %0.2f %s", time.taken,
                        attr(time.taken, "units"))
    appendTextview(tv, "\n", time.msg)
    appendLog(time.msg)
    setStatusBar(Rtxt("A survival model has been generated."), time.msg)
  }
  return(crs$survival)
}

showModelSurvivalExists <- function(state=!is.null(crs$survival))
{
  # If a survival model exists then make sensitive the relevant
  # buttons that require the model to exist. For the Survival model
  # this will be the plot functions.

  if (state && class(crs$survival) == "coxph")
  {
    theWidget("model_survival_plots_label")$setSensitive(TRUE)
    theWidget("model_survival_plot_survival_button")$setSensitive(TRUE)
    theWidget("model_survival_plot_residual_button")$setSensitive(TRUE)
  }

  theWidget("score_class_radiobutton")$
  setActive(class(crs$survival) == "survreg")
  theWidget("score_probability_radiobutton")$
  setSensitive(class(crs$survival) == "coxph")
}

plotSurvivalModel <- function()
{
  startLog(Rtxt("Survival chart."))

  plot.cmd <- paste('plot(survfit(crs$survival), xlab=crs$target,',
                    'ylab="Survival Probability", col=3)\n',
                    genPlotTitleCmd('Survival Chart', crs$target, 'to',
                                    crs$risk), sep="")
  appendLog(Rtxt("Plot the survival chart for",
                 "the most recent survival model."), plot.cmd)
  newPlot()
  eval(parse(text=plot.cmd))
}

plotResidualModels <- function()
{
  startLog(Rtxt("Survival model residuals plot."))

  # 100417 Use the max number per page count from the plots page.
  
  pmax <- theWidget("plots_per_page_spinbutton")$getValue()

  plot.cmd <- paste('temp <- cox.zph(crs$survival)',
                    sprintf("pmax <- %d", pmax),
                    "pcnt <- 0",
                    'nr <- nrow(temp$var)',
                    "if (nr < pmax) pmax <- nr",
                    'for (vnum in 1:nr)',
                    '{',
                    '  if (pcnt %% pmax == 0) newPlot(pmax)',
                    '  pcnt <- pcnt + 1',
                    '  plot(temp, var=vnum)',
                    '  abline(0, 0, lty=3)',
                    '  # A linear fit.',
                    '  abline(lm(temp$y[,vnum] ~ temp$x)$coefficients, lty=4, col=3)',
                    '}',
                    sep="\n")
  appendLog(Rtxt("Plot the scaled Schoenfeld residuals of proportional hazards."),
            plot.cmd)
  eval(parse(text=plot.cmd))
}

########################################################################
# Export

exportSurvivalModel <- function()
{
  # Make sure we have a model first!

  if (noModelAvailable(crs$survival, crv$SURVIVAL)) return(FALSE)

  startLog(Rtxt("Export survival model."))

  save.name <- getExportSaveName(crv$SURVIVAL)
  if (is.null(save.name)) return(FALSE)
  ext <- tolower(get.extension(save.name))

  # Generate appropriate code.
  
  pmml.cmd <- sprintf("pmml(crs$survival%s)",
                      ifelse(length(crs$transforms) > 0,
                             ", transforms=crs$transforms", ""))

  if (ext == "xml")
  {
    appendLog(Rtxt("Export survival regression as PMML."),
              sprintf('saveXML(%s, "%s")', pmml.cmd, save.name))
    XML::saveXML(eval(parse(text=pmml.cmd)), save.name)
  }
  
  setStatusBar("The", toupper(ext), "file", save.name, "has been written.")
}

########################################################################
# Evaluate

genPredictSurvival <- function(dataset)
{
  # Generate a command to obtain the prediction results when applying
  # the model to new data. For the coxph model we return the median of
  # the expected survival time. 091115 Note that the code to extract
  # this is somewhat convoluted until I generate the actual formula
  # that is used to calculate the median. This will also be needed to
  # convert to C code.
  
  is.coxph <- class(crs$survival) == "coxph"

  if (is.coxph)
    cmd <- sprintf(paste("crs$pr <- survival:::survmean(survfit(crs$survival,",
                         '%s), scale=1, rmean="none")[[1]][,5]'), dataset)
  else
    cmd <- sprintf("crs$pr <- predict(crs$survival, %s)", dataset)

  return(cmd)
}

genResponseSurvival <- function(dataset)
{
  # Generate a command to obtain the response when applying the model
  # to new data.
  
  return(genPredictSurvival(dataset))
}

genProbabilitySurvival <- function(dataset)
{
  # Generate a command to obtain the probability when applying the
  # model to new data. For the coxph model we predict the risk. This
  # is not a probability, but is a number relative to 1, so that
  # greater than 1 has a higher risk that the average of the event
  # occuring, and below 1 has a lower risk of the event occuring than
  # the average.
  
  is.coxph <- class(crs$survival) == "coxph"

  return(sprintf("crs$pr <- predict(crs$survival, %s%s)", dataset,
                 ifelse(is.coxph, ', type="risk"', "")))
}

