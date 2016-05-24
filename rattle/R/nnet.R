# Gnome R Data Miner: GNOME interface to R for Data Mining
#
# Time-stamp: <2015-07-13 20:25:20 gjw>
#
# NNET OPTION 061230
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

#######################################################################
#
# COMMENTS
#
# 080514 For classification tasks we use multinom which calls nnet. As
# in the documentation in R, multinom fits multinomial log-linear
# models via neural networks. If you use glm instead for the case of
# binary classificatio we get only a very slightly different
# model. See Wikipedia for details of GLM.
#
# TODO
#
# Figure out the actual framework for two class prediction
#
# use multinom from nnet package - but that is just same as glm.
#
##        mymn <- multinom(TARGET_Adjusted ~ .,
##                         data=crs$dataset[crs$sample,c(2:10,13)])
##        crs$pr <- predict(mymn, crs$dataset[-crs$sample, c(2:10,13)]) # (0,1)
##        table(crs$pr, crs$dataset[-crs$sample, c(2:10,13)]$TARGET_Adjusted,
##              dnn=c("Predicted", "Actual"))
##        crs$pr <- predict(mymn, crs$dataset[-crs$sample, c(2:10,13)],
##                          type="prob")
##        crs$eval <- evaluateRisk(crs$pr,
##                                 crs$dataset[-crs$sample,
##                                             c(2:10,13)]$TARGET_Adjusted,
##                                 crs$dataset[-crs$sample,
##                                             c(2:10,13,12)]$RISK_Adjustment)
##        plotRisk(crs$eval$Caseload, crs$eval$Precision,
##                 crs$eval$Recall, crs$eval$Risk)
#

########################################################################
#
# CALLBACKS
#

########################################################################
#
# NNET
#

executeModelNNet <- function()
{
  # 090820 nnet is only supported in Rattle at the moment for numeric
  # and binomial targets, through the Rattle interface. The
  # multinomial target is handled through multinom under Linear.
  
  # Initial setup. 
  
  TV <- "nnet_textview"

  # Obtain user interface model options.

  size <- theWidget("nnet_hidden_nodes_spinbutton")$getValue()
  
  # Load the package into the library

  startLog("Neural Network")
  lib.cmd <-  "library(nnet, quietly=TRUE)"
  if (! packageIsAvailable("nnet", Rtxt("build a neural network"))) return(FALSE)
  appendLog(Rtxt("Build a neural network model using the nnet package."), lib.cmd)
  eval(parse(text=lib.cmd))

  # Build the formula for the model.

#  if (binomialTarget() && ! is.numeric(crs$dataset[[crs$target]]))
    # 081010 Needs to be numeric, but I also subtract 1 so we get a
    # 0/1 target? 090820 But should it? It should use whatever nnet
    # does when it has a categoric target. Then we also need to remove
    # linout? Though again, check why linout is needed at all.
#    frml <- sprintf("as.numeric(%s)-1 ~ .", crs$target)
#  else

  # 091023 When the target is binary, the prediction uses type=class,
  # but that won't work if the target is binary and numeric. So we
  # test that case here and convert it to a factor. Note that we only
  # support numeric targets for nnet.

  if (binomialTarget())
    frml <- sprintf("as.factor(%s) ~ .", crs$target)
  else
    frml <- sprintf("%s ~ .", crs$target)
  
  # Variables to be included --- a string of indicies.
  
  # 20110102 included <- getIncludedVariables()
  included <- "c(crs$input, crs$target)"

  # Some convenience booleans

  sampling <- not.null(crs$sample)
  including <- not.null(included)
  subsetting <- sampling || including

  # Time the model building.
  
  start.time <- Sys.time()
  
  # Build a model. 091114 Note that we use a seed so that we get the
  # same model each time. Otherwise it can be disconcerting to the
  # user to see the model changing each time they click Execute, or
  # each time they come into the application. 101204 We should provide
  # the seed in the GUI and also maxiter in the GUI, but for now use
  # 199 since that seems to get at least a model that predicts some
  # "yes" for the weather dataset! Not a good foundation for setting
  # the seed, but fix it later.

  model.cmd <- paste(sprintf("set.seed(%d)\n", 199), # crv$seed),
                     "crs$nnet <- ",
                     ifelse(numericTarget() || binomialTarget(),
                            "nnet", "multinom"),
                     "(", frml, ",\n    data=crs$dataset",
                     if (subsetting) "[",
                     if (sampling) "crs$sample",
                     if (subsetting) ",",
                     if (including) included,
                     if (subsetting) "]",
                     if (! is.null(crs$weights))
                        sprintf(",\n    weights=(%s)%s",
                                crs$weights,
                                ifelse(sampling, "[crs$train]", "")),
                     # TODO 080427 How to choose a good value for size?
                     # TODO 090808 Why linout for a binomial target?
#                     if (numericTarget() || binomialTarget())
#                     sprintf(", size=%d, linout=TRUE, skip=TRUE", size),
                     sprintf(",\n    size=%d", size),
                     if (numericTarget()) ", linout=TRUE",
                     ", skip=TRUE",
                     ", MaxNWts=10000",
                     # 101204 When maxit=1000 the weather model only
                     # predicts 0. Set back to 100 then it "works."
                     # 100 is the default. Maybe need to allow as an
                     # option.
                     ", trace=FALSE, maxit=100",
                     ")", sep="")

  appendLog("Build the NNet model.", model.cmd)
  result <- try(eval(parse(text=model.cmd)), silent=TRUE)
  time.taken <- Sys.time() - start.time
  if (inherits(result, "try-error"))
  {
    if (any(grep("too many", result)))
    {
      errorDialog(sprintf(Rtxt("The call to 'nnet' has failed.",
                               "This appears to be due to there being too many",
                               "links and nodes in the neural network architecture.",
                               "You may want to reduce the number of nodes in the",
                               "hidden layer, or reduce the number of input variables,",
                               "particularly any categoric variables.",
                               "The actual error was:\n\n%s"),
                          result))
      setTextview(TV)
    }
else
errorReport(model.cmd, result)
    return(FALSE)
  }
  
  # Print the results of the modelling.

  if (numericTarget() || binomialTarget())
    print.cmd <- paste('cat(sprintf("A %s network with %d weights.\\n",',
                       '    paste(crs$nnet$n, collapse="-"),',
                       '    length(crs$nnet$wts)))',
                       'cat(sprintf("Inputs: %s.\\n",',
                       '    paste(crs$nnet$coefnames, collapse=", ")))',
                       'cat(sprintf("Output: %s.\\n",',
                       '    names(attr(crs$nnet$terms, "dataClasses"))[1]))',
                       'cat(sprintf("Sum of Squares Residuals: %.4f.\\n",',
                       '    sum(residuals(crs$nnet) ^ 2)))',
                       'cat("\\n")',
                       "print(summary(crs$nnet))", "cat('\\n')", sep="\n")
  else
    print.cmd <- "print(crs$nnet)"

  appendLog(Rtxt("Print the results of the modelling."), print.cmd)
  resetTextview(TV)
  setTextview(TV,
              sprintf(Rtxt("Summary of the Neural Net model (built using %s):"),
                      ifelse(numericTarget() || binomialTarget(),
                             Rtxt("nnet"), Rtxt("multinom"))),
              "\n\n",
              collectOutput(print.cmd))

  # Now that we have a model, make sure appropriate actions are sensitive.
  
  # Finish up.

  reportTimeTaken(TV, time.taken, commonName(crv$NNET))

  return(TRUE)
}

exportNNetModel <- function()
{
  # Make sure we have a model first! 090812 DRY move all
  # export<model>Tab functions to use this test instead of their
  # individual tests. 090812 DRY Unify much more the export<model>Tab
  # functions - they are all mostly the same, so don't repeat
  # yourself.

  if (noModelAvailable(crs$nnet, crv$NNET)) return(FALSE)

  startLog(paste(Rtxt("Export"), commonName(crv$RPART)))

  save.name <- getExportSaveName(crv$NNET)
  if (is.null(save.name)) return(FALSE)
  ext <- tolower(get.extension(save.name))
  
  # 090812 DRY - make regression and rpart all the same mold.

  ## REMOVE 090812 - Not consistent with rpart and regression so do we
  ## need this?
  ##
  ## # Require the pmml package
  ##
  ## lib.cmd <- "library(pmml, quietly=TRUE)"
  ## if (! packageIsAvailable("pmml", "export neural net")) return(FALSE)
  ## appendLog("Load the PMML package to export a neural net.", lib.cmd)
  ## # Load the package unless we already have a pmml defined (through source).
  ## if (! exists("pmml")) eval(parse(text=lib.cmd))

  pmml.cmd <- sprintf("pmml(crs$nnet%s)",
                      ifelse(length(crs$transforms),
                             ", transforms=crs$transforms", ""))
  if (ext == "xml")
  {
    appendLog(Rtxt("Export the Neural Net model as PMML."),
              sprintf('saveXML(%s, "%s")', pmml.cmd, save.name))
    XML::saveXML(eval(parse(text=pmml.cmd)), save.name)
  }
  else if (ext == "c")
  {
    if (isWindows()) save.name <- tolower(save.name)

    model.name <- sub("\\.c", "", basename(save.name))

    export.cmd <- generateExportPMMLtoC(model.name, save.name, "nnet_textview")
    
    appendLog(sprintf(Rtxt("Export %s as a C routine."), commonName(crv$NNET)),
              sprintf('pmml.cmd <- "%s"\n\n', pmml.cmd),
              export.cmd)

    eval(parse(text=export.cmd))

    
    ## appendLog(Rtxt("Export the Neural Net model as C code for WebFocus."),
    ##           sprintf('cat(pmmltoc(toString(%s), "%s", %s, %s, %s), file="%s")',
    ##                   pmml.cmd, model.name, 
    ##                   attr(save.name, "includePMML"),
    ##                   attr(save.name, "includeMetaData"),
    ##                   attr(save.name, "exportClass"),
    ##                   save.name))
    ## cat(pmmltoc(toString(eval(parse(text=pmml.cmd))), model.name,
    ##             attr(save.name, "includePMML"),
    ##             ifelse(attr(save.name, "includeMetaData"),
    ##                    getTextviewContent("nnet_textview"),
    ##                    "\"Not Included\""),
    ##             attr(save.name, "exportClass")), file=save.name)
  }
  setStatusBar(sprintf(Rtxt("The model has been exported to '%s'."), save.name))
}

# 140206 This is modelled on Riper's print.summary.nnet and is
# backward compatible and overrides the original when loaded which I
# think is safe since it is a print method, intended for display
# rather than down stream processing.

print.summary <- function(...) UseMethod("print.summary")

print.summary.nnet <- function(x, ...)
{
  eval(parse(text="library(nnet)")) # No longer for Log - remove it
  cat(Rtxt("Neural Network build options:"))
  tconn <- diff(x$nconn)
  ssep <- ""
  if (tconn[length(tconn)] > x$n[2L]+1L)
    {cat(ssep, "skip-layer connections"); ssep=";"}
  if (x$nunits > x$nsunits && !x$softmax)
    {cat(ssep, "linear output units"); ssep=";"}
  if (x$entropy)
    {cat(ssep, "entropy fitting"); ssep=";"}
  if (x$softmax)
    {cat(ssep, "softmax modelling"); ssep=";"}
  if (x$decay[1L] > 0)
    {cat(ssep, "decay=", x$decay[1L], sep=""); ssep=";"}
  cat(".\n\nIn the following table:\n",
      "  b  represents the bias associated with a node\n",
      "  h1 represents hidden layer node 1\n",
      "  i1 represents input node 1 (i.e., input variable 1)\n",
      "  o  represents the output node\n")
  wts <- format(round(coef(x),2))
  lapply(split(wts, rep(1L:x$nunits, tconn)),
         function(x)
         {
           cat(sprintf("\nWeights for node %s:\n",
                       sub("^.*->", "", names(x)[1])))
           print(x, quote=FALSE)
         })
  invisible(x)
}

