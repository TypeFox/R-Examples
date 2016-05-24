# Rattle: A GUI for Data Mining in R
#
# RPART TAB
#
# Time-stamp: <2015-07-26 11:55:23 gjw>
#
# Copyright (c) 2009-2014 Togaware Pty Ltd except as noted:
#
# drawTreeNodes() is a modification of draw.tree() from the maptree
# package written by Denis White. No copyright is claimed by Denis
# White but assumed by default to be Copyright (c) Denis White, and it
# is released on CRAN under an Unlimited license for any purpose.
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
# CALLBACKS
#

on_priors_entry_changed <- function(action, window)
{
  setStatusBar()
}

on_rpart_loss_comboboxentry_set_focus_child <- function(action, window)
{
  # When the rpart loss combobox is entered, retrieve the list of
  # matrix objects in R and populate the choices. This is the simplest
  # approach to handling the loss matrix at present and may be
  # sufficient.

  # Generate a list of suitable (matrix) R objects

  ml <- unlist(sapply(ls(sys.frame(0)),
                      function(x)
                      {
                        cmd <- sprintf("is.matrix(%s)",x)
                        var <- try(ifelse(eval(parse(text=cmd), sys.frame(0)),
                                          x, NULL), silent=TRUE)
                        if (inherits(var, "try-error"))
                          var <- NULL
                        return(var)
                      }))
  if (not.null(ml))
  {
    action$getModel()$clear()
    lapply(ml, action$appendText)
  }
}

on_rpart_plot_button_clicked <- function(button)
{

  # Make sure there is an rpart object first.

  if (is.null(crs$rpart))
  {
    errorDialog("E122: This is an unexpected error.", crv$support.msg)
    return()
  }

  # If there is only a root node then there is nothing to plot.

  if (theWidget("model_tree_rpart_radiobutton")$getActive() &&
      nrow(crs$rpart$frame) == 1)
  {
    errorDialog(Rtxt("The tree consists just of a root node.",
                     "There is nothing to plot."))
    return()
  }

  if (theWidget("model_tree_rpart_radiobutton")$getActive())
    if (theWidget("use_ggplot2")$getActive() # Not really ggplot2 but convenient.
        && packageIsAvailable("rpart.plot", "plot nice looking decision trees")
        && packageIsAvailable("RColorBrewer", "choose colours for tree plot"))
    {
        
      plot.cmd <- sprintf('fancyRpartPlot(crs$rpart, main="%s")',
                            genPlotTitleCmd(commonName(crv$RPART), crs$dataname,
                                            "$", crs$target, vector=TRUE)[1])
      
      log.txt <- "rpart.plot package"
    }
    else
    {
      plot.cmd <- paste("drawTreeNodes(crs$rpart)\n",
                        genPlotTitleCmd(commonName(crv$RPART),
                                        crs$dataname, "$", crs$target),
                        sep="")

      log.txt <- "maptools support functions"
    }
  else # ctree
  {
    plot.cmd <- "plot(crs$rpart)"
    log.txt <- "party package"
  }

  # Log the R command and execute.

  startLog(sprintf(Rtxt("Plot the resulting %s."), commonName(crv$RPART)))
           
  ##   plotcp.cmd <- paste("\n\n## Plot the cross validation results.\n\n",
  ##                           "plotcp(crs$rpart)\n",
  ##                           genPlotTitleCmd("Cross Validated Error",
  ##                                              crs$dataname, "$", crs$target))

  appendLog(sprintf(Rtxt("We use the %s."), log.txt), plot.cmd)

  # 150726 Do I need newPlot() here as when it is used the text in the
  # plot is obliterated wit a coloured box??? Perhaps I need to
  # eliminate all usage of newPlot() eventually.

  # newPlot()
  eval(parse(text=plot.cmd))

  ## newPlot()
  ## appendLog(plotcp.command)
  ## eval(parse(text=plotcp.command))

  setStatusBar(Rtxt("The decision tree has been plotted."))
}

on_rpart_rules_button_clicked <- function(button)
{
  ## Initial setup

  TV <- "rpart_textview"

  ## Make sure there is an rpart object first.

  if (is.null(crs$rpart))
  {
    errorDialog("E130: There is no rpart model yet.", crv$support.msg)
    return()
  }

  rules.cmd <- "asRules(crs$rpart)"
  appendLog(Rtxt("List the rules from the tree using a Rattle support function."),
          rules.cmd)
  addTextview(TV, paste(Rtxt("Tree as rules:"), "\n"), collectOutput(rules.cmd, TRUE),
              textviewSeparator())

  setStatusBar(paste(Rtxt("The corresponding rules have been listed."),
                     Rtxt("You may need to scroll the textview to view them.")))
}

on_help_rpart_activate <- function(action, window)
{
  if (showHelpPlus(Rtxt("A decision tree is the prototypical data mining tool,",
                        "used widely for its ease of interpretation. It consists of a root node",
                        "split by a single variable into two partitions. In turn, these two new",
                        "partitions become new nodes that may then each be further split on",
                        "a single (and usually",
                        "different) variable. This divide and conquering continues until no",
                        "further splitting would improve the performance of the model.",
                        "<<>>",
                        "While a choice of measures are available to select a variable to split",
                        "the dataset on, the Gini measure is used, and generally is no",
                        "different to the information measure for binary classification. To",
                        "explore the alternatives, copy the relevant code from the Log and",
                        "paste it into the R Console and change any of the options.",
                        "<<>>",
                        "Common options that a user may change from their default values are",
                        "available.",
                        "<<>>",
                        "Priors: used to boost a particularly important class, by giving it a",
                        "higher prior probability. Expects a list of numbers that sum up to 1,",
                        "and of the same length as the number of classes in the training dataset:",
                        "e.g.,0.5,0.5.",
                        "<<>>",
                        "Loss Matrix: used to weight the outcome classes differently:",
                        "e.g., 0,10,1,0.",
                        "<<>>",
                        "Other",
                        "options exist, but are not usually required. For example, 10-fold",
                        "cross validation, used in deciding how to prune to the best decision",
                        "tree, is generally regarded as the right number. Transferring the",
                        "commands from the Log tab into the R Console does give you full access",
                        "to all options.",
                        "<<>>",
                        "Decision trees work with both numeric and categoric data.",
                        "<<>>",
                        "The rpart package is used to build the decision tree.")))
  {
    popupTextviewHelpWindow("rpart", "rpart")
  }
}

on_rpart_build_radiobutton_toggled <- function(button)
{
  theWidget("rpart_tune_entry")$setSensitive(!button$getActive())
}

on_rpart_tune_radiobutton_toggled <- function(button)
{
  theWidget("rpart_tune_entry")$setSensitive(button$getActive())
}

on_rpart_best_radiobutton_toggled <- function(button)
{
  theWidget("rpart_tune_entry")$setSensitive(button$getActive())
}

########################################################################
#
# MODEL RPART
#

executeModelRPart <- function(action="build")
{
  # Initial setup

  TV <- "rpart_textview"

  num.classes <- length(levels(as.factor(crs$dataset[[crs$target]])))
  control <- NULL
  parms <- NULL

  # 100222 Use information as the default split method, as per a
  # machine learning view of the approach.

  parms <- ',\n    parms=list(split="information")'

  # Obtain the value of the tuning controls

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
                            priors,crs$target, num.classes))
        return(FALSE)
      }
    if (sum(pr) != 1)
      {
        errorDialog(sprintf(Rtxt("The supplied priors (%s)",
                                 "add up to %0.2f whereas",
                                 "they need to add up 1.00"),
                            priors, sum(pr)))
        return(FALSE)
      }
    if (is.null(parms))
      parms <- sprintf(",\n      parms=list(prior=c(%s))", priors)
    else
      parms <- gsub(")$", sprintf(",\n      prior=c(%s))", priors), parms)
  }

  # Retrieve the Min Split and check if it is different from the
  # default, and if so then use it.

  minsplit <- theWidget("rpart_minsplit_spinbutton")$getValue()
  if (minsplit != crv$rpart.minsplit.default)
  {
    if (is.null(control))
      control <- sprintf(",\n      control=rpart.control(minsplit=%d)", minsplit)
    else
      control <- gsub(")$", sprintf(",\n      minsplit=%d)", minsplit), control)
  }

  # Retrieve the Min Bucket and check if it is different from the
  # default, and if so then use it.

  minbucket <- theWidget("rpart_minbucket_spinbutton")$getValue()
  if (minbucket != crv$rpart.minbucket.default)
  {
    if (is.null(control))
      control <- sprintf(",\n      control=rpart.control(minbucket=%d)", minbucket)
    else
      control <- gsub(")$", sprintf(",\n           minbucket=%d)", minbucket), control)
  }

  # Retrieve the Max Depth and check if it is different from the
  # default, and if so then use it.

  maxdepth <- theWidget("rpart_maxdepth_spinbutton")$getValue()
  if (maxdepth != crv$rpart.maxdepth.default)
  {
    if (is.null(control))
      control <- sprintf(",\n      control=rpart.control(maxdepth=%d)", maxdepth)
    else
      control <- gsub(")$", sprintf(",\n           maxdepth=%d)", maxdepth), control)
  }

  # Retrieve the Complexity and check if it is different from the
  # default, and if so then use it.

  cp <- theWidget("model_tree_cp_spinbutton")$getValue()

  if (abs(cp-crv$rpart.cp.default) > 0.00001) ## Diff when same is 2.2352e-10!!!
  {
    if (is.null(control))
      control <- sprintf(",\n      control=rpart.control(cp=%f)", cp)
    else
      control <- gsub(")$", sprintf(",\n           cp=%f)", cp), control)
  }

  # Retrieve the Include Missing checkbutton status and if not set
  # then change default beahviour in usesurrogate.

  usesurrogate <- theWidget("model_tree_include_missing_checkbutton")$
                  getActive()
  if (! usesurrogate)
  {
    if (is.null(control))
      control <- paste(",\n    control=rpart.control(usesurrogate=0,",
                       "\n        maxsurrogate=0)")
    else
      control <- gsub(")$", paste(",\n        usesurrogate=0,",
                                  "\n        maxsurrogate=0)"), control)
  }

  # Retrieve the Cross Validation value and if different from
  # default, use it. No longer. Common wisdom is that 10 is right, so
  # in Rattle just go with that.

  # xval <- theWidget("rpart_xval_spinbutton")$getValue()
  # if (xval != crv$rpart.xval.default)
  # {
  #  if (is.null(control))
  #    control <- sprintf(", control=rpart.control(xval=%d)", xval)
  #  else
  #    control <- gsub(")$", sprintf(", xval=%d)", xval), control)
  # }

  # Retrieve the loss matrix and ensure it matches the shape of the
  # data.

  loss <- theWidget("model_tree_loss_entry")$getText()
  if (nchar(loss) > 0)
  {
    lo <- as.numeric(unlist(strsplit(loss, ",")))
    if (length(lo) != num.classes * num.classes)
    {
      errorDialog(sprintf(Rtxt("The supplied loss matrix (%s)",
                                "needs to have %d values.",
                                "Please enter that many values, comma separated."),
                          loss, num.classes*num.classes))
      return(FALSE)
    }

    # TODO: Perform other checks on the matrix here.  The loss matrix
    # must have zeros on the diagonal and positive numbers
    # elsewhere. It must be the same dimensions as the number of
    # classes.

    lo <- sprintf("matrix(c(%s), byrow=TRUE, nrow=%d)", loss, num.classes)

    if (is.null(parms))
      parms <- sprintf(",\n    parms=list(loss=%s)", lo)
    else
      parms <- gsub(")$", sprintf(",\n    loss=%s)", lo), parms)
  }

  # Build the formula for the model. Rpart has only a formula
  # interface.

  frml <- paste(crs$target, "~ .")

  # Variables to be included --- a string of indicies.

  # included <- getIncludedVariables()
  included <- "c(crs$input, crs$target)" # 20110102

  # Some convenience booleans

  sampling  <- not.null(crs$train)
  including <- not.null(included)
  subsetting <- sampling || including

  # Commands.

  lib.cmd <- "library(rpart, quietly=TRUE)"
  if (! packageIsAvailable("rpart", Rtxt("build decision trees"))) return(FALSE)

  if (action %in%  c("tune", "best"))
  {
    lib.cmd <- paste(lib.cmd, "library(e1071, quietly=TRUE)", sep="\n")
    if (! packageIsAvailable("e1071", Rtxt("tune decision trees"))) return(FALSE)
  }

  # For now, don't use any of the other parameter settings if tune or
  # best. Eventually I want to use the other parameter setting sand
  # override them with the tune options.

  if (action == "build")
  {
    ds.string <- paste("crs$dataset",
                       if (subsetting) "[",
                       if (sampling) "crs$train",
                       if (subsetting) ", ",
                       if (including) included,
                       if (subsetting) "]", sep="")

    rpart.cmd <- paste("crs$rpart <- rpart(", frml, ",\n    data=", ds.string,
                       ifelse(is.null(crs$weights), "",
                              sprintf(",\n    weights=(%s)%s",
                                      crs$weights,
                                      ifelse(sampling, "[crs$train]", ""))),
                       ',\n    method=',
                       ifelse(categoricTarget(),
                              '"class"', '"anova"'),
                       ifelse(is.null(parms), "", parms),
                       ifelse(is.null(control), "", control),
                       ")", sep="")

    print.cmd <- paste("print(crs$rpart)", "printcp(crs$rpart)",
                       'cat("\\n")', sep="\n")

    # 090126 Add error matrix. 100321 Don't add the error matricies -
    # they are more a standard evaluation and belong in the Evaluate
    # tab.

  ##   if (categoricTarget())
  ##   {
  ##     pds.string <- paste("crs$dataset",
  ##                         if (subsetting) "[",
  ##                         if (sampling) "-crs$sample",
  ##                         if (subsetting) ", ",
  ##                         if (including) included,
  ##                         if (subsetting) "]", sep="")
  ##     print.cmd <- paste(print.cmd, "\n",
  ##                        'cat("\\n',
  ##                        ifelse(sampling, "Validation ", "Training "),
  ##                        'dataset error matrix - counts\\n\\n")\n',
  ##                        "print(table(predict(crs$rpart, ",
  ##                        pds.string, ', type="class"),',
  ##                        pds.string, '$', crs$target,
  ##                        ', dnn=c("Predicted", "Actual")))\n',
  ##                        'cat("\\n")\n',
  ##                        sep="")
  ##     print.cmd <- paste(print.cmd, "\n",
  ##                        'cat("\\n',
  ##                        ifelse(sampling, "Validation ", "Training "),
  ##                        'dataset error matrix - percentages\\n\\n")\n',
  ##                        "print(round(100*table(predict(crs$rpart, ",
  ##                        pds.string, ', type="class"),',
  ##                        pds.string, '$', crs$target,
  ##                        ', dnn=c("Predicted", "Actual"))/nrow(',
  ##                        pds.string, ")))\n",
  ##                        'cat("\\n")\n',
  ##                        sep="")
  ##   }
  }
  else if (action == "tune")
  {
    rpart.cmd <- paste("crs$tune.rpart <- tune.rpart(", frml, ", data=crs$dataset",
                       if (subsetting) "[",
                       if (sampling) "crs$train",
                       if (subsetting) ",",
                       if (including) included,
                       if (subsetting) "]",
                       sprintf(", %s", tune.controls),
                       ")", sep="")

    print.cmd <- paste("print(crs$tune.rpart)", "plot(crs$tune.rpart)", sep="\n")
    }
  else if (action == "best")
  {
    # This won't work - best.rpart usese the tune.control() structure
    rpart.cmd <- paste("crs$rpart <- best.rpart(", frml, ", data=crs$dataset",
                     if (subsetting) "[",
                     if (sampling) "crs$train",
                     if (subsetting) ",",
                     if (including) included,
                     if (subsetting) "]",
                     sprintf(", %s", tune.controls),
                     ")", sep="")

    print.cmd <- paste("print(crs$rpart)", "printcp(crs$rpart)", sep="\n")
  }

  # Load the required library.

  startLog(commonName(crv$RPART))
  appendLog(packageProvides("rpart", "rpart"), lib.cmd)

  eval(parse(text=lib.cmd))

  # Set the seed so that xerror and xstd are consistent each time

  seed.cmd <- 'set.seed(crv$seed)'
  appendLog(Rtxt("Reset the random number seed to obtain the same results each time."),
            seed.cmd)
  eval(parse(text=seed.cmd))

  # Build the model.

  appendLog(sprintf(Rtxt("Build the %s model."), commonName(crv$RPART)), rpart.cmd)
  start.time <- Sys.time()
  result <- try(eval(parse(text=rpart.cmd)), silent=TRUE)
  time.taken <- Sys.time() - start.time
  if (inherits(result, "try-error"))
  {
    if (any(grep("syntax error.*weights", result)))
      errorDialog(sprintf(Rtxt("The call to 'rpart' has a syntax error in the",
                               "weights formula. The error message was:\n\n%s"),
                          result))
    else
      errorDialog(errorMessageFun("rpart", result))
    return(FALSE)
  }

  # Summary: Show the resulting model.

  appendLog(sprintf(Rtxt("Generate a textual view of the %s model."),
                    commonName(crv$RPART)), print.cmd)

  resetTextview(TV)
  setTextview(TV,
              sprintf(Rtxt("Summary of the %s model for %s (built using '%s'):"),
                      commonName(crv$RPART),
                      Rtxt("Classification"), # 080604 TODO put the right type
                      "rpart"),
              "\n\n",
              collectOutput(print.cmd))

  if (sampling) crs$smodel <- union(crs$smodel, crv$RPART)

  # Now that we have a model, make sure the rules and plot buttons are
  # visible.

  showModelRPartExists()

  # Finish up.

  reportTimeTaken(TV, time.taken, commonName(crv$RPART))

  return(TRUE)
}

showModelRPartExists <- function(state=!is.null(crs$rpart))
{
  # If an rpart model exists then show the Rules and Draw buttons on
  # the Model tab.

  if (state)
  {
    theWidget("rpart_plot_button")$show()
    theWidget("rpart_plot_button")$setSensitive(TRUE)
    if ("BinaryTree" %in% class(crs$rpart))
      theWidget("rpart_rules_button")$hide()
    else
    {
      theWidget("rpart_rules_button")$show()
      theWidget("rpart_rules_button")$setSensitive(TRUE)
    }
  }
  else
  {
    theWidget("rpart_plot_button")$hide()
    theWidget("rpart_rules_button")$hide()
  }
}

list.rule.nodes.rpart <- function(model)
{
  # The information we need is in the rpart frame
  frm <- model$frame
  # Obtain the probabilities
  nodes <- frm$yval2[,5]
  # Get the probability ordered index of leaf nodes
  ordered <- sort(frm$yval2[,5][frm[,1] == "<leaf>"],
                  decreasing=TRUE, index=TRUE)
  # Return the list of node numbers
  return(row.names(frm)[which(frm[,1] == "<leaf>")][ordered$ix])
}

# drawTreeNodes() is a modification of draw.tree() from maptree and we
# acknowledge the original author Denis White who released the code of
# the maptree package under an Unlimited license, and with no
# copyright claims found in any of the source code files.
#
# The modification is to draw rule nodes rather than sequential
# numbers to label the leaves.
#
# 080315 change size from 2.5 to 4 so the percentages are not cramped
# up---they did not use to but with a new version of R (not sure which
# one) they did.
#
# 080315 Rob Williams noted that when using loss (e.g., 0,4,1,0) the
# percentages being displayed are selecting the larger probabiliy, not
# the correct probability.

drawTreeNodes <- function (tree, cex = par("cex"), pch = par("pch"),
                           size = 4 * cex, col = NULL, nodeinfo = FALSE,
                           units = "", cases = "obs",
                           digits = getOption("digits"),
                           decimals = 2,
                           print.levels = TRUE, new = TRUE)
{
  if (new) plot.new()

  rtree <- length(attr(tree, "ylevels")) == 0
  tframe <- tree$frame
  rptree <- length(tframe$complexity) > 0
  node <- as.numeric(row.names(tframe))
  leafnode <- node[tframe$var == "<leaf>"]
  #proportions <- sprintf("%0.2f", tframe$yval2[,5])
  depth <- floor(log(node, base = 2) + 1e-07)
  depth <- as.vector(depth - min(depth))
  maxdepth <- max(depth)
  x <- -depth
  y <- x
  leaves <- tframe$var == "<leaf>"
  x[leaves] <- seq(sum(leaves))
  depth <- split(seq(node)[!leaves], depth[!leaves])
  parent <- match(node%/%2, node)
  left.child <- match(node * 2, node)
  right.child <- match(node * 2 + 1, node)
  for (i in rev(depth)) x[i] <- 0.5 * (x[left.child[i]] + x[right.child[i]])
  nleaves <- sum(leaves)
  nnodes <- length(node)
  if (rtree)
  {
    dev <- tframe$dev
    pcor <- rep(0, nnodes)
    for (i in 1:nnodes) if (!leaves[i]) {
      l <- dev[node == (node[i] * 2)]
      r <- dev[node == (node[i] * 2 + 1)]
      pcor[i] <- dev[i] - l - r
    }
    pcor <- round(pcor/dev[1], 3) * 100
  }
  else
  {
    crate <- rep(0, nnodes)
    trate <- 0
    if (!rptree)
    {
      for (i in 1:nnodes)
      {
        yval <- tframe$yval[i]
        string <- paste("tframe$yprob[,\"", as.character(yval),
                        "\"]", sep = "")
        crate[i] <- eval(parse(text = string))[i]
        if (leaves[i])
          trate <- trate + tframe$n[i] * crate[i]
      }
    }
    else
    {
      # 131120 Jouanne-Diedrich Holger noted the probabilities in the
      # tree were wrong. Seems like rpart's frame is now returning a
      # yval2 with 2 extra columns!
      
      nlv <- floor((ncol(tframe$yval2)-2)/2)
      for (i in 1:nnodes)
      {
        yval <- tframe$yval[i]
        # [080315 gjw] Now sort the class rates and get the largest!!!
        # But wouldn't we want to get the one corresponding to yval
        # rather than the largest? It won't necessarily be the largest
        # the when we use a loss matrix? The original code from
        # draw.tree uses the largest, but for my purposes I was
        # wanting the correct one! So I keep the index calculation
        # (for notation) but replace it in the crate assignment with
        # yval instead. With this change I don't think I affect too
        # much. The trate variable is affected, but it is only printed
        # with nodeinfo set to TRUE. I'm not sure the original is
        # correct though. This is reported as the "total classified
        # correct."
        index <- rev(order(tframe$yval2[i, 2:(nlv + 1)]))[1]
        # ORIG crate[i] <- tframe$yval2[i, (nlv + 1 + index1)]
        crate[i] <- tframe$yval2[i, (nlv + 1 + yval)]
        if (leaves[i])
        {
          trate <- trate + tframe$n[i] * crate[i]
        }
      }
    }
    crate <- round(crate, 3) * 100
    trate <- round(trate/tframe$n[1], 3) * 100
  }
  if (is.null(col))
    kol <- rainbow(nleaves)
  else if (col == "gray" | col == "grey")
    kol <- gray(seq(0.8, 0.2, length = nleaves))
  else kol <- col
  xmax <- max(x)
  xmin <- min(x)
  ymax <- max(y)
  ymin <- min(y)
  pinx <- par("pin")[1]
  piny <- par("pin")[2]
  xscale <- (xmax - xmin)/pinx
  box <- size * par("cin")[1]
  if (box == 0)
    xbh <- xscale * 0.2
  else xbh <- xscale * box/2
  chr <- cex * par("cin")[2]
  tail <- box + chr
  yscale <- (ymax - ymin)/(piny - tail)
  ytail <- yscale * tail
  if (box == 0)
    ybx <- yscale * 0.2
  else ybx <- yscale * box
  ychr <- yscale * chr
  ymin <- ymin - ytail
  xf <- 0.1 * (xmax - xmin)
  yf <- 0.1 * (ymax - ymin)
  x1 <- xmin - xf
  x2 <- xmax + xf
  y1 <- ymin - yf
  y2 <- ymax + yf
  par(usr = c(x1, x2, y1, y2))
    v <- as.character(tframe$var[1])
    if (rptree) {
        sp <- tree$splits[1, ]
        val <- sp["index"]
        if (sp["ncat"] > 1) {
            r <- sp["index"]
            string <- "attributes(tree)$xlevels$"
            string <- paste(string, v, sep = "")
            xl <- eval(parse(text = string))
            lf <- rf <- ""
            for (k in 1:sp["ncat"]) if (tree$csplit[r, k] == 1)
                lf <- paste(lf, xl[k], sep = ",")
            else rf <- paste(rf, xl[k], sep = ",")
            if (!print.levels)
                string <- v
            else if (nchar(lf) + nchar(rf) > 30) # Avoid too long
              string <- v
            else
              string <- paste(lf, "=", v, "=", rf)

        }
        else {
            if (sp["ncat"] < 0)
                op <- "< =>"
            else op <- ">= <"
            string <- paste(v, op, round(val, decimals))
        }
    }
    else {
        val <- substring(as.character(tframe$splits[1, 1]), 2)
        string <- paste(as.character(v), "<= >", round(val, decimals))
    }
    text.default(x[1], y[1], string, cex = cex)
    if (nodeinfo)
    {
      n <- tframe$n[1]
      if (rtree)
      {
        z <- round(tframe$yval[1], digits)
        r <- pcor[1]
        string <- paste(z, " ", units, "; ", n, " ", cases,
                        "; ", r, "%", sep = "")
      }
      else {
            z <- attr(tree, "ylevels")[tframe$yval[1]]
            r <- crate[1]
            string <- paste(z, "; ", n, " ", cases, "; ", r,
                "%", sep = "")
        }
        text.default(x[1], y[1] - ychr, string, cex = cex)
    }
  for (i in 2:nnodes)
  {
    ytop <- ychr * (as.integer(nodeinfo) + 1)
    if (y[i] < y[i-1])
    {
      lines(c(x[i-1], x[i]), c(y[i-1] - ytop, y[i-1] - ytop))
      lines(c(x[i], x[i]), c(y[i-1] - ytop, y[i] + ychr))
    }
    else
    {
      lines(c(x[parent[i]], x[i]), c(y[parent[i]] - ytop, y[parent[i]] - ytop))
      lines(c(x[i], x[i]), c(y[parent[i]] - ytop, y[i] + ychr))
    }
    if (!leaves[i])
    {
      v <- as.character(tframe$var[i])
      if (rptree)
      {
        k <- 1
        for (j in 1:(i-1))
        {
          m <- tframe$ncompete[j]
          if (m > 0) k <- k + m + 1
          m <- tframe$nsurrogate[j]
          if (m > 0) k <- k + m
        }
        sp <- tree$splits[k, ]
        val <- sp["index"]
        if (sp["ncat"] > 1) {
          r <- sp["index"]
          string <- "attributes(tree)$xlevels$"
          string <- paste(string, v, sep = "")
          xl <- eval(parse(text = string))
          lf <- rf <- ""
          for (k in 1:sp["ncat"])
            if (tree$csplit[r, k] == 1)
              lf <- paste(lf, xl[k], sep = ",")
            else rf <- paste(rf, xl[k], sep = ",")
          if (!print.levels)
            string <- v
          else if (nchar(lf) + nchar(rf) > 10) # Avoid too long
            string <- v
          else
            string <- paste(lf, "=", v, "=", rf)
        }
        else {
          if (sp["ncat"] < 0)
            op <- "< =>"
          else op <- expression(">= <")
          string <- paste(v, op, round(val, decimals))
        }
      }
      else {
        val <- substring(as.character(tframe$splits[i,
                                                    1]), 2)
        string <- paste(as.character(v), "< =>", round(val, decimals))
      }
      text.default(x[i], y[i], string, cex = cex)
      if (nodeinfo) {
        n <- tframe$n[i]
        if (rtree) {
          z <- round(tframe$yval[i], digits)
          r <- pcor[i]
          string <- paste(z, " ", units, "; ", n, " ",
                          cases, "; ", r, "%", sep = "")
        }
        else {
          z <- attr(tree, "ylevels")[tframe$yval[i]]
          r <- crate[i]
          string <- paste(z, "; ", n, " ", cases, "; ",
                          r, "%", sep = "")
        }
        text.default(x[i], y[i] - ychr, string, cex = cex)
      }
    }
    else
    {
      if (box == 0)
      {
        lines(c(x[i], x[i]), c(y[i], y[i] + ychr))
        lines(c(x[i] - xbh, x[i] + xbh), c(y[i], y[i]))
      }
      else
      {
        # points(x[i], y[i], pch = pch, cex = size, col = kol[x[i]])
      }
      if (rtree)
      {
        z <- round(tframe$yval[i], digits)
        text.default(x[i], y[i] - ybx, paste(z, units,
                                             sep = " "), cex = cex)
      }
      else
      {
        z <- attr(tree, "ylevels")[tframe$yval[i]]
        text.default(x[i], y[i] - ybx, z, cex = cex)
      }
      n <- tframe$n[i]
      text.default(x[i], y[i] - ybx - ychr,
                   paste(n, cases, sep = " "), cex=cex)
      if (! rtree)
        text.default(x[i], y[i] - 1.6*ybx - ychr,
                     paste(crate[i], "%", sep = ""), cex=cex)
      #paste(crate[i], "%/", proportions[i], sep = ""), cex = cex)
      if (box != 0)
      {
        # ORIG text.default(x[i], y[i], as.character(x[i]),
        text.default(x[i], y[i], as.character(leafnode[x[i]]),
                     cex = cex, col=kol[x[i]], font=2)
      }
    }
  }
  if (nodeinfo) {
    if (rtree)
      string <- paste("Total deviance explained =", sum(pcor),
                      "%")
    else string <- paste("Total classified correct =", trate,
                         "%")
    if (box == 0)
      text.default(mean(x), ymin - 3.5 * ychr, string, cex = 1.2 *
                   cex)
    else text.default(mean(x), ymin - 2.2 * ybx, string,
                      cex = 1.2 * cex)
  }
}

exportRpartModel <- function()
{
  # Make sure we have a model first!

  if (noModelAvailable(crs$rpart, crv$RPART)) return(FALSE)

  startLog(paste(Rtxt("Export"), commonName(crv$RPART)))

  save.name <- getExportSaveName(crv$RPART)
  if (is.null(save.name)) return(FALSE)
  ext <- tolower(get.extension(save.name))

  # Construct the command to produce PMML.

  pmml.cmd <- sprintf("pmml(crs$rpart%s, dataset=crs$dataset)",
                      ifelse(length(crs$transforms),
                             ", transforms=crs$transforms", ""))

  # We can't pass "\" in a filename to the parse command in MS/Windows
  # so we have to run the save/write command separately, i.e., not
  # inside the string that is being parsed.

  if (ext == "xml")
  {
    appendLog(sprintf(Rtxt("Export %s as PMML."), commonName(crv$RPART)),
              sprintf('saveXML(%s, "%s")', pmml.cmd, save.name))
    XML::saveXML(eval(parse(text=pmml.cmd)), save.name)
  }
  else if (ext == "c")
  {
    # 090103 gjw Move to a function: saveC(pmml.cmd, save.name,
    # "decision tree")

    # 090223 Why is this tolower being used? Under GNU/Linux it is
    # blatantly wrong. Maybe only needed for MS/Widnows

    if (isWindows()) save.name <- tolower(save.name)

    model.name <- sub("\\.c", "", basename(save.name))

    export.cmd <- generateExportPMMLtoC(model.name, save.name, "rpart_textview")
    
    appendLog(sprintf(Rtxt("Export %s as a C routine."), commonName(crv$RPART)),
              sprintf('pmml.cmd <- "%s"\n\n', pmml.cmd),
              export.cmd)

    eval(parse(text=export.cmd))
  }      
  setStatusBar(sprintf(Rtxt("The model has been exported to '%s'."), save.name))
}

