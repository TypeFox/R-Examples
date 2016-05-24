# Gnome R Data Miner: GNOME interface to R for Data Mining
#
# Time-stamp: <2015-11-25 08:31:04 gjw>
#
# Implement EXPLORE functionality.
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

resetExploreTab <- function(new.dataset=TRUE)
{
  cbox1 <- theWidget("pairs_color_combobox")
  cbox1$setSensitive(TRUE)
  
  if (new.dataset)
  {
    vl <- getCategoricVariables(type="names",include.target=T)
    
    if (length(vl))
    {
      cbox1$getModel()$clear()
      cbox1$appendText(" ") # no color
      lapply(vl, cbox1$appendText)
      if ( length(crs$target))
        cbox1$setActive(length(vl)) # by construction we know is the last one
    }
  }
}  

########################################################################
# EXECUTION

executeExploreTab <- function()
{
  
  # Can not explore the data if there is no dataset.

  if (noDatasetLoaded()) return()

  # 080315 Don't proceed if the variable selections have changed. For
  # example, two targets might have been selected, resulting in a
  # popup, but not yet resolved, and so many of the plots will fail.
  
  if (variablesHaveChanged(Rtxt("building a model"))) return()

  # Ensure Sample does not require executing.

  sampling <- theWidget("data_sample_checkbutton")$getActive()
  if (sampling && sampleNeedsExecute()) return()

  # We generate a string representing the subset of the dataset on
  # which the exploration is to be performed. This is then passed to
  # the individually dispatched functions.

  vars <- "c(crs$input, crs$risk, crs$target)" # 20110102 getIncludedVariables(risk=TRUE)
  dataset <- sprintf("%s[%s, %s]", "crs$dataset",
                     ifelse(sampling, "crs$sample", ""),
                     ifelse(is.null(vars),"", vars))

  # For the distribution plot, we do list all variables in the
  # interface, even if they are ignored. TODO 061006 We could instead
  # grey out the ignored ones (i.e., make them not sensitive). But
  # for now, for plots, allow all variables, even the ignored ones,
  # and thus we need a dataset that includes all variables - the
  # "avdataset".

  avdataset <- sprintf("%s[%s,]", "crs$dataset",
                     ifelse(sampling, "crs$sample", ""))
  
  vars <- "crs$numeric" # 20110102 getIncludedVariables(numonly=TRUE)
  
  
  # TODO 060606 The question here is whether NULL means all variables
  # or means none found?
  
  #if (is.null(vars))
  #  ndataset <- NULL
  #else
    ndataset <- sprintf("%s[%s, %s]", "crs$dataset",
                        ifelse(sampling, "crs$sample", ""),
                        ifelse(is.null(vars),"",vars))

  # Numeric input variables

  vars <- "crs$numeric" # 20110102 inputVariables(numonly=TRUE)
  nidataset <- sprintf("%s[%s, %s]", "crs$dataset",
                       ifelse(sampling, "crs$sample", ""),
                       ifelse(is.null(vars),"",vars))
  
  # Dispatch
  
  if (theWidget("summary_radiobutton")$getActive())
    executeExploreSummary(dataset)
  else if (theWidget("explore_distr_radiobutton")$getActive())
  {
    if (crv$appname == "RattleXX")
      executeExplorePlot2(avdataset)
    else
      executeExplorePlot(avdataset)
  }
  else if (theWidget("explore_interactive_radiobutton")$getActive())
  {
#    if (theWidget("explore_interactive_latticist_radiobutton")$getActive())
#      executeExplorePlaywith(dataset)
#    else
      if (theWidget("explore_interactive_ggobi_radiobutton")$getActive())
      executeExploreGGobi(dataset, crs$dataname)
#    else if (theWidget("explore_interactive_plotbuilder_radiobutton")$getActive())
#      executeExplorePlotBuilder()
  }
  else if (theWidget("explore_correlation_radiobutton")$getActive())
  {
    if (theWidget("explore_correlation_hier_checkbutton")$getActive())
      executeExploreHiercor(ndataset)
    else
    {
      if (theWidget("explore_correlation_na_checkbutton")$getActive())
        executeExploreCorrelation(dataset)
      else
        executeExploreCorrelation(ndataset)
    }
  }
  else if (theWidget("prcomp_radiobutton")$getActive())
    executeExplorePrcomp(nidataset)
}

executeExploreSummary <- function(dataset)
{
  TV <- "summary_textview"

  # Get the current state of the relevant buttons.
  
  use.sample  <- theWidget("data_sample_checkbutton")$getActive()
  do.summary  <- theWidget("summary_checkbutton")$getActive()
  do.describe <- theWidget("describe_checkbutton")$getActive()
  do.basics   <- theWidget("basics_checkbutton")$getActive()
  do.kurtosis <- theWidget("kurtosis_checkbutton")$getActive()
  do.skewness <- theWidget("skewness_checkbutton")$getActive()
  do.missing  <- theWidget("missing_checkbutton")$getActive()
  do.crosstab <- theWidget("explore_crosstab_checkbutton")$getActive()

  # Make sure something has been selected.
  
  if (! (do.summary || do.describe || do.basics ||
         do.kurtosis || do.skewness || do.missing ||
         do.crosstab))
  {
    infoDialog(Rtxt("No summary type has been selected.",
                    "Please click at least one checkbox."))
    return()
  }
    
  # Other useful information:
  #   is there a sample
  #   list of numeric variables
  
  sampling  <- not.null(crs$sample)

  numeric.cmd <- sprintf(paste("seq(1,ncol(%s))",
                               "[as.logical(sapply(%s, is.numeric))]",
                               sep=""), dataset, dataset)
  nvars <- simplifyNumberList(eval(parse(text=numeric.cmd)))

  # Start the trace to the log.
  
  startLog()
  theWidget(TV)$setWrapMode("none")
  resetTextview(TV)

  # Construct and execute the requested commands.

  if (do.summary)
  {
    # Find the number of observations with any missing value for the
    # non-ignored variables.
    
    missing.cmd <- sprintf('length(attr((na.omit(%s)), "na.action"))', dataset)
    result <- try(missing <- eval(parse(text=missing.cmd)), silent=TRUE)
    if (inherits(result, "try-error")) missing <- 0
    
    # Use Hmisc's contents to summarise the data frame, if Hmisc is
    # available.

    contents.cmd <- ""
    if (packageIsAvailable("Hmisc", Rtxt("describe the contents of a data frame")))
    {
      lib.cmd <- "library(Hmisc, quietly=TRUE)"
      appendLog(packageProvides("Hmisc", "contents"), lib.cmd)
      eval(parse(text=lib.cmd))
      contents.cmd <- sprintf("contents(%s)", dataset)
    }
    summary.cmd <- sprintf("summary(%s)", dataset)
    
    appendLog(Rtxt("Obtain a summary of the dataset."), contents.cmd, "\n", summary.cmd)
    appendTextview(TV,
                   Rtxt("Below we summarise the dataset."), "\n",
                   ifelse(use.sample && sampling,
                          paste("\n",
                                Rtxt("The data is limited to the training dataset."),
                                "\n", sep=""), ""),
                   if (missing > 0)
                   sprintf(Rtxt("Note that the data contains %d observations",
                                "with missing values.",
                                "\nEnable the 'Show Missing'",
                                "check box for details."),
                           missing),
                   "\n",
                   collectOutput(contents.cmd),
                   "\n\n",
                   Rtxt("For the simple distribution tables below the 1st and 3rd Qu.",
                        "\nrefer to the first and third quartiles, indicating that 25%",
                        "\nof the observations have values of that variable which are",
                        "\nless than or greater than (respectively) the value listed."),
                   "\n\n",
                   collectOutput(summary.cmd))
  }

  if (do.describe)
  {
    ## A different summary, using Hmisc's describe.
  
    if (packageIsAvailable("Hmisc", Rtxt("describe the data")))
    {
      lib.cmd <- "library(Hmisc, quietly=TRUE)"
      appendLog(packageProvides("Hmisc", "describe"), lib.cmd)
      eval(parse(text=lib.cmd))
      
      describe.cmd <- sprintf("describe(%s)", dataset)
      appendLog(Rtxt("Generate a description of the dataset."), describe.cmd)
      appendTextview(TV,
                     Rtxt("Below is a description of the dataset."),
                     ifelse(use.sample && sampling,
                            Rtxt("\nThe data is limited to the training dataset."), ""),
                     "\n\n",
                     collectOutput(describe.cmd, TRUE, width=200))
    }
  }

  if (do.basics || do.kurtosis || do.skewness)
  {
    ## These all require the fBasics library, so check only once.
    
    if (packageIsAvailable("fBasics", Rtxt("calculate basic stats, skew and kurtosis")))
    {
      lib.cmd <- "library(fBasics, quietly=TRUE)"
      
      if (do.basics)
      {
        appendLog(packageProvides("basicStats", "fBasics"), lib.cmd)
        eval(parse(text=lib.cmd))
        basics.cmd <- sprintf("lapply(%s[,%s], basicStats)", dataset,
                              ifelse(is.null(nvars), "", nvars))
        appendLog(Rtxt("Generate a description of the numeric data."), basics.cmd)
        appendTextview(TV,
                       Rtxt("Basic statistics for each numeric variable",
                            "of the dataset."),
                       "\n\n",
                       collectOutput(basics.cmd, TRUE))
      }
      
      if (do.kurtosis)
      {
        appendLog(packageProvides("kurtosis", "fBasics"), lib.cmd)
        eval(parse(text=lib.cmd))
        kurtosis.cmd <- sprintf("kurtosis(%s[,%s], na.rm=TRUE)", dataset,
                                ifelse(is.null(nvars), "", nvars))

        appendLog(Rtxt("Summarise the kurtosis of the numeric data."), kurtosis.cmd)
        appendTextview(TV,
                       Rtxt("Kurtosis for each numeric variable of the dataset.",
                            "\nLarger values mean sharper peaks and flatter tails.",
                            "\nPositive values indicate an acute peak around the mean.",
                            "\nNegative values indicate a smaller peak around",
                            "the mean."),
                       "\n\n",
                       collectOutput(kurtosis.cmd, TRUE))
      }

      if (do.skewness)
      {
        appendLog(packageProvides("skewness", "fBasics"), lib.cmd)
        eval(parse(text=lib.cmd))
        skewness.cmd <- sprintf("skewness(%s[,%s], na.rm=TRUE)", dataset,
                                ifelse(is.null(nvars), "", nvars))

        appendLog(Rtxt("Summarise the skewness of the numeric data."), skewness.cmd)
        appendTextview(TV,
                       Rtxt("Skewness for each numeric variable of the dataset.",
                            "\nPositive means the right tail is longer."),
                       "\n\n",
                       collectOutput(skewness.cmd, TRUE))
      }
    }
  }
      
  if (do.missing)
  {
    ## Add in a summary of the missing values.
  
    if (packageIsAvailable("mice", Rtxt("summarise missing values")))
    {
      ## Load the mice package into the library

      lib.cmd <- "library(mice, quietly=TRUE)"
      appendLog(packageProvides("md.pattern", "mice"), lib.cmd)
      eval(parse(text=lib.cmd))
      
      ## Variables to be included, as a string of indicies.
  
      included <- "c(crs$input, crs$target)" # 20110102 getIncludedVariables()
      including <- not.null(included)

      ## Add radio buttons to choose: Full, Train, Test dataset to summarise.
  
      ## Build the summary command

      summary.cmd <- paste("md.pattern(crs$dataset[,",
                           if (including) included,
                           "])", sep="")

      appendLog(Rtxt("Generate a summary of the missing values in the dataset."),
                summary.cmd)
      appendTextview(TV,
                     Rtxt("Missing Value Summary"),
                     "\n\n",
                     collectOutput(summary.cmd, TRUE))
    }
  }

  if (do.crosstab && length(getCategoricVariables()))
  {
    if (packageIsAvailable("descr", Rtxt("produce a cross tabulation")))
    {
      lib.cmd <- "library(descr, quietly=TRUE)"
      appendLog(packageProvides("CrossTable", "descr"), lib.cmd)
      eval(parse(text=lib.cmd))
      crosstab.cmd <- paste(sprintf("for (i in %s)", getCategoricVariables()),
                            "\n{",
                            "\n",
                            sprintf(paste(" cat(sprintf('%s',",
                                          "names(crs$dataset)[i], crs$target))"),
                                    Rtxt("CrossTab of %s by target variable %s\\n\\n")),
                            "\n  print(CrossTable(crs$dataset[[i]],",
                            "crs$dataset[[crs$target]], expected=TRUE, format='SAS'))",
                            "\n  cat(paste(rep('=', 70), collapse=''), '\n\n')",
                            "\n}")
      appendLog(Rtxt("Generate cross tabulations for categoric data."),
                crosstab.cmd)
      appendTextview(TV,
                     Rtxt("Cross tabulations:"),
                     "\n\n",
                     collectOutput(crosstab.cmd, TRUE))
      
    }
  }
  
  # Report completion to the user through the Status Bar.
  
  setStatusBar(Rtxt("Data summary generated."))
}

getVariableIndicies <- function(variables)
{
  indicies <- NULL
  if (not.null(variables))
    indicies <- unlist(lapply(variables, match, colnames(crs$dataset)))
  return(indicies)
}

calcInitialDigitDistr <- function(l,
                                  digit=1,
                                  len=1,
                                  sp=c("none", "positive", "negative"))
{

  # From a list of numbers return a vector of first digit
  # frequencies. If DIGIT is given, then return the distribution for
  # that digit, rather than the default first digit. 130821 xtend to
  # allow LEN to be specified as the number of digits (1 or 2) to
  # calculate the distribution for.

  # The default SP is none, meaning that both positive and negative
  # numbers are considered (ignoring the sign). Otherwise we return
  # the distribution for either only the positive numbers in the list
  # or for only the negative numbers in the list.

  if (sp == "positive")
    l <- l[l>0]
  else if (sp == "negative")
    l <- l[l<0]

  # Ignore all zeros.

  l <- l[l!=0]
  
  # If we don't have any numbers in the distrbution, return a list of
  # zeros.
  
  if (length(l) == 0)
  {
    if (digit == 1)
    {
      result <- rep(0, 9)
      names(result) <- 1:9
    }
    else
    {
      result <- rep(0, 10)
      names(result) <- 0:9
    }
    return(result)
  }

  # Note that we remove decimal points (i.e., the decimal dot itself,
  # not any digits) from real numbers.
  
  ds <- data.frame(digit=as.numeric(substr(gsub("\\.", "",
                     as.character(abs(l))), digit, digit+len-1)),
                   value=1)
#[071201  ds <- data.frame(digit=as.numeric(gsub("(.).*", "\\1",
#                     as.character(abs(l)))),
#                   value=1)
  
  # Ignore any zeros
  
  if (digit == 1) ds <- ds[ds$digit!=0,]

  # Add in any mising digits as value=0

  if (len == 1)
  {
    missing <- setdiff(ifelse(digit>1,0,1):9, unique(ds[,1]))

    # Don't do it for 2 digit - we want to just plot the actual ones,
    # not all possible combinations.
    #
    # setdiff(as.character(paste(ifelse(digit>1, 0, 1):9,
    #                            rep(0:9, ifelse(digit>1, 10, 9)),
    #                            sep="")),
    #         as.character(unique(ds[,1])))
    #
    if (length(missing) > 0)
      ds <- rbind(ds, data.frame(digit=missing, value=0))
  }
  dsb <- by(ds, as.factor(ds$digit), function(x) sum(x$value))
  return(as.matrix(dsb)[,1]/sum(as.matrix(dsb)[,1]))
}

#130825 Is this actually used anywhere?
#
## plotBenfordsLaw <- function(l)
## {
##   if (! packageIsAvailable("gplots", Rtxt("plot Benford's law"))) return()
##   library(gplots, quietly=TRUE)
  
##   actual <- calcInitialDigitDistr(l)
  
##   x  <- 1:9
##   expect <- log10(1 + 1/x)
  
##   nds <- t(as.matrix(data.frame(expect=expect, actual=actual)))

##   ttl <- genPlotTitleCmd(Rtxt("Benford's Law"), vector=TRUE)
  
##   barplot2(nds, beside=TRUE, main=ttl[1], sub=ttl[2],
##            xlab=Rtxt("Initial Digit"), ylab=Rtxt("Probability"))
## }

#-----------------------------------------------------------------------
# 130824 Implement Benfords plot using ggplot2 for Advanced Graphics,
# and add in the length option.

digitDistr <- function(x, digit=1, len=1, name="freq")
{
  # From the numbers in x return a data frame of digit
  # frequencies. The specified DIGIT is analysed (default is the first
  # digit). For first digit analysis, LEN can be provided as the
  # number of digits to analyse.

  # Ignore zeros.
  
  x <- x[x!=0]
  
  # If we don't have any numbers in the distrbution, return a list of
  # zeros.
  
  if (!length(x))
  {
    if (digit == 1)
      x <- data.frame(digit=1:9, freq=rep(0, 9))
    else
      x <- data.frame(digit=0:9, freq=rep(0, 10))
  }
  else
  {
    # Remove decimal points, retaining only digits.
  
    x <- data.frame(digit=as.numeric(substr(gsub("\\.", "",
                      as.character(abs(x))), digit, digit+len-1)),
                    value=1)
    
    # Add in any mising digits as value=0
    
    if (len == 1)
      missing <- setdiff(ifelse(digit>1,0,1):9, unique(x[,1]))
    else
      missing <- setdiff(as.character(paste(ifelse(digit>1, 0, 1):9,
                                            rep(0:9, ifelse(digit>1, 10, 9)),
                                            sep="")),
                         as.character(unique(x[,1])))
    
    if (length(missing) > 0)
      x <- rbind(x, data.frame(digit=missing, value=0))
    
    x <- plyr::ddply(x, plyr::.(digit), function(y) sum(y$value))
    x$V1 <- x$V1/sum(x$V1)
  }

  names(x)[2] <- name
  return(x)
}

benfordDistr <- function(digit, len=1)
{
  if (digit != 1 && len != 1)
    stop("len must be 1 for other than first digit")
  
  if (digit == 1)
  {
    digits <- seq(10^(len-1), (10^len)-1)
    expect <- log10(1 + 1/digits)
  }
  else
  {
    digits <- 0:9    
    expect <- unlist(lapply(digits, function(x)
                            sum(log10(1 + 1/(10*(seq(10^(digit-2),
                                                     (10^(digit-1))-1))
                                             + x)))))
  }
  
  result <- data.frame(digit=digits, Benford=expect)
  return(result)
}

plotDigitFreq <- function(ds)
{
  dsm <- reshape::melt(ds, id.vars="digit")
  len <- nchar(as.character(ds[1,1]))
  
  p <- ggplot2::ggplot(dsm, ggplot2::aes_string(x="digit", y="value", colour="variable"))
  p <- p + ggplot2::geom_line()
  if (len < 3) p <- p + ggplot2::geom_point()
  if (substr(as.character(ds[1,1]), 1, 1)=="1")
    p <- p + ggplot2::scale_x_continuous(breaks=seq(10^(len-1), (10^len)-1, 10^(len-1)))
  else
    p <- p + ggplot2::scale_x_continuous(breaks=seq(0, 9, 1))
  p <- p + ggplot2::ylab("Frequency") + ggplot2::xlab("Digits")
  p <- p + ggplot2::theme(legend.title=ggplot2::element_blank())
  return(p)
}

executeBenfordPlot2 <- function(dname, benplots, target, targets, stratify, sampling, pmax)
{
  if (! packageIsAvailable("ggplot2", Rtxt("plot a Benford's distribution")) ||
      ! packageIsAvailable("reshape", Rtxt("arrange data for Benford's distribution")))
    return(FALSE)
  
  startLog("Benford's Law")

  lib.cmd <- "library(ggplot2, quietly=TRUE)"
  appendLog(packageProvides("ggplot2", "ggplot"), lib.cmd)
  eval(parse(text=lib.cmd))

  lib.cmd <- "library(reshape, quietly=TRUE)"
  appendLog(packageProvides("reshape", "melt"), lib.cmd)
  eval(parse(text=lib.cmd))

  absbutton <- theWidget("benford_abs_radiobutton")$getActive()
  posbutton <- theWidget("benford_pos_radiobutton")$getActive()
  negbutton <- theWidget("benford_neg_radiobutton")$getActive()
  digit <- theWidget("benford_digits_spinbutton")$getValue()
  len   <- theWidget("benford_num_digits_spinbutton")$getValue()

  if (len > 1 && digit != 1)
  {
    errorDialog(Rtxt("Multiple digit analysis only supported from",
                     "the first digit. For number of digits > 1",
                     "please choose 1 as the starting digit."))
    return()
  }
  
  for (var in benplots)
  {
    newPlot()
    
    # We don't actually eval this command - just for information in
    # the Log.
    initialise.cmd <- paste(ifelse(length(target),
                                   paste('target <- "', target, '"\n', sep=""),
                                   ""),
                            'var    <- "', var, '"\n',
                            "digit  <- ", digit, "\n",
                            "len    <- ", len, sep="")
    appendLog(Rtxt("Initialies the parameters."), initialise.cmd)

    ds.cmd <- paste("ds <- merge(benfordDistr(digit, len),\n",
                    "            digitDistr(", dname, '[var], digit, len, "All"))',
                    sep="")
    if (length(target))
      ds.cmd <- paste(ds.cmd, "\n",
                      "for (i in unique(", dname, "[[target]]))\n",
                      "  ds <- merge(ds, digitDistr(", dname, "[", dname,
                      "[target]==i, var], digit, len, i))", sep="")
    appendLog("Build the dataset", ds.cmd)
    eval(parse(text=ds.cmd))

    title <- paste("Digital Analysis of",
                   switch(digit, "First", "Second", "Third",
                          paste(digit, "th", sep="")),
                   ifelse(len>1, "2 Digits", "Digit"), "\nof", var,
                   ifelse(length(target), paste("by", target), ""))
    plot.cmd <- paste("p <- plotDigitFreq(ds)\n",
                      'p <- p + ggtitle("', title, '")\n',
                      'print(p)', sep="")
    appendLog("Plot the digital distribution", plot.cmd)
    eval(parse(text=plot.cmd))
  }    
}
  

if (FALSE)
{
  target <- "TARGET_Adjusted"
  var    <- "Income"
  digit <- 1
  len <- 1
  ds <- read.csv(system.file("csv/audit.csv", package="rattle"))
  x <- merge(benfordDistr(digit, len),
             digitDistr(ds[var], digit, len, "All"))
  for (i in unique(ds[[target]]))
    x <- merge(x, digitDistr(ds[ds[target]==i, var], digit, len, i))

  plotDigitFreq(x)
}

#-----------------------------------------------------------------------
# Generate a title for plots.

generateTitleText <- function(var, target, sampling, doby)
{
  if (sampling)
    if (doby)
      title.txt <- sprintf(Rtxt("Distribution of %s (sample)\\nby %s"),
                           var, target)
    else
      title.txt <- sprintf(Rtxt("Distribution of %s (sample)"), var)
  else
    if (doby)
      title.txt <- sprintf(Rtxt("Distribution of %s\\nby %s"), var, target)
    else
      title.txt <- sprintf(Rtxt("Distribution of %s"), var)
  return(title.txt)
}

#-----------------------------------------------------------------------
# Original version.

executeExplorePlot <- function(dataset,
                               boxplots = getSelectedVariables("boxplot"),
                               hisplots = getSelectedVariables("hisplot"),
                               cumplots = getSelectedVariables("cumplot"),
                               benplots = getSelectedVariables("benplot"),
                               barplots = getSelectedVariables("barplot"),
                               dotplots = getSelectedVariables("dotplot"),
                               mosplots = getSelectedVariables("mosplot"),
                               paiplots = getSelectedVariables("paiplot"),
                               stratify=TRUE, sampling=NULL,
                               target=crs$target, newplot=TRUE)
{
  # Plot the data. The DATASET is a character string that defines the
  # dataset to use. Information about what variables to plot and the
  # kind of plots is obtained from the continuous_treeview and the
  # categorical_treeview which are displayed in the Explore tab's
  # Distribution option. The appropriate plots are displayed. 090323
  # By having the list of varaiables to display for each type of plot
  # set in the parameter list we can call this function to plot
  # variables from the command line, as in using Sweave. The function
  # remains an internal Rattle function though - only for those who
  # know!

  # 150916 Move to allowing the stratification to be specified here
  # rather than relying on going back to the data tab to change.
    
  target <- theWidget("pairs_color_combobox")$getActiveText()
  if (length(target) && target == " ") target <- NULL
    
  # Has the advanced graphics (ggplot2) option be enabled?

  advanced.graphics <- theWidget("use_ggplot2")$getActive()

  # Obtain the selection of variables.

  nboxplots <- length(boxplots)
  nhisplots <- length(hisplots)
  nbenplots <- length(benplots)
  nbarplots <- length(barplots)
  ndotplots <- length(dotplots)
  nmosplots <- length(mosplots)
  npaiplots <- as.integer(length(paiplots) >= 2) # only one plot is generated, 2 variables minimum are needed

  if (npaiplots == 0 && length(paiplots)>0)
  {
    errorDialog(Rtxt("The pairs plot requires at least two variables."))
    return()
  }
  
  total.plots <- nboxplots + nhisplots + length(cumplots) +
    nbenplots + nbarplots + ndotplots + nmosplots + npaiplots
  
  pmax <- theWidget("plots_per_page_spinbutton")$getValue()
  pcnt <- 0

  # Iterate over all target values if a target is defined and has
  # less than 10 values. The plots will then also display the
  # distributions per target value.

  # 091011 Move to using the value of crs$target instead of getting it
  # from the interface. Eventually, pass target in as an argument so
  # we can be independent of the GUI?

  # target <- getSelectedVariables("target")

  if (length(target))
    targets <- levels(as.factor(crs$dataset[[target]]))
  else
    targets <- NULL

  if (length(targets) > 10) targets <- NULL

  # For now, let's plot always, since I was wondering why the Benford
  # plot was not showing all the targets!
  
##   if (length(targets) > 10)
##   {
##     target <- NULL
##     targets <- NULL
##   }
  
  # 091011 Check if there are mosaic plots requested but there is not
  # target - notify that the mosaic plots will not be plotted.

  if (nmosplots > 0 && is.null(target))
  {
    infoDialog(Rtxt("A mosaic plot can not be displayed without identifying a target",
                    "variable. The requested mosaic plots will be ignored."))
    total.plots <- total.plots - nmosplots
    mosplots <- NULL
    nmosplots <- 0
  }

  # Don't waste real estate if we are plotting less than number
  # allowed per page.
  
  if (total.plots < pmax) pmax <- total.plots
  
  # Check for sampling.

  if (is.null(sampling))
  {
    use.sample <- theWidget("data_sample_checkbutton")$getActive()
    sampling  <- use.sample && not.null(crs$sample)
  }

  # Record other options.

  annotate <- theWidget("explot_annotate_checkbutton")$getActive()
  
  # Split the data, first for all values.

  bind.cmd <- sprintf('rbind(data.frame(dat=%s[,"%%s"], grp="All")', dataset)

  for (i in seq_along(targets))
  {
    bind.cmd <- sprintf("%s,\n            data.frame(dat=%s",
                        bind.cmd, dataset)
    
    bind.cmd <- sprintf('%s[crs$dataset%s$%s=="%s","%%s"], grp="%s")',
                        bind.cmd,
                        ifelse(sampling, "[crs$sample,]", ""),
                        target, targets[i], targets[i])
  }
  
  # Finish off the command to create the dataset for plotting.
  
  bind.cmd <- sprintf("%s)", bind.cmd)

  # Build a list of generic datasets. This describes how to get the
  # relevant rows from the dataset for All the data, then each of the
  # levels of a target. Each contains a "%s" which is replace gor
  # specific chosen variables at the time of using this construct to
  # obtain the data for the plot. The form is:
  #
  # All = crs$dataset$%s
  #
  # or if sampling is enabled:
  #
  # All = crs$dataset[crs$sample,]$%s  
  #
  # For each level:
  #
  # '0' = crs$dataset[crs$dataset$Adjusted=="0",]$%s
  #
  # or if sampling is enabled:
  #
  # '0' = crs$dataset[crs$sample,][crs$dataset[crs$sample,]$Adjusted=="0",]$%s
  #
  # This is a newer alternative to identifying the dataset
  # segments. We build this list of target and a specification of the
  # correspending data subset. Eventually move all plotting to use
  # this approach rather than using bind.cmd.

  genericDataSet <- data.frame(All=sprintf('%s$%%s', dataset))
  for (i in seq_along(targets))
  {
    tmpDataSet <- data.frame(New=sprintf('%s[crs$dataset%s$%s=="%s",]$%%s',
                               dataset,
                               ifelse(sampling, "[crs$sample,]", ""),
                               target, targets[i]))
    colnames(tmpDataSet) <-  c(targets[i])
    genericDataSet <- cbind(genericDataSet, tmpDataSet)
  }

  # 120205 Even newer way of dealing with datasets: use "with" to
  # avoid using "crs$" over and over again. Results in much cleaner
  # looking code.

  new.dataset <- gsub("crs\\$", "", dataset)
  
  # Generate a plot for each variable. If there are too many
  # variables, ask the user if we want to continue.

  if (total.plots > 10 && pmax == 1)
    if (! questionDialog(sprintf(Rtxt("We are about to generate %d",
                                      "individual plots. That's quite a few.",
                                      "You could select fewer variables, or you",
                                      "can change the number of plots per page,",
                                      "but you can also proceed if you like.",
                                      "\n\nWould you like to proceed?"),
                                 total.plots)))
      return()

  #---------------------------------------------------------------------
  # 091005 If no plots are specified then generate a matrix of scatter
  # plots following the example in the pairs function.
  
  if (total.plots == 0)
    if (advanced.graphics)
      displayPairsPlot2(dataset)
    else
      displayPairsPlot(dataset)

  #---------------------------------------------------------------------

  if (nboxplots > 0)
  {
    # Show a box plot for numeric data. A box plot shows the
    # distribution of numeric data graphically. The box iteself
    # extends from the lower to the upper quartiles with the median
    # drawn in the box. The lines then extend to the maximum and
    # minimum points that are no more than 1.5 times the interquartile
    # range from the median. Outliers are then also plotted as
    # points. The notches indicate significant differences, in that if
    # nocthes do not overlap, then the distribution medians are
    # significantly different.")

    if (advanced.graphics &&
        packageIsAvailable("ggplot2", Rtxt("plot using ggplot2")))
    {
      executeBoxPlot2(new.dataset, boxplots, target, targets, stratify, sampling, pmax)
    }
    else
    {

      # 080918 Use the colorspace package to get a better colour
      # map. See
      # http://epub.wu-wien.ac.at/dyn/virlib/wp/eng/showentry?ID=epub-wu-01_c87
      # and
      # http://statmath.wu.ac.at/~zeileis/papers/Zeileis+Hornik+Murrell-2009.pdf
    
    if (packageIsAvailable("colorspace"))
      cols <- "col=colorspace::rainbow_hcl(%d)," # 090524, start = 270, end = 150),"
    else
      cols <- "col=rainbow(%d),"
    
    plot.cmd <- paste('bp <<- boxplot(formula=dat ~ grp, data=ds,',
                      sprintf(cols, length(targets)+1),
                      ifelse(is.null(targets), "",
                             sprintf('xlab="%s",', target)),
                      'ylab="%s",',
                      'varwidth=TRUE,',
                      'notch=TRUE)')

    # Based on an example from Jim Holtman on r-help 070406.
    
    annotate.cmd <- paste("for (i in seq(ncol(bp$stats)))",
                          "{text(x=i,",
                          "y=bp$stats[,i] - 0.02*(max(ds$dat, na.rm=TRUE)",
                          "- min(ds$dat, na.rm=TRUE)),",
                          "labels=bp$stats[,i])}",
                          "\nif (length(bp$out))",
                          "text(x=bp$group+0.1, y=bp$out, labels=bp$out, cex=0.5)")
    
    lib.cmd <- "library(doBy, quietly=TRUE)"
    
    # TODO: Try using "by" instead of needing another package to
    # provide summaryBy. Also, the new version of doBy (061006) seems
    # to be outputting extra status information that makes the R
    # Console a little chatty unneccessarily - perhaps this will
    # disappear again - it looks like debugging information!
    #
    # status:
    # lhsvar     : dat 
    # rhsvar     : grp 
    # idvar      :  
    # fun.names  : mean 
    # varPrefix  : mean 
    # newNames   : mean.dat 

    # Only use summaryBy if there is a target, because it fails if
    # there is actually only one group in the data. Might be a new
    # bug in the doBy package.
    
    if (length(targets) > 1)
      mean.cmd <- paste(sprintf("points(1:%d,", length(targets)+1),
                        "summaryBy(dat ~ grp, data=ds,",
                        "FUN=mean, na.rm=TRUE)$dat.mean,",
                        "pch=8)")
    else
      mean.cmd <- paste(sprintf("points(1:%d,", length(targets)+1),
                        "mean(ds$dat, na.rm=TRUE),",
                        "pch=8)")
    
    for (s in seq_len(nboxplots))
    {

      startLog(Rtxt("Box Plot"))
      cmd <- paste("sprintf(bind.cmd,",
                   paste(paste('"', rep(boxplots[s], length(targets)+1), '"',
                               sep=""),
                         collapse=","),
                   ")")
      cmd <- eval(parse(text=cmd))
      appendLog(sprintf(Rtxt("Subset the data for a boxplot of the variable %s",
                             " with the appropriate groups."),
                        boxplots[s]),
                paste("ds <-", cmd))
      ds <- eval(parse(text=cmd))

      if (newplot && pcnt %% pmax == 0) newPlot(pmax)
      pcnt <- pcnt + 1

      this.plot.cmd <- sprintf(plot.cmd, boxplots[s])
      appendLog(Rtxt("Plot the data, grouped appropriately."),
                gsub("<<", "<", this.plot.cmd))
      eval(parse(text=this.plot.cmd))

      # Add a value for the mean to each boxplot.
      
      if (packageIsAvailable("doBy", Rtxt("add means to box plots")))
      {
        appendLog(packageProvides("doBy", "summaryBy"), lib.cmd)
        eval(parse(text=lib.cmd))

        appendLog(Rtxt("Calculate the group means."), mean.cmd)
        eval(parse(text=mean.cmd))
      }
        
      # Optionally include annotations.

      if (annotate)
      {
        appendLog(Rtxt("Add annotations to the plot."), annotate.cmd)
        eval(parse(text=annotate.cmd))
      }        
      
      # Add a title to the plot.
      
      title.cmd <- genPlotTitleCmd(generateTitleText(boxplots[s], target, sampling,
                                                     stratify && length(targets)))
      
      appendLog(Rtxt("Add a title to the plot."), title.cmd)
      eval(parse(text=title.cmd))
    }
  }
  }

  ##--------------------------------------------------------------------
  
  if (nhisplots > 0)
  {
    # Plot a histogram for numeric data. 090523 Bob Muenchen suggested
    # using a fixed colour for the bars and use the colours as in the
    # box plots for the breakdown between class values, for
    # consistency. I have also introduced colour for the rug and
    # ensured the y limits are calcuated in case the density plot is
    # higher than the histogram. The current colours are not yet right
    # though. For binary data the distinction between the blue and
    # green is very difficult to see for thin lines. 

    if (advanced.graphics &&
        packageIsAvailable("ggplot2", Rtxt("plot using ggplot2")) &&
        packageIsAvailable("dplyr", Rtxt("select from the dataset")))
    {
      executeHistPlot2(new.dataset, hisplots, target, targets, stratify, sampling, pmax)
    }
    else
    {
    
    if (packageIsAvailable("colorspace")) # 090524 Why vcd? comes from colorspace....
      # cols <- "col=rainbow_hcl(%s, start = 270, end = 150)"
      cols <- "col=colorspace::rainbow_hcl(%s)"# 090524, start = 0, end = 150)"
    else
      cols <- "col=rainbow(%s)"

    # 090523 The preplot.cmd was considered for use in determining the
    # highest point of the desnity plots so that we can ensure we have
    # enough on the y axis to draw them all. We could use xpd=NA so as
    # not to crop, but then the y axis may not go up high enough and
    # the density may still be higher than fits in the device, and it
    # looks ugly.
    #
    # 090524 An alternative could be, as with much of Rattle, that
    # this cropping is good enough, and the user can fine tune it
    # themselves. But even with xpd=NA the lines look ugly. Also note
    # that xpd=NA makes the left and right extensions of the desnity
    # plots visible, which is not so nice.
    #
    # The original attempt, whereby the below hist (with plot=FALSE)
    # was included in a preplot.cmd, did not work because somehow this
    # plot actually interferred with multiple plots on a page.
    #
    # Note this tip that could help:
    # http://wiki.r-project.org/rwiki/doku.php?id=tips:easier:tbtable
    #

    ## hs <- hist(ds[ds$grp=="All",1], breaks="fd", plot=FALSE)
    ## dens <- density(ds[ds$grp=="All",1], na.rm=TRUE)
    ## rs <- max(hs$counts)/max(dens$y)
    ## maxy <- max(dens$y*rs)
    ## print(maxy)
    
    ## if (length(targets))
    ##   for (i in 1:length(targets))
    ##   {
    ##     dens <- density(ds[ds$grp==targets[i], 1], na.rm=TRUE)
    ##     maxy <- max(c(maxy, dens$y*rs))
    ##   }
    ## print(maxy)

    # 090524 A note on using breaks in the histogram: The default is
    # to use Sturges' rule from 1926. Rob Hyndman
    # (http://robjhyndman.com/papers/sturges.pdf) argues that Sturges'
    # rule is wrong. Sturges' rule leads to oversmoothed
    # histograms. Hyndman notes that Scott's (1979) rule and Freedman
    # and Diaconis's (1981) rule are just as simple to use as Sturges'
    # rule, but are well-founded in statistical theory. Sturges' rule
    # does not work for large n. Also see
    # http://www.math.leidenuniv.nl/~gill/teaching/statistics/histogram.pdf
    # So let's use the Freedman and Diaconis approach.
    
    plot.cmd <- paste('hs <- hist(ds[ds$grp=="All",1], main="", xlab="%s", ',
                      sprintf('ylab="%s", ', Rtxt("Frequency")),
                      # cols,
                      'col="grey90"',
                      ', ylim=c(0, %s)', #', ceiling(maxy), ')',
                      ', breaks="fd", border=TRUE)\n',
                      'dens <- density(ds[ds$grp=="All",1], na.rm=TRUE)\n',
                      # 090523 Now done in preplot.cmd - not any more? 091021
                      'rs <- max(hs$counts)/max(dens$y)\n',
                      'lines(dens$x, dens$y*rs, type="l", ',
                      sprintf(cols, length(targets)+1), '[1])',
                      sep="")
    if (stratify && length(targets))
    {
      plot.cmd <- paste(plot.cmd, "\n",
                        paste(sprintf(paste('dens <- density(ds[ds$grp=="%s",',
                                            '1], na.rm=TRUE)\n',
                                            'rs <- max(hs$counts)/max(dens$y)\n',
                                            'lines(dens$x, dens$y*rs, ',
                                            'type="l", ',
                                            '%s[%s])', sep=""),
                                      targets,
                                      sprintf(cols, length(targets)+1),
                                      (seq_along(targets)+1)),
                                      #eval(parse(text=sprintf(cols,
                                      #             length(targets))))),
                              collapse="\n"),
                        sep="")
    }

    # 110213 Remove colour from the rug plot as it is misleading when
    # values are identical and there is a lot of overplotting.
    
#    if (stratify && length(targets))
#    {
#      rug.cmd <- paste(sprintf('rug(ds[ds$grp=="%s", 1], %s[%s])', targets,
#                               sprintf(cols, length(targets)+1),
#                               seq_along(targets)+1), collapse="\n")
#    }
#    else
#    {
      rug.cmd <- 'rug(ds[ds$grp=="All",1])'
#    }

    # If the data looks more categoric then do a more usual hist
    # plot. TODO 080811 Add in a density plot - just need to get the
    # maximum frequency as hs$count above. BUT the density makes no
    # sense, because the bars are the actual data, there is no
    # grouping.

    #if (packageIsAvailable("vcd"))
    #  cols2 <- "col=colorspace::rainbow_hcl(30, start = 270, end = 150)"
    #else
    #  cols2 <- "col=rainbow(30)"

    altplot.cmd <- paste('plot(as.factor(round(ds[ds$grp=="All", 1], ',
                         'digits=2)), col="grey90")\n',
                         #'dens <- density(ds[ds$grp=="All",1], na.rm=TRUE)\n',
                         #'rs<- max(summary(as.factor(round(ds[ds$grp=="All",',
                         #'1], digits=2))))/max(dens$y)\n',
                         #'lines(dens$x, dens$y*rs, type="l")',
                         sep="")

    for (s in seq_len(nhisplots))
    {
      startLog(Rtxt("Plot a Histogram"))
      
      cmd <- paste("sprintf(bind.cmd,",
                   paste(paste('"', rep(hisplots[s], length(targets)+1), '"',
                               sep=""),
                         collapse=","),
                   ")")
      cmd <- eval(parse(text=cmd))
      appendLog(sprintf(Rtxt("Generate just the data for a histogram of",
                             "the variable '%s'."),
                        hisplots[s]),
                paste("ds <-", cmd))
      ds <- eval(parse(text=cmd))

      # 090524 Perform the max y calculation here for each plot
      
      hs <- hist(ds[ds$grp=="All",1], breaks="fd", plot=FALSE)
      dens <- density(ds[ds$grp=="All",1], na.rm=TRUE)
      rs <- max(hs$counts)/max(dens$y)
      maxy <- max(dens$y*rs)
      
      if (length(targets))
        for (i in 1:length(targets))
        {
          dens <- density(ds[ds$grp==targets[i], 1], na.rm=TRUE)
          maxy <- max(c(maxy, dens$y*rs))
        }

      if (newplot && pcnt %% pmax == 0) newPlot(pmax)
      pcnt <- pcnt + 1
      
      # Determine whether to plot a histogram of the numeric data or
      # as a factor (is.integer and unique <= 20).

      dsmin <- eval(parse(text="min(ds[ds$grp=='All',1], na.rm=TRUE)"))
      dsmax <- eval(parse(text="max(ds[ds$grp=='All',1], na.rm=TRUE)"))
      dsuni <- eval(parse(text="unique(ds[ds$grp=='All',1], na.rm=TRUE)"))

      # 080925 Determine the likely number of bars for the plot. This
      # does not always seem to get it correct.
      
      nbars <- nclass.FD(na.omit(ds[ds$grp=="All",1]))

      if (length(dsuni) <= 20 && dsmax - dsmin <= 20)
      {
        appendLog(Rtxt("Plot the data."), altplot.cmd)
        eval(parse(text=altplot.cmd))
      }
      else
      {
        # 090523 REMOVE The nbars is no longer required because we are now
        # using a single colour.
        #
        # plot.cmd <- paste(preplot.cmd, sprintf(plot.cmd, nbars), sep="\n")
        #
        # 090524 REMOVE If I don't worry about trying to determine a maximum,
        # simply don't crop, through using the xpd parameter as set
        # above, then plots look ugly, so do try....
        #

        # 090524 Note the sprintf here, to dynamically specify the max
        # y value for each individual plot.

        appendLog(Rtxt("Plot the data."), sprintf(plot.cmd, hisplots[s], maxy))
        eval(parse(text=sprintf(plot.cmd, hisplots[s], round(maxy))))
        appendLog(Rtxt("Add a rug to the plot to highlight density distribution."),
                  rug.cmd)
        eval(parse(text=rug.cmd))
        if (stratify && length(targets))
        {
          # 090524 REMOVE
#          legend.cmd <- sprintf(paste('legend("topright", c(%s),',
#                                      'fill=c("black", rainbow(%s)))'),
          legend.cmd <- sprintf('legend("topright", c(%s), bty="n", %s)',
                                paste(sprintf('"%s"', c(Rtxt("All"), targets)),
                                      collapse=", "),
                                sprintf(sub("col", "fill", cols),
                                        length(targets)+1))
          appendLog(Rtxt("Add a legend to the plot."), legend.cmd)
          eval(parse(text=legend.cmd))
        }
      }
      
      title.cmd <- genPlotTitleCmd(generateTitleText(hisplots[s], target, sampling,
                                                     stratify && length(targets)))
      appendLog(Rtxt("Add a title to the plot."), title.cmd)
      eval(parse(text=title.cmd))
    }
  }
  }
  
  #---------------------------------------------------------------------
  
  if (not.null(cumplots))
  {
    # Cumulative plot for numeric data.

    nplots <- length(cumplots)

    lib.cmd <- "library(Hmisc, quietly=TRUE)"
    
    for (s in seq_len(nplots))
    {
      startLog()

      if (packageIsAvailable("colorspace"))
        col <- colorspace::rainbow_hcl(length(targets)+1) #, start = 30, end = 300)
      else
        col <- rainbow(length(targets)+1)
      
      plot.cmd <- paste('Ecdf(ds[ds$grp=="All",1],',
                        sprintf('col="%s",', col[1]),
                        'xlab="%s", lwd=2,',
                        sprintf('ylab=expression(%s <= x),', Rtxt("Proportion")),
                        'subtitles=FALSE)\n')
      if (not.null(targets))
        for (t in seq_along(targets))
        {
          plot.cmd <- paste(plot.cmd,
                            sprintf('Ecdf(ds[ds$grp=="%s",1], ', targets[t]),
                            sprintf('col="%s", lty=%d, ', col[t+1], t+1),
                            'xlab="", lwd=2, subtitles=FALSE, add=TRUE)\n',
                            sep="")
        }

      if (packageIsAvailable("colorspace"))
        cols <- "col=colorspace::rainbow_hcl(%d)" # 090524, start = 30, end = 300)"
      else
        cols <- "col=rainbow(%d)"

      if (not.null(targets))
        legend.cmd <- sprintf(paste('legend("bottomright", c(%s), bty="n", ',
                                   cols, ", lwd=2, lty=1:%d,",
                                   'inset=c(0.05,0.05))'),
                             paste(sprintf('"%s"', c("All", targets)),
                                   collapse=","),
                             length(targets)+1, length(targets)+1)
        
      cmd <- paste("sprintf(bind.cmd,",
                    paste(paste('"', rep(cumplots[s], length(targets)+1), '"',
                               sep=""),
                         collapse=","),
                   ")")
      cmd <- eval(parse(text=cmd))
      appendLog(sprintf(Rtxt("Generate just the data for an Ecdf plot of",
                             "the variable '%s'."),
                        cumplots[s]),
                paste("ds <-", cmd))
       ds <- eval(parse(text=cmd))

      if (newplot && pcnt %% pmax == 0) newPlot(pmax)
      pcnt <- pcnt + 1

      if (! packageIsAvailable("Hmisc", Rtxt("plot cumulative charts"))) break()

      appendLog(packageProvides("Hmisc", "Ecdf"), lib.cmd)
      eval(parse(text=lib.cmd))

      this.plot.cmd <- sprintf(plot.cmd, cumplots[s])
      appendLog(Rtxt("Plot the data."), this.plot.cmd)
      eval(parse(text=this.plot.cmd))
      title.cmd <- genPlotTitleCmd(generateTitleText(cumplots[s], target,
                                                     sampling, length(targets)))
      if (not.null(targets))
      {
        appendLog(Rtxt("Add a legend to the plot."), legend.cmd)
        eval(parse(text=legend.cmd))
      }

      appendLog(Rtxt("Add a title to the plot."), title.cmd)
      eval(parse(text=title.cmd))
    }
  }

  ##---------------------------------------------------------------------

  if (nbenplots > 0)
  {
    # Plot Benford's Law for numeric data.

    if (advanced.graphics &&
        packageIsAvailable("ggplot2", Rtxt("plot using ggplot2")))
    {
      executeBenfordPlot2(dataset, benplots, target, targets, stratify, sampling, pmax)
    }
    else
    {
      
    barbutton <- theWidget("benford_bars_checkbutton")$getActive()
    absbutton <- theWidget("benford_abs_radiobutton")$getActive()
    posbutton <- theWidget("benford_pos_radiobutton")$getActive()
    negbutton <- theWidget("benford_neg_radiobutton")$getActive()
    digspin <- theWidget("benford_digits_spinbutton")$getValue()
    numspin <- theWidget("benford_num_digits_spinbutton")$getValue()

    if (numspin > 1)
    {
      errorDialog(Rtxt("Traditional graphics option allows only a single digit.",
                       "Try Advanced Graphics option under Settings."))
      return()
    }
    
    
    benopts <- sprintf(', sp="%s", digit=%d',
                       ifelse(absbutton, "none",
                              ifelse(posbutton, "positive", "negative")),
                       digspin)
    
    # Using barplot2 from gplots
    
    lib.cmd <- "library(gplots, quietly=TRUE)"

    if (packageIsAvailable("colorspace"))
      cols <- "colorspace::rainbow_hcl(%d)" # 090524, start = 30, end = 300)"
    else
      cols <- "rainbow(%d)"

    # Calculate the expected distribution according to Benford's Law

    if (digspin == 1)
      if (numspin == 1)
        expect.cmd <- paste('unlist(lapply(1:9, function(x) log10(1 + 1/x)))')
      else
        expect.cmd <- paste('unlist(lapply(10:99, function(x) log10(1 + 1/x)))')
    # see http://www.mathpages.com/home/kmath302/kmath302.htm
    else if (digspin > 1) 
      expect.cmd <- sprintf(paste('unlist(lapply(0:9, function(x) {sum(log10',
                                  '(1 + 1/(10*(seq(10^(%d-2), ',
                                  '(10^(%d-1))-1)) + x)))}))'),
                            digspin, digspin)

    # Construct the command to plot the distribution.

    xlab.txt <- switch(as.character(digspin),
                       "1"=Rtxt("Distribution of the First Digit"),
                       "2"=Rtxt("Distribution of the Second Digit"),
                       "3"=Rtxt("Distribution of the Third Digit"),
                       sprintf(Rtxt("Distribution of the %dth Digit"), digspin))

    if (barbutton)
    {
      plot.cmd <- paste('barplot2(ds, beside=TRUE, ',
                        'col=c("black", ', sprintf(cols, length(targets)+1), "), ",
                        'xlab="', xlab.txt, '", ylab="', Rtxt("Probability"), '")', sep="")
    }
    else
    {
      plot.cmd <- paste('plot(', ifelse(digspin==1, "1", "0"),
                        ':9, ds[1,], type="b", pch=19, col=',
                        '"black"', ', ',
# 090524                        sprintf(cols, 1), ', ',
                       'ylim=c(0,max(ds)), axes=FALSE, ',
                       'xlab="', xlab.txt,
                        '", ylab="', Rtxt("Probability"), '")\n',
                        'axis(1, at=',
                        ifelse(digspin==1, "1", "0"),
                        ':9)\n', 'axis(2)\n',
                       sprintf(paste('points(%d:9, ds[2,],',
                                     'col=%s[1], pch=19, type="b")\n'),
                               ifelse(digspin==1, 1, 0),
                               ifelse(is.null(target),
                                      sprintf(cols, 1),
                                      sprintf(cols, length(targets)+1))),
                       sep="")

      if (not.null(targets))
        for (i in seq_along(targets))
        {
          plot.cmd <- sprintf(paste('%s\npoints(%d:9, ds[%d,],',
                                   'col=%s[%d], pch=%d, type="b")'),
                             plot.cmd, ifelse(digspin==1, 1, 0), i+2,
                             sprintf(cols, length(targets)+1),
                             i+1, 19)
        }
    }
    if (packageIsAvailable("gplots", Rtxt("plot a bar chart for Benford's Law")))
    {
      startLog("Benford's Law")
      
      appendLog(packageProvides("gplots", "barplot2"), lib.cmd)
      eval(parse(text=lib.cmd))
      
      appendLog(Rtxt("Generate the expected distribution for Benford's Law."),
               paste("expect <-", expect.cmd))
      expect <- eval(parse(text=expect.cmd))

      if (is.null(targets) && ! barbutton)
      {
        # Plot all Benford's plots on the one line graph

        bc <- sub(Rtxt("All"), "%s", substr(bind.cmd, 7, nchar(bind.cmd)-1))
        new.bind.cmd <- substr(bind.cmd, 1, 6)
        data.cmd <- 't(as.matrix(data.frame(expect=expect'
        plot.cmd <- paste('plot(1:9, ds[1,], type="b", ',
                         'pch=19, col=',
                          sprintf(cols, 1),
                          ', ',
                         'ylim=c(0,max(ds)), axes=FALSE, ',
                         'xlab="', Rtxt("Initial Digit"), '", ylab="',
                          Rtxt("Probability"), '")\n',
                         'axis(1, at=1:9)\n', 'axis(2)\n',
                         sep="")
        for (s in seq_len(nbenplots))
        {
          new.bind.cmd <- paste(new.bind.cmd, 
                           sprintf(bc, benplots[s], benplots[s]),
                           ",\n     ",
                           sep="")
          data.cmd <- paste(data.cmd, ",\n     ",
                           sprintf(paste('"%s"=calcInitialDigitDistr',
                                         '(ds[ds$grp=="%s", 1]%s)', sep=""),
                                   benplots[s], benplots[s], benopts),
                           sep="")
          plot.cmd <- paste(plot.cmd,
                           sprintf(paste('points(1:9, ds[%d,],',
                                         'col=%s[%d], pch=19, type="b")\n'),
                                   s+1, sprintf(cols, nbenplots+1), s+1),
                           sep="")
        }
        new.bind.cmd <- paste(substr(new.bind.cmd, 1,
                                     nchar(new.bind.cmd)-7), ")",
                            sep="")
        data.cmd <- paste(data.cmd, ")))", sep="")

        legend.cmd <- sprintf(paste('legend("%s", c(%s), bty="n", ',
                                   'fill=%s)'),
                              ifelse(digspin>2, "botright", "topright"),
                              paste(sprintf('"%s"',
                                            c("Benford", benplots)),
                                    collapse=","),
                              sprintf(cols, nbenplots+1))

        appendLog(Rtxt("Generate the required data."),
                 paste("ds <-", new.bind.cmd))
        ds <- eval(parse(text=new.bind.cmd))

        appendLog(Rtxt("Generate the data specifically for the plot."),
                  paste("ds <-", data.cmd))
        ds <- eval(parse(text=data.cmd))

        if (newplot && pcnt %% pmax == 0) newPlot(pmax)
        pcnt <- pcnt + 1

        par(xpd=TRUE)
        
        appendLog(Rtxt("Display the plot."), plot.cmd)
        eval(parse(text=plot.cmd))

        appendLog(Rtxt("Add a legend to the plot."), legend.cmd)
        eval(parse(text=legend.cmd))

        if (sampling)
          if (length(targets))
            title.txt <- sprintf(Rtxt("Benford's Law (sample)\nby %s"), target)
          else
            title.txt <- Rtxt("Benford's Law (sample)")
        else
          if (length(targets))
            title.txt <- sprintf(Rtxt("Benford's Law\nby %s"), target)
          else
            title.txt <- Rtxt("Benford's Law")
          
        title.cmd <- genPlotTitleCmd(title.txt)

        appendLog(Rtxt("Add a title to the plot."), title.cmd)
        eval(parse(text=title.cmd))
      }
      else
      {
        # Plot multiple graphs since we have a target, and will split
        # each graph according to the target values.
        
        for (s in seq_len(nbenplots))
        {
          #
          # Record the sizes of the subsets for the legend
          #
          sizes.cmd <- paste('sizes <<- (function(x)(paste(names(x), " (",',
                             ' x, ")", sep="")))(by(ds, ds$grp, nrow))')
          
          data.cmd <- paste('t(as.matrix(data.frame(expect=expect,\n    ',
                           'All=calcInitialDigitDistr(ds[ds$grp=="All", 1]',
                            benopts, ')')
        
          if (not.null(targets))
            for (t in seq_along(targets))
              data.cmd <- paste(data.cmd, ",\n     ",
                               sprintf('"%s"=', targets[t]),
                               'calcInitialDigitDistr(ds[ds$grp==',
                               sprintf('"%s", ', targets[t]), '1]',
                                benopts, ')',
                               sep="")
          data.cmd <- paste(data.cmd, ")))", sep="")

          if (not.null(targets))
            if (barbutton)
              legend.cmd <- sprintf(paste('legend("topright", c(%s), bty="n", ',
# 090524                                         'fill=heat.colors(%d), title="%s")'),
                                         'fill=c("black", %s))'),
                                   paste(sprintf('"%s"',
                                                 c(Rtxt("Benford"), Rtxt("All"), targets)),
                                         collapse=","),
                                   sprintf(cols, length(targets)+1))
            else
              legend.cmd <- sprintf(paste('legend("%s", c(%s), inset=.05, bty="n",',
                                          'fill=c("black", %s))'),
                                    ifelse(digspin>2, "bottomright",
                                           "topright"),
                                    paste('"', Rtxt("Benford's"), '", sizes', sep=""),
                                    sprintf(cols, length(targets)+1))
          else
            if (barbutton)
              legend.cmd <- paste('legend("topright", c("', Rtxt("Benford"),
                                  '", "', Rtxt("All"), '"), bty="n", ',
                                 'fill=heat.colors(2))', sep="")
            else
              legend.cmd <- paste('legend("topright", c("', Rtxt("Benford"),
                                  '", "', Rtxt("All"), '"), bty="n", ',
                                  'fill=', sprintf(cols, 2), sep="")
          
          cmd <- paste("sprintf(bind.cmd,",
                       paste(paste('"', rep(benplots[s], length(targets)+1),
                                   '"', sep=""), collapse=","),
                       ")")
          cmd <- eval(parse(text=cmd))
          
          appendLog(sprintf(Rtxt("Generate the data for the plot of the variable '%s'."),
                            benplots[s]),
                   paste("ds <-", cmd))
          ds <- eval(parse(text=cmd))

          appendLog(Rtxt("Generate legend entries with subset sizes."),
                    gsub("<<-", "<-", sizes.cmd))
          eval(parse(text=sizes.cmd))
          
          appendLog(Rtxt("Generate the frequency of the initial digits."),
                   paste("ds <-", data.cmd))
          ds <- eval(parse(text=data.cmd))

          nan.cmd <- "ds[is.nan(ds)] <- 0"
          appendLog(Rtxt("Ensure rows with no digits are treated as zeros."), nan.cmd)
          ds[is.nan(ds)] <- 0
          
          if (newplot && pcnt %% pmax == 0) newPlot(pmax)
          pcnt <- pcnt + 1

          par(xpd=TRUE)
          
          appendLog(Rtxt("Display the plot."), plot.cmd)
          eval(parse(text=plot.cmd))
          
          appendLog(Rtxt("Add a legend to the plot."), legend.cmd)
          eval(parse(text=legend.cmd))

          # 100312 This is messy because we need the full strings for
          # translations.
          
          if (sampling)
            if (posbutton)
              if (length(targets))
                title.txt <- sprintf(Rtxt("Benford's Law: %s (sample) (positive values)\nby %s"),
                                     benplots[s], target)
              else
                title.txt <- sprintf(Rtxt("Benford's Law: %s (sample) (positive values)"),
                                     benplots[s])
            else if (negbutton)
              if (length(targets))
                title.txt <- sprintf(Rtxt("Benford's Law: %s (sample) (negative values)\nby %s"),
                                     benplots[s], target)
              else
                title.txt <- sprintf(Rtxt("Benford's Law: %s (sample) (negative values)"),
                                     benplots[s])
            else
              if (length(targets))
                title.txt <- sprintf(Rtxt("Benford's Law: %s (sample)\nby %s"),
                                     benplots[s], target)
              else
                title.txt <- sprintf(Rtxt("Benford's Law: %s (sample)"),
                                     benplots[s])
          else
            if (posbutton)
              if (length(targets))
                title.txt <- sprintf(Rtxt("Benford's Law: %s (positive values)\nby %s"),
                                     benplots[s], target)
              else
                title.txt <- sprintf(Rtxt("Benford's Law: %s (positive values)"),
                                     benplots[s])
            else if (negbutton)
              if (length(targets))
                title.txt <- sprintf(Rtxt("Benford's Law: %s (negative values)\nby %s"),
                                     benplots[s], target)
              else
                title.txt <- sprintf(Rtxt("Benford's Law: %s (negative values)"),
                                     benplots[s])
            else
              if (length(targets))
                title.txt <- sprintf(Rtxt("Benford's Law: %s\nby %s"),
                                     benplots[s], target)
              else
                title.txt <- sprintf(Rtxt("Benford's Law: %s"),
                                     benplots[s])
            
          
          title.cmd <- genPlotTitleCmd(title.txt)
          
          appendLog(Rtxt("Add a title to the plot."), title.cmd)
          eval(parse(text=title.cmd))
        }
      }
    }
  }
  }
  
  

  #---------------------------------------------------------------------

  if (nbarplots > 0)
  {
    # Plot a frequency plot for a categoric variable.

    # Use barplot2 from gplots.
    
    lib.cmd <- "library(gplots, quietly=TRUE)"

    # Construct a generic data command built using the genericDataSet
    # values. To generate a barplot we use the output of the summary
    # command on each element in the genericDataSet, and bring them
    # together into a single structure. The resulting generic.data.cmd
    # will have a number of "%s"s (one for the whole dataset, then one
    # for each level) from the original genericDataSet string that
    # will be replaced with the name of each variable as it is being
    # plotted.

    generic.data.cmd <- paste(lapply(genericDataSet,
                                   function(x) sprintf("summary(na.omit(%s))", x)),
                            collapse=",\n    ")
    generic.data.cmd <- sprintf("rbind(%s)", generic.data.cmd)

    # If the gplots package is available then generate a plot for each
    # chosen vairable.
    
    if (packageIsAvailable("gplots", Rtxt("plot a bar chart")))
    {
      startLog()
      appendLog(packageProvides("gplots", "barplot2"), lib.cmd)
      eval(parse(text=lib.cmd))

      for (s in seq_len(nbarplots))
      {
        startLog(Rtxt("Bar Plot"))

        # Construct and evaluate a command string to generate the
        # data for the plot.

        ds.cmd <- paste(sprintf("sprintf('%s',", generic.data.cmd),
                        paste(paste('"', rep(barplots[s], length(targets)+1),
                                    '"', sep=""), collapse=","), ")")
        ds.cmd <- eval(parse(text=ds.cmd))

        appendLog(Rtxt("Generate the summary data for plotting."),
                 paste("ds <-", ds.cmd))
        ds <- eval(parse(text=ds.cmd))
        
        ## Construct and evaluate the command to plot the
        ## distribution.  Determine maxium value so that the y axis
        ## can extend to it. We save the output from barplot2 in order
        ## to add numbers to the plot.
    
        if (newplot && pcnt %% pmax == 0) newPlot(pmax)
        pcnt <- pcnt + 1

        #if (is.null(target))
        #  ord.cmd <- 'order(ds[1,])'
        #else
          ord.cmd <- 'order(ds[1,], decreasing=TRUE)'
        appendLog(Rtxt("Sort the entries."), paste("ord <-", ord.cmd))
        ord <- eval(parse(text=ord.cmd))

        cols <- sprintf(ifelse(packageIsAvailable("colorspace"),
                               "colorspace::rainbow_hcl(%s)", # 090524, start = 270, end = 150)",
                               "rainbow(%s)"),
                       length(targets)+1) 
        
        maxFreq <- max(ds)
        plot.cmd <- sprintf(paste('barplot2(ds[,ord], beside=TRUE, ',
                                  'ylab="', Rtxt("Frequency"), '", xlab="%s", ',
                                  'ylim=c(0, %d), col=%s)', sep=""),
                            barplots[s], round(maxFreq+maxFreq*0.20), cols)
        appendLog(Rtxt("Plot the data."), paste("bp <- ", plot.cmd))
        bp <- eval(parse(text=plot.cmd))

        ## Construct and evaluate a command to add text to the top of
        ## the bars in the bar chart. Only do this if there are not
        ## too many values for the category, otherwise the numbers
        ## look bad. I could, alternatively, scale the font?

        if (ncol(bp) <= 5)
        {
          text.cmd <- sprintf("text(bp, ds[,ord]+%d, ds[,ord])",
                             round(maxFreq*0.040))
          appendLog(Rtxt("Add the actual frequencies."), text.cmd)
          eval(parse(text=text.cmd))
        }

        ## Construct and evaluate a command to add a legend to the
        ## plot, but only if there is a target, optherwise it is
        ## obvious.
        
        if (not.null(targets))
        {
          legend.cmd <- sprintf(paste('legend("topright", bty="n", c(%s), ',
                                     "fill=%s)"),
                               paste(sprintf('"%s"', c(Rtxt("All"), targets)),
                                     collapse=","),
                               cols)
          appendLog(Rtxt("Add a legend to the plot."), legend.cmd)
          eval(parse(text=legend.cmd))
        }
        
        # Construct and evaluate a command to add the title to the
        # plot.
        
        title.cmd <- genPlotTitleCmd(generateTitleText(barplots[s], target,
                                                       sampling, length(targets)))

        appendLog(Rtxt("Add a title to the plot."), title.cmd)
        eval(parse(text=title.cmd))
      }
    }
  }
  
########################################################################
#
#   Pairs plot
#
########################################################################

  if (npaiplots > 0 )
  {
    if (packageIsAvailable("GGally","Display a Pairs plot with ggpairs"))
    {
      executePairsPlotSelect2(dataset, paiplots, target, targets, stratify, sample, pmax)
    }
  }

  if (ndotplots > 0)
  {
    
    # Construct a generic data command built using the genericDataSet
    # values. To generate a barplot we use the output of the summary
    # command on each element in the genericDataSet, and bring them
    # together into a single structure. The resulting generic.data.cmd
    # will have a number of "%s"s (one for the whole dataset, then one
    # for each level) from the original genericDataSet string that
    # will be replaced with the name of each variable as it is being
    # plotted.

    generic.data.cmd <- paste(lapply(genericDataSet,
                                   function(x) sprintf("summary(na.omit(%s))", x)),
                            collapse=",\n    ")
    generic.data.cmd <- sprintf("rbind(%s)", generic.data.cmd)

    # This should have been removed at some stage! We seem to be using
    # dotchart from grpahics now.
    #
    #    appendLog("Use dotplot from lattice for the plots.", lib.cmd)
    #    eval(parse(text=lib.cmd))

    for (s in seq_len(ndotplots))
    {

      startLog(Rtxt("Dot Plot"))

      # Construct and evaluate a command string to generate the data
      # for the plot.

      ds.cmd <- paste(sprintf("sprintf('%s',", generic.data.cmd),
                     paste(paste('"', rep(dotplots[s], length(targets)+1),
                                 '"', sep=""), collapse=","), ")")
      ds.cmd <- eval(parse(text=ds.cmd))
      appendLog(Rtxt("Generate the summary data for the plot."),
               paste("ds <-", ds.cmd))
      ds <- eval(parse(text=ds.cmd))

      # Construct and evaluate the command to determine the order in
      # which to print the catgories, from larges to smallest.

      if (is.null(target))
        ord.cmd <- 'order(ds[1,])'
      else
        ord.cmd <- 'order(ds[1,], decreasing=TRUE)'
      appendLog(Rtxt("Sort the entries."),
               paste("ord <-", ord.cmd))
      ord <- eval(parse(text=ord.cmd))
        
      # Construct and evaluate the command to plot the distribution.
    
      if (newplot && pcnt %% pmax == 0) newPlot(pmax)
      pcnt <- pcnt + 1
      
      titles <- genPlotTitleCmd(generateTitleText(dotplots[s], target,
                                                  sampling, length(targets)),
                                vector=TRUE)

      cols <- sprintf(ifelse(packageIsAvailable("colorspace"),
                             "colorspace::rainbow_hcl(%s)", # 090524, start = 270, end = 150)",
                             "rainbow(%s)"),
                      length(targets)+1) 

      plot.cmd <- sprintf(paste('dotchart(%s, main="%s", sub="%s", ',
                                'col=rev(%s),%s ',
                                'xlab="', Rtxt("Frequency"), '", ylab="%s", ',
                                'pch=c(%s, 19))',
                                sep=""),
                          # 090525 reverse the row order to get the
                          # order I want in the dot chart - start with
                          # All, and then the rest. It is not clear
                          # wht dotplots does this.
                          "ds[nrow(ds):1,ord]", titles[1], titles[2], cols,
                          ifelse(is.null(target), "", ' labels="",'), dotplots[s],
                          sprintf("1:%s", length(targets)))
      appendLog(Rtxt("Plot the data."), plot.cmd)
      eval(parse(text=plot.cmd))

      if (not.null(target))
      {
        legend.cmd <- sprintf(paste('legend("bottomright", bty="n",',
                                   'c(%s), col=%s,',
                                   'pch=c(19, %s))'),
                              paste(sprintf('"%s"', c(Rtxt("All"), targets)),
                                    collapse=","),
                              cols,
                              sprintf("%s:1", length(targets)))
        appendLog(Rtxt("Add a legend."), legend.cmd)
        eval(parse(text=legend.cmd))
      }
    }
  }

  #---------------------------------------------------------------------

  for (s in seq_len(nmosplots))
  {

    startLog(Rtxt("Mosaic Plot"))

    # Construct and evaluate a command string to generate the
    # data for the plot.

    if (is.null(target))
      ds.cmd <- sprintf("table(crs$dataset%s$%s)",
                        ifelse(sampling, "[crs$sample,]", ""),
                        mosplots[s])
    else
      ds.cmd <- paste(sprintf(paste("table(crs$dataset%s$%s,",
                                      "crs$dataset%s$%s)"),
                              ifelse(sampling, "[crs$sample,]", ""), mosplots[s],
                              ifelse(sampling, "[crs$sample,]", ""), target))
    appendLog(Rtxt("Generate the table data for plotting."),
              paste("ds <-", ds.cmd))
    ds <- eval(parse(text=ds.cmd))

    # Construct and evaluate the command to determin the order in
    # which to print the catgories, from larges to smallest.

      if (is.null(target))
        ord.cmd <- 'order(ds, decreasing=TRUE)'
      else
        ord.cmd <- "order(apply(ds, 1, sum), decreasing=TRUE)"
    appendLog(Rtxt("Sort the entries."), paste("ord <-", ord.cmd))
    ord <- eval(parse(text=ord.cmd))
    
    # Construct and evaluate the command to plot the
    # distribution.
    
    if (newplot && pcnt %% pmax == 0) newPlot(pmax)
    pcnt <- pcnt + 1

    if (sampling)
      if (length(target))
        title.txt <- sprintf(Rtxt("Mosaic of %s (sample)\nby %s"), mosplots[s], target)
      else
        title.txt <- sprintf(Rtxt("Mosaic of %s (sample)"), mosplots[s])
    else
      if (length(target))
        title.txt <- sprintf(Rtxt("Mosaic of %s\nby %s"), mosplots[s], target)
      else
        title.txt <- sprintf(Rtxt("Mosaic of %s"), mosplots[s])

    titles <- genPlotTitleCmd(title.txt, vector=TRUE)

    if (packageIsAvailable("colorspace"))
      cols <- "color=colorspace::rainbow_hcl(%d)" # 090524, start = 270, end = 150)"
    else
      cols <- "color=rainbow(%d)"

    plot.cmd <- sprintf(paste('mosaicplot(ds[ord%s], main="%s", sub="%s", ',
                              cols, '%s, cex=0.7, xlab="%s", ylab="%s")',
                              sep=""),
                        ifelse(is.null(target), "", ","),
                        titles[1], titles[2], length(targets)+1,
                        ifelse(is.null(target), "", "[-1]"),
                        mosplots[s], target)
    appendLog(Rtxt("Plot the data."), plot.cmd)
    eval(parse(text=plot.cmd))
  }
  
  # Update the status bar.
  
  if (total.plots > 1)
    setStatusBar(sprintf(Rtxt("All %d plots have been generated."), total.plots))
  else if (total.plots ==  1)
    setStatusBar(Rtxt("One plot has been generated."))
  else
    setStatusBar(Rtxt("A pairs plot has been generated."))
}

# 120205 The following is migrating into the original
# executeExplorePlot, and splitting out the different plots into
# separate functions.

executeExplorePlot2 <- function(dataset,
                               boxplots = getSelectedVariables("boxplot"),
                               hisplots = getSelectedVariables("hisplot"),
                               cumplots = getSelectedVariables("cumplot"),
                               benplots = getSelectedVariables("benplot"),
                               barplots = getSelectedVariables("barplot"),
                               dotplots = getSelectedVariables("dotplot"),
                               mosplots = getSelectedVariables("mosplot"),
                               paiplots = getSelectedVariables("paiplot"),
                               stratify=TRUE, sampling=NULL,
                               target=crs$target)
{
  # 091214 Plot the data using ggplot2.
  #
  # 091214 Box plots have been convereted to ggplot2
  #
  # The DATASET is a character string that defines the dataset to
  # use. Information about what variables to plot and the kind of
  # plots is obtained from the continuous_treeview and the
  # categorical_treeview which are displayed in the Explore tab's
  # Distribution option. The appropriate plots are displayed. 090323
  # By having the list of varaiables to display for each type of plot
  # set in the parameter list we can call this function to plot
  # variables from the command line, as in using Sweave. The function
  # remains an internal Rattle function though - only for those who
  # know!

  cat("Experimental ggplot2 Version\n")

  # Obtain the selection of variables.

  nboxplots <- length(boxplots)
  nhisplots <- length(hisplots)
  nbenplots <- length(benplots)
  nbarplots <- length(barplots)
  ndotplots <- length(dotplots)
  nmosplots <- length(mosplots)
  npaiplots <- as.integer(length(paiplots) >= 2) # only one plot is generated, 2 variables minimum are needed

  if (npaiplots == 0 && length(paiplots)>0)
  {
    infoDialog(Rtxt("The Pairs plot requires at least two variables."))
  }
  
  total.plots <- nboxplots + nhisplots + length(cumplots) +
    nbenplots + nbarplots + ndotplots + nmosplots + npaiplots
  
  pmax <- theWidget("plots_per_page_spinbutton")$getValue()
  pcnt <- 0

  # Iterate over all target values if a target is defined and has
  # less than 10 values. The plots will then also display the
  # distributions per target value.

  # 091011 Move to using the value of crs$target instead of getting it
  # from the interface. Evenetually, pass target in as an argument so
  # we can be independent of the GUI?

  # target <- getSelectedVariables("target")

  if (length(target))
    targets <- levels(as.factor(crs$dataset[[target]]))
  else
    targets <- NULL

  if (length(targets) > 10) targets <- NULL

  # For now, let's plot always, since I was wondering why the Benford
  # plot was not showing all the targets!
  
##   if (length(targets) > 10)
##   {
##     target <- NULL
##     targets <- NULL
##   }
  
  # 091011 Check if there are mosaic plots requested but there is not
  # target - notify that the mosaic plots will not be plotted.

  if (nmosplots > 0 && is.null(target))
  {
    infoDialog(Rtxt("A mosaic plot can not be displayed without identifying a target",
                    "variable. The requested mosaic plots will be ignored."))
    total.plots <- total.plots - nmosplots
    mosplots <- NULL
    nmosplots <- 0
  }

  # Don't waste real estate if we are plotting less than number
  # allowed per page.
  
  if (total.plots < pmax) pmax <- total.plots
  
  # Check for sampling.

  if (is.null(sampling))
  {
    use.sample <- theWidget("data_sample_checkbutton")$getActive()
    sampling  <- use.sample && not.null(crs$sample)
  }

  # Record other options.

  annotate <- theWidget("explot_annotate_checkbutton")$getActive()
  
  # Split the data, first for all values.

  bind.cmd <- sprintf('rbind(data.frame(%%s=%s[,"%%s"], %s="All")', dataset, target)

  for (i in seq_along(targets))
  {
    bind.cmd <- sprintf("%s,\n            data.frame(%%s=%s",
                        bind.cmd, dataset)
    
    bind.cmd <- sprintf('%s[crs$dataset%s$%s=="%s","%%s"], %s="%s")',
                        bind.cmd,
                        ifelse(sampling, "[crs$sample,]", ""),
                        target, targets[i], target, targets[i])
  }
  
  # Finish off the command to create the dataset for plotting.
  
  bind.cmd <- sprintf("%s)", bind.cmd)

  # Build a list of generic datasets. This describes how to get the
  # relevant rows from the dataset for All the data, then each of the
  # levels of a target. Each contains a "%s" which is replace gor
  # specific chosen variables at the time of using this construct to
  # obtain the data for the plot. The form is:
  #
  # All = crs$dataset$%s
  #
  # or if sampling is enabled:
  #
  # All = crs$dataset[crs$sample,]$%s  
  #
  # For each level:
  #
  # '0' = crs$dataset[crs$dataset$Adjusted=="0",]$%s
  #
  # or if sampling is enabled:
  #
  # '0' = crs$dataset[crs$sample,][crs$dataset[crs$sample,]$Adjusted=="0",]$%s
  #
  # This is a newer alternative to identifying the dataset
  # segments. We build this list of target and a specification of the
  # correspending data subset. Eventually move all plotting to use
  # this approach rather than using bind.cmd.

  genericDataSet <- data.frame(All=sprintf('%s$%%s', dataset))
  for (i in seq_along(targets))
  {
    tmpDataSet <- data.frame(New=sprintf('%s[crs$dataset%s$%s=="%s",]$%%s',
                               dataset,
                               ifelse(sampling, "[crs$sample,]", ""),
                               target, targets[i]))
    colnames(tmpDataSet) <-  c(targets[i])
    genericDataSet <- cbind(genericDataSet, tmpDataSet)
  }

  # Generate a plot for each variable. If there are too many
  # variables, ask the user if we want to continue.

  if (total.plots > 10 && pmax == 1)
    if (! questionDialog(sprintf(Rtxt("We are about to generate %d",
                                      "individual plots. That's quite a few.",
                                      "You could select fewer variables, or you",
                                      "can change the number of plots per page,",
                                      "but you can also proceed if you like.",
                                      "\n\nWould you like to proceed?"),
                                 total.plots)))
      return()

  #---------------------------------------------------------------------
  # 091005 If no plots are specified then generate a matrix of scatter
  # plots following the example in the pairs function.
  
  if (total.plots == 0)
    displayPairsPlot(dataset)
    

  #---------------------------------------------------------------------

  if (nboxplots > 0)
  {
    # Show a box plot for numeric data. A box plot shows the
    # distribution of numeric data graphically. The box iteself
    # extends from the lower to the upper quartiles with the median
    # drawn in the box. The lines then extend to the maximum and
    # minimum points that are no more than 1.5 times the interquartile
    # range from the median. Outliers are then also plotted as
    # points.

    # 080918 Use the colorspace package to get a better colour map. See
    # http://epub.wu-wien.ac.at/dyn/virlib/wp/eng/showentry?ID=epub-wu-01_c87 and
    # http://statmath.wu.ac.at/~zeileis/papers/Zeileis+Hornik+Murrell-2009.pdf
    
    if (packageIsAvailable("colorspace"))
      cols <- "fill=colorspace::rainbow_hcl(%d)," # 090524, start = 270, end = 150),"
    else
      cols <- "fill=rainbow(%d),"

    # 091214 Generate the compnents of the plot and then bring them
    # togther. The specific variable name ends up being a %s, ready
    # for substitution within the loop over the variables selected for
    # a box plot.
    
    ggplot.cmd  <- sprintf("ggplot(ds, aes(%s, %%s))", target)
    boxplot.cmd <- sprintf("geom_boxplot(%s)",
                           sprintf(cols, length(targets)+1))

    title.txt <- generateTitleText("%%s", target, sampling, stratify && length(targets))
    
    title.cmd <- sprintf(paste('ggtitle("%s")', title.txt))

    sub.cmd <- sprintf('labs(x="%s\\n\\n%%s")', target)

    plot.cmd <- sprintf("print(%s + %s + %s + %s)", ggplot.cmd,
                        boxplot.cmd, title.cmd, sub.cmd)

    ## 091214 I don't see how to add annotations yet. So skip all this for now.
    
    ## # Based on an example from Jim Holtman on r-help 070406.
    
    ## annotate.cmd <- paste("for (i in seq(ncol(bp$stats)))",
    ##                       "{text(i,",
    ##                       "bp$stats[,i] - 0.02*(max(ds$dat, na.rm=TRUE)",
    ##                       "- min(ds$dat, na.rm=TRUE)),",
    ##                       "labels=bp$stats[,i])}")
    
    ## lib.cmd <- "library(doBy, quietly=TRUE)"
    
    ## # TODO: Try using "by" instead of needing another package to
    ## # provide summaryBy. Also, the new version of doBy (061006) seems
    ## # to be outputting extra status information that makes the R
    ## # Console a little chatty unneccessarily - perhaps this will
    ## # disappear again - it looks like debugging information!
    ## #
    ## # status:
    ## # lhsvar     : dat 
    ## # rhsvar     : grp 
    ## # idvar      :  
    ## # fun.names  : mean 
    ## # varPrefix  : mean 
    ## # newNames   : mean.dat 

    ## # Only use summaryBy if there is a target, because it fails if
    ## # there is actually only one group in the data. Might be a new
    ## # bug in the doBy package.
    
    ## if (length(targets) > 1)
    ##   mean.cmd <- paste(sprintf("points(1:%d,", length(targets)+1),
    ##                     "summaryBy(dat ~ grp, data=ds,",
    ##                     "FUN=mean, na.rm=TRUE)$dat.mean,",
    ##                     "pch=8)")
    ## else
    ##   mean.cmd <- paste(sprintf("points(1:%d,", length(targets)+1),
    ##                     "mean(ds$dat, na.rm=TRUE),",
    ##                     "pch=8)")
    
    for (s in seq_len(nboxplots))
    {

      startLog(Rtxt("Box Plot"))
      cmd <- paste("sprintf(bind.cmd,",
                   paste(paste('"', rep(boxplots[s], 2*(length(targets)+1)), '"',
                               sep=""),
                         collapse=", "),
                   ")")
      cmd <- eval(parse(text=cmd))
      appendLog(sprintf(Rtxt("Subset the data for a boxplot of the variable %s",
                             " with the appropriate groups."),
                        boxplots[s]),
                paste("ds <-", cmd))
      ds <- eval(parse(text=cmd))

      # if (pcnt %% pmax == 0) newPlot(pmax)
      # pcnt <- pcnt + 1

      newPlot()

      this.plot <- sprintf(plot.cmd, boxplots[s], boxplots[s],
                           genPlotTitleCmd("", vector=TRUE)[2])
      
      appendLog(Rtxt("Display the plot."), this.plot)
      eval(parse(text=this.plot))

      ## 091214 Annotations are not yet supported.
      
      ## # Add a value for the mean to each boxplot.
      
      ## if (packageIsAvailable("doBy", "add means to box plots"))
      ## {
      ##   appendLog("Use the doBy package to group the data for means.",
      ##            lib.cmd)
      ##   eval(parse(text=lib.cmd))

      ##   appendLog("Calculate the group means.", mean.cmd)
      ##   eval(parse(text=mean.cmd))
      ## }
        
      ## # Optionally include annotations.

      ## if (annotate)
      ## {
      ##   appendLog("Add annotations to the plot.", annotate.cmd)
      ##   eval(parse(text=annotate.cmd))
      ## }        
      
    }
  }
#  }
  

  ##--------------------------------------------------------------------
  
  if (nhisplots > 0)
  {
    # Plot a histogram for numeric data. 090523 Bob Muenchen suggested
    # using a fixed colour for the bars and use the colours as in the
    # box plots for the breakdown between class values, for
    # consistency. I have also introduced colour for the rug and
    # ensured the y limits are calcuated in case the density plot is
    # higher than the histogram. The current colours are not yet right
    # though. For binary data the distinction between the blue and
    # gree is very difficult to see for thin lines.

    #if (packageIsAvailable("RColorBrewer"))
    #  cols <- 'col=brewer.pal(%s, "Set1")'
    #else
    if (packageIsAvailable("colorspace")) # 090524 Why vcd? comes from colorspace....
      # cols <- "col=rainbow_hcl(%s, start = 270, end = 150)"
      cols <- "col=colorspace::rainbow_hcl(%s)"# 090524, start = 0, end = 150)"
    else
      cols <- "col=rainbow(%s)"

    # 090523 The preplot.cmd was considered for use in determining the
    # highest point of the desnity plots so that we can ensure we have
    # enough on the y axis to draw them all. We could use xpd=NA so as
    # not to crop, but then the y axis may not go up high enough and
    # the density may still be higher than fits in the device, and it
    # looks ugly.
    #
    # 090524 An alternative could be, as with much of Rattle, that
    # this cropping is good enough, and the user can fine tune it
    # themselves. But even with xpd=NA the lines look ugly. Also note
    # that xpd=NA makes the left and right extensions of the desnity
    # plots visible, which is not so nice.
    #
    # The original attempt, whereby the below hist (with plot=FALSE)
    # was included in a preplot.cmd, did not work because somehow this
    # plot actually interferred with multiple plots on a page.
    #
    # Note this tip that could help:
    # http://wiki.r-project.org/rwiki/doku.php?id=tips:easier:tbtable
    #

    ## hs <- hist(ds[ds$grp=="All",1], breaks="fd", plot=FALSE)
    ## dens <- density(ds[ds$grp=="All",1], na.rm=TRUE)
    ## rs <- max(hs$counts)/max(dens$y)
    ## maxy <- max(dens$y*rs)
    ## print(maxy)
    
    ## if (length(targets))
    ##   for (i in 1:length(targets))
    ##   {
    ##     dens <- density(ds[ds$grp==targets[i], 1], na.rm=TRUE)
    ##     maxy <- max(c(maxy, dens$y*rs))
    ##   }
    ## print(maxy)

    # 090524 A note on using breaks in the histogram: The default is
    # to use Sturges' rule from 1926. Rob Hyndman
    # (http://robjhyndman.com/papers/sturges.pdf) argues that Sturges'
    # rule is wrong. Sturges' rule leads to oversmoothed
    # histograms. Hyndman notes that Scott's (1979) rule and Freedman
    # and Diaconis's (1981) rule are just as simple to use as Sturges'
    # rule, but are well-founded in statistical theory. Sturges' rule
    # does not work for large n. Also see
    # http://www.math.leidenuniv.nl/~gill/teaching/statistics/histogram.pdf
    # So let's use the Freedman and Diaconis approach.
    
    plot.cmd <- paste('hs <- hist(ds[ds$grp=="All",1], main="", xlab="%s", ',
                      # cols,
                      'col="grey90"',
                      ', ylim=c(0, %s)', #', ceiling(maxy), ')',
                      ', breaks="fd", border=TRUE)\n',
                      'dens <- density(ds[ds$grp=="All",1], na.rm=TRUE)\n',
                      # 090523 Now done in preplot.cmd - not any more? 091021
                      'rs <- max(hs$counts)/max(dens$y)\n',
                      'lines(dens$x, dens$y*rs, type="l", ',
                      sprintf(cols, length(targets)+1), '[1])',
                      sep="")
    if (stratify && length(targets))
    {
      plot.cmd <- paste(plot.cmd, "\n",
                        paste(sprintf(paste('dens <- density(ds[ds$grp=="%s",',
                                            '1], na.rm=TRUE)\n',
                                            'rs <- max(hs$counts)/max(dens$y)\n',
                                            'lines(dens$x, dens$y*rs, ',
                                            'type="l", ',
                                            '%s[%s])', sep=""),
                                      targets,
                                      sprintf(cols, length(targets)+1),
                                      (seq_along(targets)+1)),
                                      #eval(parse(text=sprintf(cols,
                                      #             length(targets))))),
                              collapse="\n"),
                        sep="")
    }

    if (stratify && length(targets))
    {
      rug.cmd <- paste(sprintf('rug(ds[ds$grp=="%s", 1], %s[%s])', targets,
                               sprintf(cols, length(targets)+1),
                               seq_along(targets)+1), collapse="\n")
    }
    else
    {
      rug.cmd <- 'rug(ds[ds$grp=="All",1])'
    }

    # If the data looks more categoric then do a more usual hist
    # plot. TODO 080811 Add in a density plot - just need to get the
    # maximum frequency as hs$count above. BUT the density makes no
    # sense, because the bars are the actual data, there is no
    # grouping.

    #if (packageIsAvailable("vcd"))
    #  cols2 <- "col=colorspace::rainbow_hcl(30, start = 270, end = 150)"
    #else
    #  cols2 <- "col=rainbow(30)"

    altplot.cmd <- paste('plot(as.factor(round(ds[ds$grp=="All", 1], ',
                         'digits=2)), col="grey90")\n',
                         #'dens <- density(ds[ds$grp=="All",1], na.rm=TRUE)\n',
                         #'rs<- max(summary(as.factor(round(ds[ds$grp=="All",',
                         #'1], digits=2))))/max(dens$y)\n',
                         #'lines(dens$x, dens$y*rs, type="l")',
                         sep="")

    for (s in seq_len(nhisplots))
    {
      startLog(Rtxt("Plot a Histogram"))
      
      cmd <- paste("sprintf(bind.cmd,",
                   paste(paste('"', rep(hisplots[s], length(targets)+1), '"',
                               sep=""),
                         collapse=","),
                   ")")
      cmd <- eval(parse(text=cmd))
      appendLog(sprintf(Rtxt("Generate just the data for a histogram of",
                             "the variable '%s'."),
                        hisplots[s]),
                paste("ds <-", cmd))
      ds <- eval(parse(text=cmd))

      # 090524 Perform the max y calculation here for each plot
      
      hs <- hist(ds[ds$grp=="All",1], breaks="fd", plot=FALSE)
      dens <- density(ds[ds$grp=="All",1], na.rm=TRUE)
      rs <- max(hs$counts)/max(dens$y)
      maxy <- max(dens$y*rs)
      
      if (length(targets))
        for (i in 1:length(targets))
        {
          dens <- density(ds[ds$grp==targets[i], 1], na.rm=TRUE)
          maxy <- max(c(maxy, dens$y*rs))
        }

      if (pcnt %% pmax == 0) newPlot(pmax)
      pcnt <- pcnt + 1
      
      # Determine whether to plot a histogram of the numeric data or
      # as a factor (is.integer and unique <= 20).

      dsmin <- eval(parse(text="min(ds[ds$grp=='All',1], na.rm=TRUE)"))
      dsmax <- eval(parse(text="max(ds[ds$grp=='All',1], na.rm=TRUE)"))
      dsuni <- eval(parse(text="unique(ds[ds$grp=='All',1], na.rm=TRUE)"))

      # 080925 Determine the likely number of bars for the plot. This
      # does not always seem to get it correct.
      
      nbars <- nclass.FD(na.omit(ds[ds$grp=="All",1]))

      if (length(dsuni) <= 20 && dsmax - dsmin <= 20)
      {
        appendLog(Rtxt("Plot the data."), altplot.cmd)
        eval(parse(text=altplot.cmd))
      }
      else
      {
        # 090523 REMOVE The nbars is no longer required because we are now
        # using a single colour.
        #
        # plot.cmd <- paste(preplot.cmd, sprintf(plot.cmd, nbars), sep="\n")
        #
        # 090524 REMOVE If I don't worry about trying to determine a maximum,
        # simply don't crop, through using the xpd parameter as set
        # above, then plots look ugly, so do try....
        #

        # 090524 Note the sprintf here, to dynamically specify the max
        # y value for each individual plot.

        appendLog(Rtxt("Plot the data."), sprintf(plot.cmd, hisplots[s], maxy))
        eval(parse(text=sprintf(plot.cmd, hisplots[s], round(maxy))))
        appendLog(Rtxt("Add a rug to the plot to highlight density distribution."),
                  rug.cmd)
        eval(parse(text=rug.cmd))
        if (stratify && length(targets))
        {
          # 090524 REMOVE
#          legend.cmd <- sprintf(paste('legend("topright", c(%s),',
#                                      'fill=c("black", rainbow(%s)))'),
          legend.cmd <- sprintf('legend("topright", c(%s), bty="n", %s)',
                                paste(sprintf('"%s"', c(Rtxt("All"), targets)),
                                      collapse=", "),
                                sprintf(sub("col", "fill", cols),
                                        length(targets)+1))
          appendLog(Rtxt("Add a legend to the plot."), legend.cmd)
          eval(parse(text=legend.cmd))
        }
      }

      if (sampling)
        if (stratify && length(targets))
          title.msg <- sprintf(Rtxt("Distribution of %s (sample)\nby %s"), hisplots[s], target)
        else
          title.msg <- sprintf(Rtxt("Distribution of %s (sample)"), hisplots[s])
      else
        if (stratify && length(targets))
          title.msg <- sprintf(Rtxt("Distribution of %s\nby %s"), hisplots[s], target)
        else
          title.msg <- sprintf(Rtxt("Distribution of %s"), hisplots[s])
      
      title.cmd <- genPlotTitleCmd(title.msg)
      appendLog(Rtxt("Add a title to the plot."), title.cmd)
      eval(parse(text=title.cmd))
    }
  }

  #---------------------------------------------------------------------
  
  if (not.null(cumplots))
  {
    # Cumulative plot for numeric data.

    nplots <- length(cumplots)

    lib.cmd <- "library(Hmisc, quietly=TRUE)"
    
    for (s in seq_len(nplots))
    {
      startLog()

      if (packageIsAvailable("colorspace"))
        col <- colorspace::rainbow_hcl(length(targets)+1) #, start = 30, end = 300)
      else
        col <- rainbow(length(targets)+1)
      
      plot.cmd <- paste('Ecdf(ds[ds$grp=="All",1],',
                        sprintf('col="%s",', col[1]),
                        'xlab="%s",',
                        'ylab=expression(Proportion <= x),',
                        'subtitles=FALSE)\n')
      if (not.null(targets))
        for (t in seq_along(targets))
        {
          plot.cmd <- paste(plot.cmd,
                            sprintf('Ecdf(ds[ds$grp=="%s",1], ', targets[t]),
                            sprintf('col="%s", lty=%d, ', col[t+1], t+1),
                            'xlab="", subtitles=FALSE, add=TRUE)\n',
                            sep="")
        }

      if (packageIsAvailable("colorspace"))
        cols <- "col=colorspace::rainbow_hcl(%d)" # 090524, start = 30, end = 300)"
      else
        cols <- "col=rainbow(%d)"

      if (not.null(targets))
        legend.cmd <- sprintf(paste('legend("bottomright", c(%s), bty="n", ',
                                   cols, ", lty=1:%d,",
                                   'inset=c(0.05,0.05))'),
                             paste(sprintf('"%s"', c(Rtxt("All"), targets)),
                                   collapse=","),
                             length(targets)+1, length(targets)+1)
        
      cmd <- paste("sprintf(bind.cmd,",
                    paste(paste('"', rep(cumplots[s], length(targets)+1), '"',
                               sep=""),
                         collapse=","),
                   ")")
      cmd <- eval(parse(text=cmd))
      appendLog(sprintf(Rtxt("Generate just the data for an Ecdf plot of",
                             "the variable '%s'."),
                        cumplots[s]),
                paste("ds <-", cmd))
       ds <- eval(parse(text=cmd))

      if (pcnt %% pmax == 0) newPlot(pmax)
      pcnt <- pcnt + 1

      if (! packageIsAvailable("Hmisc", Rtxt("plot cumulative charts"))) break()

      appendLog(packageProvides("Hmisc", "Ecdf"), lib.cmd)
      eval(parse(text=lib.cmd))

      this.plot.cmd <- sprintf(plot.cmd, cumplots[s])
      appendLog(Rtxt("Plot the data."), this.plot.cmd)
      eval(parse(text=this.plot.cmd))

      if (sampling)
        if (stratify && length(targets))
          title.msg <- sprintf(Rtxt("Cummulative %s (sample)\nby %s"), cumplots[s], target)
        else
          title.msg <- sprintf(Rtxt("Cummulative %s (sample)"), cumplots[s])
      else
        if (stratify && length(targets))
          title.msg <- sprintf(Rtxt("Cummulative %s\nby %s"), cumplots[s], target)
        else
          title.msg <- sprintf(Rtxt("Cummulative %s"), cumplots[s])
      
      title.cmd <- genPlotTitleCmd(title.msg)

      if (not.null(targets))
      {
        appendLog(Rtxt("Add a legend to the plot."), legend.cmd)
        eval(parse(text=legend.cmd))
      }

      appendLog(Rtxt("Add a title to the plot."), title.cmd)
      eval(parse(text=title.cmd))
    }
  }

  ##---------------------------------------------------------------------

  if (nbenplots > 0)
  {
    ## Plot Benford's Law for numeric data.

    barbutton <- theWidget("benford_bars_checkbutton")$getActive()
    absbutton <- theWidget("benford_abs_radiobutton")$getActive()
    posbutton <- theWidget("benford_pos_radiobutton")$getActive()
    negbutton <- theWidget("benford_neg_radiobutton")$getActive()
    digspin <- theWidget("benford_digits_spinbutton")$getValue()

    benopts <- sprintf(', sp="%s", digit=%d',
                       ifelse(absbutton, "none",
                              ifelse(posbutton, "positive", "negative")),
                       digspin)
    
    # Using barplot2 from gplots
    
    lib.cmd <- "library(gplots, quietly=TRUE)"

    if (packageIsAvailable("colorspace"))
      cols <- "colorspace::rainbow_hcl(%d)" # 090524, start = 30, end = 300)"
    else
      cols <- "rainbow(%d)"

    # Calculate the expected distribution according to Benford's Law

    if (digspin == 1)
      expect.cmd <- paste('unlist(lapply(1:9, function(x) log10(1 + 1/x)))')
    # see http://www.mathpages.com/home/kmath302/kmath302.htm
    else if (digspin > 1) 
      expect.cmd <- sprintf(paste('unlist(lapply(0:9, function(x) {sum(log10',
                                  '(1 + 1/(10*(seq(10^(%d-2), ',
                                  '(10^(%d-1))-1)) + x)))}))'),
                            digspin, digspin)

    # Construct the command to plot the distribution.

    xlab.txt <- switch(as.character(digspin),
                       "1"=Rtxt("Distribution of the First Digit"),
                       "2"=Rtxt("Distribution of the Second Digit"),
                       "3"=Rtxt("Distribution of the Third Digit"),
                       sprintf(Rtxt("Distribution of the %dth Digit"), digspin))

    if (barbutton)
    {
      plot.cmd <- paste('barplot2(ds, beside=TRUE, ',
                        'col=c("black", ', sprintf(cols, length(targets)+1), "), ",
                        'xlab="', xlab.txt, '", ylab="', Rtxt("Probability"), '")', sep="")
    }
    else
    {
      plot.cmd <- paste('plot(', ifelse(digspin==1, "1", "0"),
                        ':9, ds[1,], type="b", pch=19, col=',
                        '"black"', ', ',
# 090524                        sprintf(cols, 1), ', ',
                       'ylim=c(0,max(ds)), axes=FALSE, ',
                       'xlab="', xlab.txt,
                        '", ylab="', Rtxt("Probability"), '")\n',
                        'axis(1, at=',
                        ifelse(digspin==1, "1", "0"),
                        ':9)\n', 'axis(2)\n',
                       sprintf(paste('points(%d:9, ds[2,],',
                                     'col=%s[1], pch=19, type="b")\n'),
                               ifelse(digspin==1, 1, 0),
                               ifelse(is.null(target),
                                      sprintf(cols, 1),
                                      sprintf(cols, length(targets)+1))),
                       sep="")
      if (not.null(targets))
        for (i in seq_along(targets))
        {
          plot.cmd <- sprintf(paste('%s\npoints(%d:9, ds[%d,],',
                                   'col=%s[%d], pch=%d, type="b")'),
                             plot.cmd, ifelse(digspin==1, 1, 0), i+2,
                             sprintf(cols, length(targets)+1),
                             i+1, 19)
        }
    }
    if (packageIsAvailable("gplots", Rtxt("plot a bar chart for Benford's Law")))
    {
      startLog(Rtxt("Benford's Law"))
      
      appendLog(packageProvides("gplots", "barplot2"), lib.cmd)
      eval(parse(text=lib.cmd))
      
      appendLog(Rtxt("Generate the expected distribution for Benford's Law."),
               paste("expect <-", expect.cmd))
      expect <- eval(parse(text=expect.cmd))

      if (is.null(targets) && ! barbutton)
      {
        # Plot all Benford's plots on the one line graph

        startLog()

        bc <- sub(Rtxt("All"), "%s", substr(bind.cmd, 7, nchar(bind.cmd)-1))
        new.bind.cmd <- substr(bind.cmd, 1, 6)
        data.cmd <- 't(as.matrix(data.frame(expect=expect'
        plot.cmd <- paste('plot(1:9, ds[1,], type="b", ',
                         'pch=19, col=',
                          sprintf(cols, 1),
                          ', ',
                         'ylim=c(0,max(ds)), axes=FALSE, ',
                         'xlab="Initial Digit", ylab="Probability")\n',
                         'axis(1, at=1:9)\n', 'axis(2)\n',
                         sep="")
        for (s in seq_len(nbenplots))
        {
          new.bind.cmd <- paste(new.bind.cmd, 
                           sprintf(bc, benplots[s], benplots[s]),
                           ",\n     ",
                           sep="")
          data.cmd <- paste(data.cmd, ",\n     ",
                           sprintf(paste('"%s"=calcInitialDigitDistr',
                                         '(ds[ds$grp=="%s", 1]%s)', sep=""),
                                   benplots[s], benplots[s], benopts),
                           sep="")
          plot.cmd <- paste(plot.cmd,
                           sprintf(paste('points(1:9, ds[%d,],',
                                         'col=%s[%d], pch=19, type="b")\n'),
                                   s+1, sprintf(cols, nbenplots+1), s+1),
                           sep="")
        }
        new.bind.cmd <- paste(substr(new.bind.cmd, 1,
                                     nchar(new.bind.cmd)-7), ")",
                            sep="")
        data.cmd <- paste(data.cmd, ")))", sep="")

        legend.cmd <- sprintf(paste('legend("%s", c(%s), bty="n", ',
                                   'fill=%s)'),
                              ifelse(digspin>2, "botright", "topright"),
                              paste(sprintf('"%s"',
                                            c("Benford", benplots)),
                                    collapse=","),
                              sprintf(cols, nbenplots+1))

        appendLog(Rtxt("Generate the required data."),
                 paste("ds <-", new.bind.cmd))
        ds <- eval(parse(text=new.bind.cmd))

        appendLog(Rtxt("Generate the data specifically for the plot."),
                  paste("ds <-", data.cmd))
        ds <- eval(parse(text=data.cmd))

        if (pcnt %% pmax == 0) newPlot(pmax)
        pcnt <- pcnt + 1

        par(xpd=TRUE)
        
        appendLog(Rtxt("Display the plot."), plot.cmd)
        eval(parse(text=plot.cmd))

        appendLog(Rtxt("Add a legend to the plot."), legend.cmd)
        eval(parse(text=legend.cmd))
        
        if (sampling)
          title.cmd <- genPlotTitleCmd(sprintf("Benford's Law (sample)%s",
                                               ifelse(length(targets),
                                                  paste("\nby", target), "")))
        else
          title.cmd <- genPlotTitleCmd(sprintf("Benford's Law%s",
                                               ifelse(length(targets),
                                                      paste("\nby", target), "")))

        appendLog(Rtxt("Add a title to the plot."), title.cmd)
        eval(parse(text=title.cmd))

      }
      else
      {
        # Plot multiple graphs since we have a target, and will split
        # each graph according to the target values.
        
        for (s in seq_len(nbenplots))
        {
          startLog()
          #
          # Record the sizes of the subsets for the legend
          #
          sizes.cmd <- paste('sizes <<- (function(x)(paste(names(x), " (",',
                             ' x, ")", sep="")))(by(ds, ds$grp, nrow))')
          
          data.cmd <- paste('t(as.matrix(data.frame(expect=expect,\n    ',
                           'All=calcInitialDigitDistr(ds[ds$grp=="All", 1]',
                            benopts, ')')
        
          if (not.null(targets))
            for (t in seq_along(targets))
              data.cmd <- paste(data.cmd, ",\n     ",
                               sprintf('"%s"=', targets[t]),
                               'calcInitialDigitDistr(ds[ds$grp==',
                               sprintf('"%s", ', targets[t]), '1]',
                                benopts, ')',
                               sep="")
          data.cmd <- paste(data.cmd, ")))", sep="")

          if (not.null(targets))
            if (barbutton)
              legend.cmd <- sprintf(paste('legend("topright", c(%s), bty="n", ',
# 090524                                         'fill=heat.colors(%d), title="%s")'),
                                         'fill=c("black", %s))'),
                                   paste(sprintf('"%s"',
                                                 c(Rtxt("Benford"), Rtxt("All"), targets)),
                                         collapse=","),
                                   sprintf(cols, length(targets)+1))
            else
              legend.cmd <- sprintf(paste('legend("%s", c(%s), inset=.05, bty="n",',
                                         'fill=c("black", %s))'),
                                    ifelse(digspin>2, "bottomright",
                                           "topright"),
                                    paste('"', Rtxt("Benford's"), '", sizes', sep=""),
                                   sprintf(cols, length(targets)+1))
          else
            if (barbutton)
              legend.cmd <- paste('legend("topright", c("Benford", "All"), bty="n",',
                                 'fill=heat.colors(2))')
            else
              legend.cmd <- paste('legend("topright", c("Benford", "All"), bty="n", ',
                                  'fill=', sprintf(cols, 2), sep="")
          
          cmd <- paste("sprintf(bind.cmd,",
                       paste(paste('"', rep(benplots[s], length(targets)+1),
                                   '"', sep=""), collapse=","),
                       ")")
          cmd <- eval(parse(text=cmd))
          
          appendLog(sprintf(Rtxt("Generate the data for the plot of the variable '%s'."),
                            benplots[s]),
                   paste("ds <-", cmd))
          ds <- eval(parse(text=cmd))

          appendLog(Rtxt("Generate legend entries with subset sizes."),
                    gsub("<<-", "<-", sizes.cmd))
          eval(parse(text=sizes.cmd))
          
          appendLog(Rtxt("Generate the frequency of the initial digits."),
                   paste("ds <-", data.cmd))
          ds <- eval(parse(text=data.cmd))

          nan.cmd <- "ds[is.nan(ds)] <- 0"
          appendLog(Rtxt("Ensure rows with no digits are treated as zeros."), nan.cmd)
          ds[is.nan(ds)] <- 0
          
          if (pcnt %% pmax == 0) newPlot(pmax)
          pcnt <- pcnt + 1

          par(xpd=TRUE)
          
          appendLog(Rtxt("Display the plot."), plot.cmd)
          eval(parse(text=plot.cmd))
          
          appendLog(Rtxt("Add a legend to the plot."), legend.cmd)
          eval(parse(text=legend.cmd))
          
          if (sampling)
            title.cmd <- genPlotTitleCmd(sprintf(paste("Benford's Law:",
                                                       "%s (sample)%s%s"),
                                                 benplots[s],
                ifelse(posbutton, " (positive values)",
                       ifelse(negbutton, " (negative values)", "")),
                                                 ifelse(length(targets),
                                                        paste("\nby", target), "")))
          else
            title.cmd <- genPlotTitleCmd(sprintf("Benford's Law: %s%s%s",
                                                 benplots[s],
                ifelse(posbutton, " (positive values)",
                       ifelse(negbutton, " (negative values)", "")),
                                                 ifelse(length(targets),
                                                        paste("\nby", target), "")))
          appendLog(Rtxt("Add a title to the plot."), title.cmd)
          eval(parse(text=title.cmd))
        }
      }
    }
  }

  #---------------------------------------------------------------------

  if (nbarplots > 0)
  {
    # Plot a frequency plot for a categoric variable.

    # Use barplot2 from gplots.
    
    lib.cmd <- "library(gplots, quietly=TRUE)"

    # Construct a generic data command built using the genericDataSet
    # values. To generate a barplot we use the output of the summary
    # command on each element in the genericDataSet, and bring them
    # together into a single structure. The resulting generic.data.cmd
    # will have a number of "%s"s (one for the whole dataset, then one
    # for each level) from the original genericDataSet string that
    # will be replaced with the name of each variable as it is being
    # plotted.

    generic.data.cmd <- paste(lapply(genericDataSet,
                                   function(x) sprintf("summary(na.omit(%s))", x)),
                            collapse=",\n    ")
    generic.data.cmd <- sprintf("rbind(%s)", generic.data.cmd)

    # If the gplots package is available then generate a plot for each
    # chosen vairable.
    
    if (packageIsAvailable("gplots", Rtxt("plot a bar chart")))
    {
      startLog()
      appendLog(packageProvides("gplots", "barplot2"), lib.cmd)
      eval(parse(text=lib.cmd))

      for (s in seq_len(nbarplots))
      {
        startLog(Rtxt("Bar Plot"))

        # Construct and evaluate a command string to generate the
        # data for the plot.

        ds.cmd <- paste(sprintf("sprintf('%s',", generic.data.cmd),
                        paste(paste('"', rep(barplots[s], length(targets)+1),
                                    '"', sep=""), collapse=","), ")")
        ds.cmd <- eval(parse(text=ds.cmd))

        appendLog(Rtxt("Generate the summary data for plotting."),
                 paste("ds <-", ds.cmd))
        ds <- eval(parse(text=ds.cmd))
        
        ## Construct and evaluate the command to plot the
        ## distribution.  Determine maxium value so that the y axis
        ## can extend to it. We save the output from barplot2 in order
        ## to add numbers to the plot.
    
        if (pcnt %% pmax == 0) newPlot(pmax)
        pcnt <- pcnt + 1

        #if (is.null(target))
        #  ord.cmd <- 'order(ds[1,])'
        #else
          ord.cmd <- 'order(ds[1,], decreasing=TRUE)'
        appendLog(Rtxt("Sort the entries."), paste("ord <-", ord.cmd))
        ord <- eval(parse(text=ord.cmd))

        cols <- sprintf(ifelse(packageIsAvailable("colorspace"),
                               "colorspace::rainbow_hcl(%s)", # 090524, start = 270, end = 150)",
                               "rainbow(%s)"),
                       length(targets)+1) 
        
        maxFreq <- max(ds)
        plot.cmd <- sprintf(paste('barplot2(ds[,ord], beside=TRUE,',
                                  'ylab="', Rtxt("Frequency"), '", xlab="%s",',
                                  'ylim=c(0, %d), col=%s)', sep=""),
                            barplots[s], round(maxFreq+maxFreq*0.20), cols)
        appendLog(Rtxt("Plot the data."), paste("bp <- ", plot.cmd))
        bp <- eval(parse(text=plot.cmd))

        ## Construct and evaluate a command to add text to the top of
        ## the bars in the bar chart. Only do this if there are not
        ## too many values for the category, otherwise the numbers
        ## look bad. I could, alternatively, scale the font?

        if (ncol(bp) <= 5)
        {
          text.cmd <- sprintf("text(bp, ds[,ord]+%d, ds[,ord])",
                             round(maxFreq*0.040))
          appendLog(Rtxt("Add the actual frequencies."), text.cmd)
          eval(parse(text=text.cmd))
        }

        ## Construct and evaluate a command to add a legend to the
        ## plot, but only if there is a target, optherwise it is
        ## obvious.
        
        if (not.null(targets))
        {
          legend.cmd <- sprintf(paste('legend("topright", bty="n", c(%s), ',
                                     "fill=%s)"),
                               paste(sprintf('"%s"', c(Rtxt("All"), targets)),
                                     collapse=","),
                               cols)
          appendLog("Add a legend to the plot.", legend.cmd)
          eval(parse(text=legend.cmd))
        }
        
        ## Construct and evaluate a command to add the title to the
        ## plot.
        
        title.cmd <- genPlotTitleCmd(sprintf("Distribution of %s%s%s",
                                             barplots[s],
                                             ifelse(sampling," (sample)",""),
                                             ifelse(length(targets),
                                                    paste("\nby", target), "")))
        appendLog(Rtxt("Add a title to the plot."), title.cmd)
        eval(parse(text=title.cmd))
      }
    }
  }

### REMOVE 080925 - Until work out multiple plots on one device issue.
###   if (nbarplots > 0)
###   {
###     # Plot a frequency plot for a categoric variable.

###     # 080817 Use barchart from lattice instead of barplot2 from
###     # ggplots.

###     lib.cmd <- "library(lattice, quietly=TRUE)"
    
###     # Construct a generic data command built using the genericDataSet
###     # values. To generate a barplot we use the output of the summary
###     # command on each element in the genericDataSet, and bring them
###     # together into a single structure. The resulting generic.data.cmd
###     # will have a number of "%s"s (one for the whole dataset, then one
###     # for each level) from the original genericDataSet string that
###     # will be replaced with the name of each variable as it is being
###     # plotted.

###     generic.data.cmd <- paste(lapply(genericDataSet,
###                                      function(x) sprintf("summary(%s)", x)),
###                               collapse=",\n    ")
###     generic.data.cmd <- sprintf("cbind(%s)", generic.data.cmd)

###     # If the lattice package is available then generate a plot for
###     # each chosen vairable.
    
###     if (packageIsAvailable("lattice", "display a bar chart"))
###     {
###       startLog()
###       appendLog("Load lattice for the barchart function.", lib.cmd)
###       eval(parse(text=lib.cmd))

###       for (s in 1:nbarplots)
###       {
###         startLog()

###         # Construct and evaluate a command string to generate the
###         # data for the plot.

###         ds.cmd <- paste(sprintf("sprintf('%s',", generic.data.cmd),
###                        paste(paste('"', rep(barplots[s], length(targets)+1),
###                                  '"', sep=""), collapse=","), ")")
###         ds.cmd <- eval(parse(text=ds.cmd))
###         appendLog(sprintf("Generate the summary data for plotting %s.", barplots[s]),
###                  paste("ds <-", ds.cmd))
###         ds <- eval(parse(text=ds.cmd))

###         names.cmd <- sprintf('colnames(ds) <- c(%s)',
###                              ifelse(length(targets)==0, '"Frequency"',
###                                     paste('"Frequency"',
###                                           paste(sprintf('"%s"', targets),
###                                                 collapse=", "),
###                                           sep=", ")))
###         appendLog("Set the appropriate column names.", names.cmd)
###         eval(parse(text=names.cmd))

###         # We don't have multiple plots on the one plot implemented yet
###         # - should we? I would guess there is a simple way to do this
###         # with lattice.
        
###         #if (pcnt %% pmax == 0) newPlot(pmax)
###         #pcnt <- pcnt + 1
###         newPlot(pmax)

###         # Construct and evaluate the command to determine the order in
###         # which to print the catgories, from smallest (at the bottom)
###         # to largest.

###         ord.cmd <- 'order(ds[,1])'
###         appendLog("Sort the entries.", paste("ord <-", ord.cmd))
###         ord <- eval(parse(text=ord.cmd))

###         plot.cmd <- sprintf(paste('print(barchart(ds[ord,%s]',
###                                   'xlab="Frequency"',
###                                   ifelse(length(targets)==0,
###                                          'groups=NULL', # Just to have something!
###                                          sprintf(paste('auto.key=list(title="%s",',
###                                                        'cex=0.75,', 'columns=%d)'),
###                                                  target, 2)),
###                                   sprintf('sub="%s"', genPlotTitleCmd(vector=TRUE)),
###                                   'main="Distribution of %s%s"))', sep=", "),
###                             ifelse(length(targets)==0, "", "-1"),
###                             barplots[s],
###                             ifelse(sampling," (sample)",""))
                            
###         appendLog("Plot the data.", plot.cmd)
###         eval(parse(text=plot.cmd))

###       }
###     }
###   }

  ##---------------------------------------------------------------------

### REMOVE 080925 - Until work out multiple plots on one device issue.
###   if (ndotplots > 0)
###   {
    
###     # 080817 Use dotplot(lattice) instead of dotchart. 080925 But not
###     # yet since it uses a different mechanism to get multiple plots on
###     # one device and I've not set that up yet.

###     # lib.cmd <- "library(lattice, quietly=TRUE)"

###     # Construct a generic data command built using the genericDataSet
###     # values. To generate a barplot we use the output of the summary
###     # command on each element in the genericDataSet, and bring them
###     # together into a single structure. The resulting generic.data.cmd
###     # will have a number of "%s"s (one for the whole dataset, then
###     # one for each level) from the original genericDataSet string
###     # that will be replaced with the name of each variable as it is
###     # being plotted.

###     generic.data.cmd <- paste(lapply(genericDataSet,
###                                    function(x) sprintf("summary(%s)", x)),
###                             collapse=",\n    ")
###     generic.data.cmd <- sprintf("cbind(%s)", generic.data.cmd)

###     # If the lattice package is available then generate a plot for
###     # each chosen vairable.
    
###     if (packageIsAvailable("lattice", "display a dot plot"))
###     {
###       startLog()
###       appendLog("Load lattice for the dotplot function.", lib.cmd)
###       eval(parse(text=lib.cmd))

###       for (s in 1:ndotplots)
###       {
###         startLog()

###         # Construct and evaluate a command string to generate the data
###         # for the plot.

###         ds.cmd <- paste(sprintf("sprintf('%s',", generic.data.cmd),
###                         paste(paste('"', rep(dotplots[s], length(targets)+1),
###                                     '"', sep=""), collapse=","), ")")
###         ds.cmd <- eval(parse(text=ds.cmd))
###         appendLog(sprintf("Generate the summary data for plotting %s.", dotplots[s]),
###                   paste("ds <-", ds.cmd))
###         ds <- eval(parse(text=ds.cmd))

###         names.cmd <- sprintf('colnames(ds) <- c(%s)',
###                              ifelse(length(targets)==0, '"Frequency"',
###                                     paste('"Frequency"',
###                                           paste(sprintf('"%s"', targets),
###                                                 collapse=", "),
###                                           sep=", ")))
###         appendLog("Set the appropriate column names.", names.cmd)
###         eval(parse(text=names.cmd))

###         # Construct and evaluate the command to determine the order in
###         # which to print the catgories, from smallest (at the bottom)
###         # to largest.

###         ord.cmd <- 'order(ds[,1])'
###         appendLog("Sort the entries.", paste("ord <-", ord.cmd))
###         ord <- eval(parse(text=ord.cmd))

###         # Construct and evaluate the command to plot the distribution.
    
###         #if (pcnt %% pmax == 0) newPlot(pmax)
###         #pcnt <- pcnt + 1
###         newPlot(pmax)
      
###         plot.cmd <- sprintf(paste('print(dotplot(ds[ord,%s]',
###                                   'xlab="Frequency"',
###                                   'type=c("p", "h", "a")',
###                                   ifelse(length(targets)==0,
###                                          'groups=NULL', # Just to have something!
###                                          sprintf(paste('auto.key=list(title="%s",',
###                                                        'cex=0.75,', 'columns=%d)'),
###                                                  target, 2)),
###                                   sprintf('sub="%s"', genPlotTitleCmd(vector=TRUE)),
###                                   'main="Distribution of %s%s"))', sep=", "),
###                             ifelse(length(targets)==0, "", "-1"),
###                             dotplots[s],
###                             ifelse(sampling," (sample)",""))
###         appendLog("Plot the data.", plot.cmd)
###         eval(parse(text=plot.cmd))
###       }
###     }
###   }

  if (ndotplots > 0)
  {
    
    # Construct a generic data command built using the genericDataSet
    # values. To generate a barplot we use the output of the summary
    # command on each element in the genericDataSet, and bring them
    # together into a single structure. The resulting generic.data.cmd
    # will have a number of "%s"s (one for the whole dataset, then one
    # for each level) from the original genericDataSet string that
    # will be replaced with the name of each variable as it is being
    # plotted.

    generic.data.cmd <- paste(lapply(genericDataSet,
                                   function(x) sprintf("summary(na.omit(%s))", x)),
                            collapse=",\n    ")
    generic.data.cmd <- sprintf("rbind(%s)", generic.data.cmd)

    # This should have been removed at some stage! We seem to be using
    # dotchart from grpahics now.
    #
    #    appendLog("Use dotplot from lattice for the plots.", lib.cmd)
    #    eval(parse(text=lib.cmd))

    for (s in seq_len(ndotplots))
    {

      startLog(Rtxt("Dot Plot"))

      # Construct and evaluate a command string to generate the data
      # for the plot.

      ds.cmd <- paste(sprintf("sprintf('%s',", generic.data.cmd),
                     paste(paste('"', rep(dotplots[s], length(targets)+1),
                                 '"', sep=""), collapse=","), ")")
      ds.cmd <- eval(parse(text=ds.cmd))
      appendLog(Rtxt("Generate the summary data for the plot."),
               paste("ds <-", ds.cmd))
      ds <- eval(parse(text=ds.cmd))

      # Construct and evaluate the command to determine the order in
      # which to print the catgories, from larges to smallest.

      if (is.null(target))
        ord.cmd <- 'order(ds[1,])'
      else
        ord.cmd <- 'order(ds[1,], decreasing=TRUE)'
      appendLog(Rtxt("Sort the entries."),
               paste("ord <-", ord.cmd))
      ord <- eval(parse(text=ord.cmd))
        
      # Construct and evaluate the command to plot the distribution.
    
      if (pcnt %% pmax == 0) newPlot(pmax)
      pcnt <- pcnt + 1
      
      titles <- genPlotTitleCmd(sprintf("Distribution of %s%s%s",
                                        dotplots[s],
                                        ifelse(sampling," (sample)",""),
                                         ifelse(length(targets),
                                                  paste("\nby", target), "")),
                                vector=TRUE)

      cols <- sprintf(ifelse(packageIsAvailable("colorspace"),
                             "colorspace::rainbow_hcl(%s)", # 090524, start = 270, end = 150)",
                             "rainbow(%s)"),
                      length(targets)+1) 

      plot.cmd <- sprintf(paste('dotchart(%s, main="%s", sub="%s",',
                               'col=rev(%s),%s',
                               'xlab="', Rtxt("Frequency"), '", ylab="%s", pch=19)',
                                sep=""),
                          # 090525 reverse the row order to get the
                          # order I want in the dot chart - start with
                          # All, and then the rest. It is not clear
                          # wht dotplots does this.
                         "ds[nrow(ds):1,ord]", titles[1], titles[2], cols,
                         ifelse(is.null(target), "", ' labels="",'), dotplots[s])
      appendLog("Plot the data.", plot.cmd)
      eval(parse(text=plot.cmd))

      if (not.null(target))
      {
        legend.cmd <- sprintf(paste('legend("bottomright", bty="n",',
                                   'c(%s), col=%s,',
                                   'pch=19)'),
                              paste(sprintf('"%s"', c(Rtxt("All"), targets)),
                                    collapse=","),
                              cols)
        appendLog("Add a legend.", legend.cmd)
        eval(parse(text=legend.cmd))
      }
    }
  }

  #---------------------------------------------------------------------

  for (s in seq_len(nmosplots))
  {

    startLog(Rtxt("Mosaic Plot"))

    # Construct and evaluate a command string to generate the
    # data for the plot.

    if (is.null(target))
      ds.cmd <- sprintf("table(crs$dataset%s$%s)",
                        ifelse(sampling, "[crs$sample,]", ""),
                        mosplots[s])
    else
      ds.cmd <- paste(sprintf(paste("table(crs$dataset%s$%s,",
                                      "crs$dataset%s$%s)"),
                              ifelse(sampling, "[crs$sample,]", ""), mosplots[s],
                              ifelse(sampling, "[crs$sample,]", ""), target))
    appendLog(Rtxt("Generate the table data for plotting."),
              paste("ds <-", ds.cmd))
    ds <- eval(parse(text=ds.cmd))

    # Construct and evaluate the command to determin the order in
    # which to print the catgories, from larges to smallest.

      if (is.null(target))
        ord.cmd <- 'order(ds, decreasing=TRUE)'
      else
        ord.cmd <- "order(apply(ds, 1, sum), decreasing=TRUE)"
    appendLog(Rtxt("Sort the entries."), paste("ord <-", ord.cmd))
    ord <- eval(parse(text=ord.cmd))
    
    # Construct and evaluate the command to plot the
    # distribution.
    
    if (pcnt %% pmax == 0) newPlot(pmax)
    pcnt <- pcnt + 1

    if (is.null(target))
      titles <- genPlotTitleCmd(sprintf("Mosaic of %s",
                                        mosplots[s],
                                        ifelse(sampling," (sample)","")),
                                vector=TRUE)
    else
      titles <- genPlotTitleCmd(sprintf("Mosaic of %s %s\nby %s",
                                        mosplots[s],
                                        ifelse(sampling," (sample)",""),
                                        target),
                                vector=TRUE)

    if (packageIsAvailable("colorspace"))
      cols <- "color=colorspace::rainbow_hcl(%d)" # 090524, start = 270, end = 150)"
    else
      cols <- "color=rainbow(%d)"

    plot.cmd <- sprintf(paste('mosaicplot(ds[ord%s], main="%s", sub="%s", ',
                              cols, '%s, cex=0.7, xlab="%s", ylab="%s")',
                              sep=""),
                        ifelse(is.null(target), "", ","),
                        titles[1], titles[2], length(targets)+1,
                        ifelse(is.null(target), "", "[-1]"),
                        mosplots[s], target)
    appendLog(Rtxt("Plot the data."), plot.cmd)
    eval(parse(text=plot.cmd))
  }
  
  # Update the status bar.
  
  if (total.plots > 1)
    setStatusBar(sprintf(Rtxt("All %d plots have been generated."), total.plots))
  else if (total.plots ==  1)
    setStatusBar(Rtxt("One plot has been generated."))
  else
    setStatusBar(Rtxt("A pairs plot has been generated."))
}

displayPairsPlot <- function(dataset)
{
  # 091005 Given a dataset (a string that evaluates to a data frame)
  # generate a pairs plot with the upper panel being a scatter plot
  # matrix, the diagonal being histograms, and the lower panel showing
  # correlations.
  #
  # If there are more than 6 input/target/risk variables, then
  # randomly sample down to 6 since larger than that is harder to
  # read.
  #
  # 100131 Try out the pairs.panels function from the psych package.

  newPlot(1)
  
  # We use a couple of helper functions from the pairs help page.
  
  pre.cmd <- 'panel.hist <- function(x, ...)
{
   usr <- par("usr"); on.exit(par(usr))
   par(usr = c(usr[1:2], 0, 1.5) )
   h <- hist(x, plot = FALSE)
   breaks <- h$breaks; nB <- length(breaks)
   y <- h$counts; y <- y/max(y)
   rect(breaks[-nB], 0, breaks[-1], y, col="cyan", ...)
}
panel.cor <- function(x, y, digits=2, prefix="", cex.cor, ...)
{
   usr <- par("usr"); on.exit(par(usr))
   par(usr = c(0, 1, 0, 1))
   r <- cor(x, y, use="complete")
   txt <- format(c(r, 0.123456789), digits=digits)[1]
   txt <- paste(prefix, txt, sep="")
   if(missing(cex.cor)) cex.cor <- 0.8/strwidth(txt)
   text(0.5, 0.5, txt, cex = cex.cor * r)
}'

  # Determine the variables to display.

  vars <- c(getVariableIndicies(crs$input),
            getVariableIndicies(crs$target),
            getVariableIndicies(crs$risk))
  if (length(vars) <= 6)
    vars <- ""
  else
  {
    vars <- simplifyNumberList(sort(sample(vars, 6)))
  }

  startLog(Rtxt("Scatter Plot"))

  if (packageIsAvailable("psych"))
  {
    # 100131 Use the pairs.panels function from the psych package
    # which does the panels functions for us.

    lib.cmd <- "library(psych, quietly=TRUE)"
    appendLog(packageProvides("psych", "pairs.panels"), lib.cmd)
    eval(parse(text=lib.cmd))
  
    plot.cmd <- paste(sprintf("pairs.panels(%s[%s],", dataset, vars), "scale=TRUE)")
  }
  else
  {
    # The plot command is the pairs function from the graphics package
    # (and hence always available).

    plot.cmd <- paste(sprintf("pairs(%s[%s],", dataset, vars),
                      "diag.panel=panel.hist,",
                      "upper.panel=panel.cor,",
                      "lower.panel=panel.smooth)")

    appendLog(Rtxt("Define support functions for the plot."), pre.cmd)
    eval(parse(text=pre.cmd))
  }

   appendLog(Rtxt("Display a pairs (scatter) plot. Note random selection of variables",
                  "if there are more than 6."), plot.cmd)
  eval(parse(text=plot.cmd))
}

displayPairsPlot2 <- function(dataset)
{
  # 110409 Given a dataset (a string that evaluates to a data frame)
  # generate a pairs plot with the lower panel being a scatter plot
  # matrix and the upper panel showing correlations.

  newPlot(1)
  
  # Determine the variables to display. 110409 We can only display
  # numeric variables.

  vars <-  intersect(c(getVariableIndicies(crs$input),
                       getVariableIndicies(crs$target),
                       getVariableIndicies(crs$risk)),
                     getVariableIndicies(crs$numeric))

  # If there are more than 6 numeric input/target/risk variables, then
  # randomly sample down to 6 since larger than that is harder to
  # read.

  if (length(vars) <= 6)
    vars <- simplifyNumberList(vars)
  else
    vars <- simplifyNumberList(sort(sample(vars, 6)))

  # Inform the log that we will now produce a scatter plot.
  
  startLog(Rtxt("Scatter or Pairs Plot (Plot a Correlation Matrix)"))

  if (! packageIsAvailable("Deducer", Rtxt("plot correlations"))) return()
  lib.cmd <- "library(Deducer, quietly=TRUE)"
  appendLog(packageProvides("Deducer", "ggcorplot"), lib.cmd)
  eval(parse(text=lib.cmd))

  ds <- sprintf("na.omit(%s[%s])", dataset, vars)
  plot.cmd <- sprintf(paste("ggcorplot(cor.matrix(%s), data=%s,",
                            "var_text_size=2, cor_text_limits=c(2,10))"),
                      ds, ds)

  appendLog(Rtxt("Display a pairs (scatter) plot.",
                 "Note random selection of variables",
                 "if there are more than 6."),
            plot.cmd)

  eval(parse(text=sprintf("print(%s)", plot.cmd)))
}
  
executeExploreGGobi <- function(dataset, name=NULL)
{
  # Based on code from Marco Lo

  # 081128 Obtain info from a collection of radio buttons as to
  # whether to brush the data.  E.g., if hclust then:
  #
  # brush.cmd <- "glyph_colour(gg[1]) <- cutree(crs$hclust, 10)"
  #
  # Note also the embed=TRUE option of the ggobi display
  # function. This allows the display to be embedded within a RGtk
  # window, and hence seemlessly become part of Rattle.  
  
  # Construct the commands.

  lib.cmd <- "library(rggobi, quietly=TRUE)"
  ggobi.cmd <- paste('crs$gg <- ggobi(', dataset, # 140206 Remove <<-
                     ifelse(not.null(name), sprintf(', name="%s"', name), ""),
                     ')')

  # Start logging and executing the R code.
  
  if (! packageIsAvailable("rggobi", Rtxt("explore the data using GGobi"))) return()

  startLog(Rtxt("GGobi Data Exploration"))
  appendLog(packageProvides("rggobi", "rggobi"), lib.cmd)
  eval(parse(text=lib.cmd))
  appendLog(Rtxt("Launch the GGobi data visualization application."),
            gsub("<<-", "<-", ggobi.cmd))
  eval(parse(text=ggobi.cmd))
  
  setStatusBar(Rtxt("GGobi executed."))
}

executeExploreCorrelation <- function(dataset, newplot=TRUE)
{
  TV <- "correlation_textview"

  if (is.null(dataset))
  {
    errorDialog(Rtxt("Correlations are calculated only for numeric data.",
                     "No numeric data was found in the dataset."))
    return()
  }

  # Obtain user interface settings.

  advanced.graphics <- theWidget("use_ggplot2")$getActive()
  ordered <- theWidget("explore_correlation_ordered_checkbutton")$getActive()

  # Map the GUI options to the right names - in case the translation
  # (like Japanese) has translated the options. 100331 The
  # getActiveText here is returning garbage. So tried using getActive
  # to get the index instead, but then found that Encoding is the
  # right approach.
  
  method <- theWidget("explore_correlation_method_combobox")$getActiveText()
  Encoding(method) <- "UTF-8"
  ##method <- theWidget("explore_correlation_method_combobox")$getActive() + 1
  method.orig <- method
  method.opts <- paste("Pearson", "Kendall", "Spearman", sep="\n")
  method.rtxt <- Rtxt (method.opts)
  method.opts <- strsplit(method.opts, "\n")[[1]]
  method.rtxt <- strsplit(method.rtxt, "\n")[[1]]
  method <-tolower(method.opts[which(method == method.rtxt)])
  ##method.orig <- method.rtxt[method]
  ##method <-tolower(method.opts[method])

  # Warn if there are too many variables. An alternative is to offer
  # to just plot the variables with the highest amount of correlation.
  
  nvars <- eval(parse(text=sprintf("ncol(%s)", dataset)))
  if (nvars > crv$max.vars.correlation &&
      ! questionDialog(sprintf(Rtxt("You have requested a Correlation plot.",
                                    "\n\nWith %d variables the plot may take",
                                    "some time to display, and the display will",
                                    "be cramped.",
                                    "Consider identifying up to only %d",
                                    "input variables.\n\n",
                                    "Would you like to continue anyhow?"),
                               nvars, crv$max.vars.correlation)))
    return(FALSE)
  
  # Construct the commands.

  # Deal with showing the missing values plot.
  
  nas <- theWidget("explore_correlation_na_checkbutton")$getActive()
  if (nas)
  {
    naids <- NULL
    naids.cmd <- sprintf('naids <- attr(na.omit(t(%s)), "na.action")\n',
                         dataset)
    eval(parse(text=naids.cmd))
    if (is.null(naids))
    {
      errorDialog(Rtxt("The data contains no missing values, and so no",
                       "missing value correlation plot can be generated."))
      return()
    }
    if (length(naids) == 1)
    {
      errorDialog(Rtxt("The data contains only one variable with missing values,",
                       "and so no missing value correlation plot can be generated."))
      return()
    }
  }

  lib.cmd <- sprintf("library(%s, quietly=TRUE)",
                     ifelse(advanced.graphics, "corrplot", "ellipse"))
  crscor.cmd  <- sprintf('%scrs$cor <- cor(%s, use="pairwise", method="%s")',
                         ifelse(nas, naids.cmd, ""),
                         ifelse(nas,
                                sprintf("is.na(%s[naids])", dataset),
                                dataset),
                         method)
  if (ordered)
    crsord.cmd  <- paste("crs$ord <- order(crs$cor[1,])",
                         "crs$cor <- crs$cor[crs$ord, crs$ord]",
                         sep="\n")
    
  print.cmd   <- "print(crs$cor)"
  if (nas)
  {
    print.cmd <- paste(print.cmd,
                       "\ncat('\\n", Rtxt("Count of missing values:"), "\\n')\n",
                       sprintf("print(apply(is.na(%s[naids]),2,sum))",
                               dataset),
                       "\ncat('\\n", Rtxt("Percent missing values:"), "\\n')\n",
                       sprintf(paste("print(100*apply(is.na(%s[naids]),",
                                     "2,sum)/nrow(%s))"),
                               dataset, dataset),
                       sep="")
    
  }
  MAXCORVARS <- 15
  par.cmd <- ifelse(nrow(crs$cor) > MAXCORVARS, "opar <- par(cex=0.5)\n", "")
  opar.cmd <- ifelse(nrow(crs$cor) > MAXCORVARS, "\npar(opar)", "")

  # 100416 I don't know why I need to do this, but without it I get an
  # error from the sprintf.

  Encoding(method.orig) <- "unknown"
  
  if (nas)
    title.txt <- sprintf(Rtxt("Correlation of Missing Values\\n%s using %s"),
                         crs$dataname, method.orig)
  else
    title.txt <- sprintf(Rtxt("Correlation %s using %s"), crs$dataname, method.orig)
                       
  plot.cmd    <- paste(par.cmd,
                       ifelse(advanced.graphics,
                              "corrplot(crs$cor, mar=c(0,0,1,0))\n",
                              paste("plotcorr(crs$cor, ",
                                    'col=colorRampPalette(c("red", "white", "blue"))(11)',
                                    '[5*crs$cor + 6])\n', sep="")),
                       genPlotTitleCmd(title.txt),
                       opar.cmd,
                       sep="")
  
  # Start logging and executing the R code.

  if (! packageIsAvailable(ifelse(advanced.graphics, "corrplot", "ellipse"),
                           Rtxt("display a correlation plot"))) return()
     
  startLog(Rtxt("Generate a correlation plot for the variables."))
  resetTextview(TV)

  if (advanced.graphics)
    appendLog(packageProvides("corrplot", "corrplot"), lib.cmd)
  else
    appendLog(packageProvides("ellipse", "plotcorr"), lib.cmd)
  eval(parse(text=lib.cmd))

  appendLog(Rtxt("Correlations work for numeric variables only."), crscor.cmd)
  if (ordered) appendLog(Rtxt("Order the correlations by their strength."), crsord.cmd)
  appendLog(Rtxt("Display the actual correlations."), print.cmd)
  appendLog(Rtxt("Graphically display the correlations."), plot.cmd)

  appendTextview(TV,
                 sprintf(ifelse(nas,
                                Rtxt("Missing values correlation summary using",
                                     "the '%s' covariance." ),
                                Rtxt("Correlation summary using the '%s' covariance.")),
                         method.orig),
                 "\n\n",
                 Rtxt("Note that only correlations between numeric variables",
                      "are reported."),
                 "\n\n",
                 collectOutput(paste(crscor.cmd,
                                     if (ordered) crsord.cmd,
                                     print.cmd,
                                     sep="\n")))

  if (newplot) newPlot()
  eval(parse(text=paste(crscor.cmd,
               if (ordered) crsord.cmd,
               plot.cmd,
               sep="\n")))
  
  ## Report completion to the user through the Status Bar.
  
  setStatusBar(Rtxt("Correlation plot and summary generated."))
}

executeExploreHiercor <- function(dataset)
{
  if (is.null(dataset))
  {
    errorDialog(Rtxt("Correlations are calculated only for numeric data.",
                     "\n\nNo numeric variables were found in the dataset",
                     "from amongst those that are not ignored.",
                     "\n\nYou may want to use the transform tab to transform",
                     "your categoric data into numeric data."))
    return()
  }

  # Obtain user interface settings. 100407 Map the translated string
  # from the GUI to the correct English version.

  # method <- tolower(theWidget("explore_correlation_method_combobox")$getActiveText())
  method <- theWidget("explore_correlation_method_combobox")$getActive() + 1
  method.opts <- paste("Pearson", "Kendall", "Spearman", sep="\n")
  method.rtxt <- Rtxt (method.opts)
  method.opts <- strsplit(method.opts, "\n")[[1]]
  method.rtxt <- strsplit(method.rtxt, "\n")[[1]]
  method.orig <- method.rtxt[method]
  method <-tolower(method.opts[method])
  
  # Check that we have sufficient data

  ncols <- eval(parse(text=sprintf("NCOL(%s)", dataset)))
  if ( ncols < 2 )
  {
    errorDialog(Rtxt("The dataset contains less than two numeric variables.",
                     "\n\nCorrelations are calculated only for numeric data.",
                     "\n\nYou may want to select more numeric variables or",
                     "use the transform tab to transform",
                     "your categoric variables into numeric variables."))
    
    return()
  }
    
  # Construct the commands.
  
  cor.cmd    <- sprintf('cc <- cor(%s, use="pairwise", method="%s")', dataset, method)
  hclust.cmd <- 'hc <- hclust(dist(cc), method="average")'
  dend.cmd   <- "dn <- as.dendrogram(hc)"

  # Modification by Ed Cox 080130 to increase margin for long variable names

  fontsize <- .75
  if (ncols>45) {fontsize <- .60}
  
  labwidth <- eval(parse(text = paste('max(unlist(lapply(colnames(',
                           dataset, '), "nchar"))) + 2', sep="")))
  plot.cmd   <- paste('op <- par(mar = c(3, 4, 3, ',
                      round(labwidth/3.5, 2), '))\n',
                      'plot(dn, horiz = TRUE, ',
                      'nodePar = list(col = 3:2, ',
                      'cex = c(2.0, 0.75), pch = 21:22, ',
                      'bg=  c("light blue", "pink"), ',
                      'lab.cex = ', fontsize, ', lab.col = "tomato"), ',
                      'edgePar = list(col = "gray", lwd = 2), ',
                      'xlab="Height"',
                      ')\n',
                      genPlotTitleCmd(paste(Rtxt("Variable Correlation Clusters"),
                                            "\n", sep=""),
                                     crs$dataname, Rtxt("using"), method.orig), '\n',
                      'par(op)\n',
                      sep="")

 # plot.cmd   <- paste('plot(dn, horiz=TRUE, ',
 #                     'nodePar=list(col=3:2, cex=c(2.0, 0.75), pch= 21:22, ',
 #                     'bg= c("light blue", "pink"), ',
 #                     'lab.cex = 0.75, lab.col = "tomato"), ',
 #                     'edgePar=list(col="gray", lwd=2)',
 #                     ')\n',
 #                     genPlotTitleCmd("Variable Correlation Clusters",
 #                                    crs$dataname),
 #                     sep="")

  # Start logging and executing the R code.

  startLog(Rtxt("Hierarchical Variable Correlation"))

  appendLog(Rtxt("Generate the correlations (numerics only)."), cor.cmd)
  eval(parse(text=cor.cmd))

  appendLog(Rtxt("Generate hierarchical cluster of variables."), hclust.cmd)
  eval(parse(text=hclust.cmd))

  appendLog(Rtxt("Generate the dendrogram."), dend.cmd)
  eval(parse(text=dend.cmd))

  appendLog(Rtxt("Now draw the dendrogram."), plot.cmd)
  newPlot()
  eval(parse(text=plot.cmd))

  # Report completion to the user through the Status Bar.
  
  setStatusBar(Rtxt("Hierarchical cluster of correlations plotted."))

}

executeExplorePrcomp <- function(dataset)
{
  TV <- "prcomp_textview"
  
  if (is.null(dataset))
  {
    errorDialog(Rtxt("Principal components are only ipmlemented for numeric data.",
                     "No numeric variables were found in the dataset",
                     "from amongst those that are not ignored."))
    return()
  }

  # User interface options
  
  svd <- theWidget("explore_prcomp_pr_radiobutton")$getActive()
  
  # Construct the commands.
  
  prcomp.cmd  <- sprintf(paste('pc <<- %s(na.omit(%s),',
                               'scale=TRUE, center=TRUE, tol=0)'),
                         ifelse(svd, "prcomp", "princomp"),
                         dataset)
  print.cmd   <- "pc"
  summary.cmd <- "summary(pc)"
  plot.cmd    <- paste('plot(pc, main="")',
                       genPlotTitleCmd(Rtxt("Principal Components Importance"),
                                      crs$dataname),
                       if (svd) # 090328 Add X axis labels for prcomp version.
                       paste("axis(1, at=seq(0.7, ncol(pc$rotation)*1.2, 1.2),",
                             "labels=colnames(pc$rotation), lty=0)"),
                       sep="\n")
  biplot.cmd  <- paste('biplot(pc, main="")',
                       genPlotTitleCmd(Rtxt("Principal Components"),
                                      crs$dataname),
                       sep="\n")

  ## Start logging and executing the R code.

  startLog()
  resetTextview(TV)
  
  appendLog(Rtxt("Principal Components Analysis (on numerics only)."),
          gsub("<<-", "<-", prcomp.cmd))
  eval(parse(text=prcomp.cmd))

  appendTextview(TV, Rtxt("Note that principal components on only the numeric\n",
                          "variables is calculated, and so we can not use this\n",
                          "approach to remove categoric variables from ",
                          "consideration.\n\n",
                          "Any numeric variables with relatively large rotation\n",
                          "values (negative or positive) in any of the first few\n",
                          "components are generally variables that you may wish\n",
                          "to include in the modelling."))

  appendLog(Rtxt("Show the output of the analysis."), print.cmd)
  appendTextview(TV, collectOutput(print.cmd, TRUE))
  
  appendLog(Rtxt("Summarise the importance of the components found."), summary.cmd)
  appendTextview(TV, collectOutput(summary.cmd, TRUE))

  newPlot(1)
  appendLog(Rtxt("Display a plot showing the relative importance of the components."),
          plot.cmd)
  eval(parse(text=plot.cmd))
  
  newPlot(1)
  appendLog(Rtxt("Display a plot showing the two most principal components."),
          biplot.cmd)
  eval(parse(text=biplot.cmd))
  
  ## Report completion to the user through the Status Bar.
  
  setStatusBar(Rtxt("A principal components analysis has been completed."))

}

## Removed - Felix decided not to update as per CRAN request.
##
## executeExplorePlaywith <- function(dataset)
## {
##   if (! packageIsAvailable("latticist", Rtxt("interactively explore data"))) return()
##   if (! packageIsAvailable("playwith", Rtxt("interactively explore data"))) return()

##   startLog(Rtxt("Explore Data"))

##   lib.cmd <- "library(latticist)"
##   appendLog(packageProvides("latticist", "latticist"), lib.cmd)
##   eval(parse(text=lib.cmd))

##   latopts <- ""
##   if (length(crs$target))
##     latopts <- sprintf(', spec=list(groups = "%s")', crs$target)
##   plot.cmd <- sprintf("latticist(%s%s)", dataset, latopts)
##   appendLog(Rtxt("Start up latticist."), plot.cmd)
##   eval(parse(text=plot.cmd))
## }

## executeExplorePlotBuilder <- function()
## {
##   # 8 Mar 2012 Currently don't know how to tell Plot builder the
##   # default dataset to use. Nor how to extract from plot builder the
##   # actual ggplot2 command that is generated - would like to capture
##   # that and place it in the Log.
  
##   if (! packageIsAvailable("Deducer", Rtxt("interactively develop a ggplot2 plot using Plot builder"))) return()

##   startLog(Rtxt("Use the Plot Builder dialog from the Deducer package."))

##   lib.cmd <- "library(Deducer, quietly=TRUE)"
##   appendLog(packageProvides("Deducer", "deducer"), lib.cmd)
##   eval(parse(text=lib.cmd))

##   # This use of ds is a hack - not sure how else to do this as Plot
##   # Builder does not notice crs$dataset, presumably only checking for
##   # data frames in .GlobalEnv

##   # 121210 As of now, remove the assign to global env - it is a bad
##   # idea and against CRAN policy. But then Plot Builder will not
##   # notice the dataset? Perhaps at least check first if 'ds' exists
##   # and if not then assign in the global environment. Ripley suggests
##   # (121210) using environments but not sure he understands the
##   # specific issue here.

##   # 130126 Deducer in
##   # Deducer/javasrc/deducer/data/DataViewerController.java, function
##   # refreshData, uses get.objects('data.frame',
##   # includeInherited=FALSE) to get the datasets to select from. This
##   # only uses the global environment. So perhaps an acceptable
##   # solution is to inform the user and to ask their permission.

##   var.name <- paste("ds", format(Sys.time(), format="%y%m%d%H%M%S"), sep="_")
##   if (!questionDialog(sprintf(Rtxt("To use Plot Builder the Rattle dataset needs to",
##                                    "be available in the global environment (the",
##                                    "user's workspace). To do this Rattle will copy",
##                                    "its internal dataset to the variable '%s'.",
##                                    "This will overwrite any variable of the same",
##                                    "name. The copy of the dataset will be removed",
##                                    "after you exit from Plot Builder.\n\n",
##                                    "Are you okay with Rattle doing this?"),
##                               var.name)))
##     return()
##   # I would really rather not do an assign, but is this the only
##   # option - CRAN do not accept this and I must remove this on
##   # submitting to CRAN.
##   assign.cmd <- "assign(var.name, crs$dataset, envir=.GlobalEnv)"
##   appendLog(Rtxt("Place dataset into Global Environment for PlotBuilder."),
##             assign.cmd)
##   eval(parse(text=assign.cmd))

##   plot.cmd <- 'deducer(cmd="Plot builder")'
##   appendLog(Rtxt("Start up Plot builder dialog."), plot.cmd)
##   eval(parse(text=plot.cmd))
##   eval(parse(text=sprintf("rm(%s, envir=.GlobalEnv)", var.name)))
## }

########################################################################
# CALLBACKS

#-----------------------------------------------------------------------
# DISPLAY APPRORIATE TAB
#
# When an Option radio button is selected, display the appropriate tab
#
# 090328 Move to using sub tabs to house the various widgets, rather
# than turning them on and off!


on_summary_radiobutton_toggled <- function(button)
{
  if (button$getActive()) crv$EXPLORE$setCurrentPage(crv$EXPLORE.SUMMARY.TAB)
  setStatusBar()
}

on_explore_distr_radiobutton_toggled <- function(button)
{
  if (button$getActive()) crv$EXPLORE$setCurrentPage(crv$EXPLORE.PLOT.TAB)
  setStatusBar()
}

on_explore_interactive_radiobutton_toggled <- function(button)
{
  if (button$getActive()) crv$EXPLORE$setCurrentPage(crv$EXPLORE.INTERACTIVE.TAB)
  setStatusBar()
}

on_explore_correlation_radiobutton_toggled <- function(button)
{
  if (button$getActive()) crv$EXPLORE$setCurrentPage(crv$EXPLORE.CORRELATION.TAB)
  setStatusBar()
}

on_prcomp_radiobutton_toggled <- function(button)
{
  if (button$getActive()) crv$EXPLORE$setCurrentPage(crv$EXPLORE.PRCOMP.TAB)
  setStatusBar()
}

#-----------------------------------------------------------------------
# DISTRIBUTIONS

cat_toggled <- function(cell, path.str, model)
{
  ## A categoric variable's radio button has been toggled in the
  ## Explore tab's Distribution option. Handle the choice.

  ## The data passed in is the model used in the treeview.

  RGtk2::checkPtrType(model, "GtkTreeModel")

  ## Extract the column number of the model that has changed.

  column <- cell$getData("column")

  ## Get the current value of the corresponding flag
  
  path <- RGtk2::gtkTreePathNewFromString(path.str) # Current row
  iter <- model$getIter(path)$iter           # Iter for the row
  current <- model$get(iter, column)[[1]]    # Get data from specific column

  ## Always invert
  
  model$set(iter, column, !current)

}

con_toggled <- function(cell, path.str, model)
{
  ## A continuous variable's radio button has been toggled in the
  ## Explore tab's Distribution option. Handle the choice.

  ## The data passed in is the model used in the treeview.

  RGtk2::checkPtrType(model, "GtkTreeModel")

  ## Extract the column number of the model that has changed.

  column <- cell$getData("column")
  
  ## Get the current value of the corresponding flag
  
  path <- RGtk2::gtkTreePathNewFromString(path.str) # Current row
  iter <- model$getIter(path)$iter           # Iter for the row
  current <- model$get(iter, column)[[1]]    # Get data from specific column

  model$set(iter, column, !current)

}

on_categorical_clear_button_clicked <- function(action, window)
{
  ## Ensure categoric all check boxes are unchecked.

  set.cursor("watch")

  ## Only clear selected rows.

  tree.selection <- theWidget("categorical_treeview")$getSelection()

  # Use the data parameter to avoid an RGtk2 bug in 2.12.1, fixed in
  # next release. 071117
  tree.selection$selectedForeach(function(model, path, iter, data)
  {
    columns <- crv$CATEGORICAL[["barplot"]]:crv$CATEGORICAL[["mosplot"]]
    for (c in columns) if (model$get(iter, c)[[1]]) model$set(iter, c, FALSE)
    return(FALSE) # Keep going through all rows
  }, TRUE)

  set.cursor()
}

on_continuous_clear_button_clicked <- function(action, window)
{
  # Ensure all continuous check boxes are unchecked.

  set.cursor("watch")

  # Only clear selected rows.

  tree.selection <- theWidget("continuous_treeview")$getSelection()

  # Use the data parameter to avoid an RGtk2 bug in 2.12.1, fixed in
  # next release. 071117
  tree.selection$selectedForeach(function(model, path, iter, data)
  {
    columns <- crv$CONTINUOUS[["boxplot"]]:crv$CONTINUOUS[["benplot"]]
    for (c in columns) if (model$get(iter, c)[[1]]) model$set(iter, c, FALSE)
    return(FALSE) # Keep going through all rows
  }, TRUE)

  set.cursor()
}

on_summary_find_button_clicked <- function(window)
{
  search.str <- theWidget("summary_find_entry")$getText()
  tv <- theWidget("summary_textview")
  start.iter <- tv$getBuffer()$getStartIter()
  summarySearch(tv, search.str, start.iter)
}

on_summary_next_button_clicked <- function(window)
{
  search.str <- theWidget("summary_find_entry")$getText()
  tv <- theWidget("summary_textview")
  last.search.pos <- tv$getBuffer()$getMark('last.search.pos')
  if (is.null(last.search.pos)) return()
  last.search.iter <- tv$getBuffer()$getIterAtMark(last.search.pos)
  summarySearch(tv, search.str, last.search.iter)
}

on_viewdata_find_button_clicked <- function(window)
{
  ## Need to get the root window of the button, and everything else is
  ## in terms of that.

  root <- window$getRootWindow()
  search.str <- root$getWidget("viewdata_find_entry")$getText()
  print(search.str)
  tv <- crs$viewdataGUI$getWidget("viewdata_textview")
  start.iter <- tv$getBuffer()$getStartIter()
  summarySearch(tv, search.str, start.iter)
}

on_viewdata_next_button_clicked <- function(window)
{
  search.str <- crs$viewdataGUI$getWidget("viewdata_find_entry")$getText()
  tv <- crs$viewdataGUI$getWidget("viewdata_textview")
  last.search.pos <- tv$getBuffer()$getMark('last.search.pos')
  if (is.null(last.search.pos)) return()
  last.search.iter <- tv$getBuffer()$getIterAtMark(last.search.pos)
  summarySearch(tv, search.str, last.search.iter)
}

summarySearch <- function(tv, search.str, start.iter)
{
  found <- start.iter$iter$forwardSearch(search.str, 0)
  tvb <- tv$getBuffer()
  if (found$retval)
  {
    tvb$selectRange(found$match.start, found$match.end)
    last.search.pos <-tvb$createMark('last.search.pos', found$match.end)

    tv$scrollToMark(last.search.pos, 0.2)
    while(RGtk2::gtkEventsPending()) RGtk2::gtkMainIterationDo(blocking=FALSE)

    setStatusBar(sprintf(Rtxt("The string '%s' was found."), search.str))
  }
  else
    setStatusBar(sprintf(Rtxt("The string '%s' was not found."), search.str))
}

#-----------------------------------------------------------------------
# CORRELATION

on_explore_correlation_hier_checkbutton_toggled <- function(button)
{
  if (button$getActive())
  {
    theWidget("explore_correlation_ordered_checkbutton")$setSensitive(FALSE)
    theWidget("explore_correlation_na_checkbutton")$setSensitive(FALSE)
  }
  else
  {
    theWidget("explore_correlation_ordered_checkbutton")$setSensitive(TRUE)
    theWidget("explore_correlation_na_checkbutton")$setSensitive(TRUE)
  }
}  


