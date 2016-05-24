# Gnome R Data Miner: GNOME interface to R for Data Mining
#
# Time-stamp: <2015-05-17 08:58:56 gjw>
#
# Test Tab
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
# Callbacks

# 090205 Now leave all Test radio buttons active, always. The user
# selects one and the Group By option may disappear if they select one
# of the Paired tests. This is more logical than using the Group By
# option to toggle the sensitivity of the Paired tests.

on_test_correlation_radiobutton_toggled <- function(button)
{
  if (button$getActive())
  {
    # 100316 Don't change the status of this button from what the user
    # might have set.
    # theWidget("test_groupby_checkbutton")$setActive(FALSE)
    theWidget("test_groupby_checkbutton")$setSensitive(FALSE)
    theWidget("test_groupby_target_label")$setSensitive(FALSE)
    theWidget("test_vars2_label")$setSensitive(TRUE)
    theWidget("test_vars2_combobox")$setSensitive(TRUE)
  }
  else
  {
    if (length(crs$target))
    {
      theWidget("test_groupby_checkbutton")$setSensitive(TRUE)
      theWidget("test_groupby_target_label")$setSensitive(TRUE)
      if (theWidget("test_groupby_checkbutton")$getActive())
      {
        theWidget("test_vars2_label")$setSensitive(FALSE)
        theWidget("test_vars2_combobox")$setSensitive(FALSE)
      }
    }
  }
}

on_test_wilcoxon_signed_radiobutton_toggled <- on_test_correlation_radiobutton_toggled


on_test_groupby_checkbutton_toggled<- function(button)
{
  # 090205 Only deal with enabling/disabling the var2 box on this bing
  # toggled.
  
  if (button$getActive())
  {
    theWidget("test_vars2_label")$setSensitive(FALSE)
    theWidget("test_vars2_combobox")$setSensitive(FALSE)
  }
  else
  {
    theWidget("test_vars2_label")$setSensitive(TRUE)
    theWidget("test_vars2_combobox")$setSensitive(TRUE)
  }
}

########################################################################
# Functionality

resetTestTab <- function(new.dataset=TRUE)
{
  cbox1 <- theWidget("test_vars1_combobox")
  cbox2 <- theWidget("test_vars2_combobox")

  if (new.dataset)
  {
    # 080921 Set up the list of numeric variables to choose from as
    # sample 1 and, optionally,  sample 2.
    
    vl <- colnames(crs$dataset)[getNumericVariables("indicies")]

    if (not.null(vl))
    {
      cbox1$getModel()$clear()
      lapply(vl, cbox1$appendText)
      cbox2$getModel()$clear()
      lapply(vl, cbox2$appendText)
    }
  }

  # If there is a target variable, and it is binary, then enable Group
  # By Target and set the target label appropriately.

  if (length(crs$target))
  {
    theWidget("test_groupby_checkbutton")$setSensitive(TRUE)
    theWidget("test_groupby_checkbutton")$setActive(TRUE)
    theWidget("test_groupby_target_label")$setSensitive(TRUE)
    theWidget("test_groupby_target_label")$setText(crs$target)
  }
  else
  {
    theWidget("test_groupby_checkbutton")$setSensitive(FALSE)
    theWidget("test_groupby_checkbutton")$setActive(FALSE)
    theWidget("test_groupby_target_label")$setSensitive(FALSE)
  }
}

executeTestTab <- function()
{

  TV <- "test_textview"

  if (noDatasetLoaded()) return(FALSE)
  
  # Obtain interface information.
  
  if (theWidget("test_groupby_checkbutton")$getActive() & !
      (theWidget("test_correlation_radiobutton")$getActive()
       || theWidget("test_wilcoxon_signed_radiobutton")$getActive()))
  {
    v1 <- v2 <- theWidget("test_vars1_combobox")$getActiveText()
    if (is.null(v1))
    {
      errorDialog(Rtxt("Please first choose a variable from which the sample",
                       "will be obtained."))
      return(FALSE)
    }

# 100328 We now ignore the group by when using one of the paired
# tests. So we don't need this message anymore.
#
#    if (theWidget("test_correlation_radiobutton")$getActive()
#        || theWidget("test_wilcoxon_signed_radiobutton")$getActive())
#    {
#      errorDialog(Rtxt ("The Correlation and Wilcoxon Signed Rank tests",
#                       "can only be applied to paired populations.",
#                       "The Group By option is not likely to give paired populations",
#                       "(requiring the same number of observations in each sample).",
#                       "Please de-select Group By and choose",
#                       "two variables that represent observations of the entity at",
#                       "two different times. Alternatively, choose a two-sample",
#                       "(non-paired) test if that is appropriate."))
#      return(FALSE)
#    }
    lvl <- levels(as.factor(crs$dataset[[crs$target]]))
    s1 <- sprintf('crs$dataset[["%s"]] == "%s"', crs$target, lvl[1])
    s2 <- sprintf('crs$dataset[["%s"]] == "%s"', crs$target, lvl[2])
    msg <- sprintf(Rtxt("The two samples being compared come from",
                        "the '%s' variable, grouped by\n",
                        "'%s', with values '%s' and '%s'"),
                   v1, crs$target, lvl[1], lvl[2])
  }
  else
  {
    v1 <- theWidget("test_vars1_combobox")$getActiveText()
    if (is.null(v1))
    {
      errorDialog(Rtxt("Please first choose a variable from which the first sample",
                       "will be obtained."))
      return(FALSE)
    }
    v2 <- theWidget("test_vars2_combobox")$getActiveText()
    if (is.null(v2))
    {
      errorDialog(Rtxt("Please first choose a variable from which the second sample",
                       "will be obtained."))
      return(FALSE)
    }
    s1 <- s2 <- ""
    msg <- sprintf(Rtxt("The two samples being compared are the",
                        "two variables, '%s' and '%s'"), v1, v2)
  }

  # Start the log for this task.

  startLog("Perform Test")

  # Ensure the package is available.

  lib.cmd <- "library(fBasics, quietly=TRUE)"
  if (! packageIsAvailable("fBasics", Rtxt("location T-Test"))) return(FALSE)
  appendLog(Rtxt("Use the fBasics package for statistical tests."), lib.cmd)
  eval(parse(text=lib.cmd))

  resetTextview(TV)

  test <- NULL
  preamble <- NULL
  options <- ""
  
  if (theWidget("test_distr_radiobutton")$getActive())
  {
    test <- "ks2Test"
    preamble <- Rtxt("Kolmogorov-Smirnov Test",
                     "\n\nThe Kolmogorov-Smirnov test is a non-parametric test of the",
                     "\nsimilarity of two distributions.",
                     "The null hypothesis is that the",
                     "\ntwo samples are drawn from the same distribution.",
                     "The two-sided and",
                     "\nthe two one-sided tests are performed.",
                     "\n\nThe STATISTIC calculated is the so called D statistic.",
                     "\nFor similar distributions the statistic converges to zero.",
                     "\n\nIf the p-value is less than 0.05",
                     "then we reject the null hypothesis and",
                     "\naccept the alternative hypothesis,",
                     "that the distributions differ, at the",
                     "95% level of confidence.")
  }
  else if (theWidget("test_ttest_radiobutton")$getActive())
  {
    test <- "locationTest"
    preamble <- Rtxt("Two-Sample t-Test",
                     "\n\nThe two-sample T-test is performed on the two specified samples. The",
                     "\nnull hypothesis is that the difference between the two means is zero.",
                     "\n\nThis test assumes that the two samples are normally distributed. If not,",
                     "\nuse the Wilcoxon Rank-Sum test.",
                     "\n\nThe confidence interval is an interval around the expected difference",
                     "\nbetween the means.",
                     "\n\nIf the p-value is less than 0.05 then we reject the null hypothesis and",
                     "\naccpet the alternative hypothesis, that the means differ, at the 95% level",
                     "\nof confidence.",
                     "\n\nTwo variants of the test are reported:",
                     "for equal and unequal variances.")
  }
  else if (theWidget("test_kw_radiobutton")$getActive())
  {
    test <- "locationTest"
    options <- ', method="kw"'
    preamble <- Rtxt("Kruskal-Wallis Test",
                     "\n\nThe Kruskal-Wallis test is performed on the two samples to test the",
                     "\nhypothesis that the difference between the two means is zero. It does not",
                     "\nassume that the two samples are normally distributed.",
                     "\n\nThe confidence interval is an interval around the expected difference",
                     "\nbetween the means.")
  }
  else if (theWidget("test_wilcoxon_radiobutton")$getActive())
  {
    test <- "wilcox.test"
    # options <- ", conf.int=TRUE" Could do this but needs more explanation
    preamble <- Rtxt("Wilcoxon Rank Sum Test",
                     "\n\nThe two-sample non-parametric Wilcoxon rank sum test (equivalent to",
                     "\nthe Mann-Whitney test) is performed on the two specified samples. The null",
                     "\nhypothesis is that the distributions are the same (i.e., there is no",
                     "\nshift in the location of the two distributions) with an alternative",
                     "\nhypothesis that they differ on location (based on median).",
                     "\n\nThis test does not assume that the two samples are normally distributed",
                     "\nbut does assume they have distributions of the same shape.",
                     "\n\nIf the p-value is less than 0.05 then we reject the null hypothesis and",
                     "\naccept the alternative hypothesis, that the two samples have different medians,",
                     "\nat the 95% level of confidence.")
  }
  else if (theWidget("test_wilcoxon_signed_radiobutton")$getActive())
  {
    test <- "wilcox.test"
    options <- ", paired=TRUE"
    if (! is.null(union(attr(na.omit(crs$dataset[,v1]), "na.action"),
                        attr(na.omit(crs$dataset[,v2]), "na.action"))))
    {
      test <- paste('miss <- union(',
                    'attr(na.omit(crs$dataset[,"', v1, '"]), "na.action"), ',
                    'attr(na.omit(crs$dataset[,"', v2, '"]), "na.action"))\n',
                    test, sep="")
      s1 <- s2 <- "-miss"
    }
    preamble <- Rtxt("Wilcoxon Signed Rank Test",
                     "\n\nThe paired sample non-parametric Wilcoxon signed rank test is",
                     "\nperformed on the two specified samples. The two samples are expected to be",
                     "\npaired (two observations for the same entity). The null hypothesis is that",
                     "\nthe distributions are the same.",
                     "\n\nThis test does not assume that the two samples are are normally distributed.",
                     "\n\nIf the p-value is less than 0.05 then we reject the null hypothesis and",
                     "\naccept the alternative hypothesis, that the distributions differ, at the",
                     "\n95% level of confidence.")
  }
  else if (theWidget("test_variance_radiobutton")$getActive())
  {
    test <- "varianceTest"
    preamble <- Rtxt("Two-Sample F-Test",
                     "\n\nThe two-sample F-test is performed on the two specified samples. The",
                     "\nnull hypothesis is that the ratio of the variances of the populations from",
                     "\nwhich they were drawn is equal to one.",
                     "\n\nThis test assumes that the two samples are normally distributed.",
                     "\n\nIf the p-value is less than 0.05 then we reject the null hypothesis and",
                     "\naccept the alternative hypothesis, that the two samples have different",
                     "\nvariances, at the 95% level of confidence.")
  }
  else if (theWidget("test_correlation_radiobutton")$getActive())
  {
    test <- "correlationTest"

    if (! is.null(union(attr(na.omit(crs$dataset[,v1]), "na.action"),
                        attr(na.omit(crs$dataset[,v2]), "na.action"))))
    {
      test <- paste('miss <- union(',
                    'attr(na.omit(crs$dataset[,"', v1, '"]), "na.action"), ',
                    'attr(na.omit(crs$dataset[,"', v2, '"]), "na.action"))\n',
                    test, sep="")
      s1 <- s2 <- "-miss"
    }

    preamble <- Rtxt("Correlation Test",
                     "\n\nThe paired sample correlation test is performed on the two specified samples.",
                     "\nThe two samples are expected to be paired (two observations for the same entity).",
                     "\nThe null hypothesis is that the two samples have no (i.e., 0) correlation.",
                     "\nPearson's product moment correlation coefficient is used.",
                     "\n\nIf the p-value is less than 0.05 then we reject the null hypothesis and",
                     "\naccept the alternative hypothesis that the samples are correlated,",
                     "\nat the 95% level of confidence.")
  }
  
  test.cmd <- sprintf(paste('%s(na.omit(crs$dataset[%s, "%s"]),',
                            'na.omit(crs$dataset[%s, "%s"])%s)'),
                      test, s1, v1, s2, v2, options)
  appendLog(Rtxt("Perform the test."), test.cmd)
  resetTextview(TV, preamble, "\n\n", msg, collectOutput(test.cmd))

  setStatusBar(Rtxt("Test completed."))
}

  
