# Gnome R Data Miner: GNOME interface to R for Data Mining
#
# Time-stamp: <2015-08-21 19:01:23 gjw>
#
# Implement kmeans functionality.
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
# CALLBACKS

# When a radio button is selected, display the appropriate tab page.

on_kmeans_radiobutton_toggled <- function(button)
{
  if (button$getActive())
  {
    crv$CLUSTER$setCurrentPage(crv$CLUSTER.KMEANS.TAB)
    if (not.null(crs$hclust))
      theWidget("kmeans_hclust_centers_checkbutton")$setSensitive(TRUE)
    else
      theWidget("kmeans_hclust_centers_checkbutton")$setSensitive(FALSE)
  }
  setStatusBar()
}

on_kmeans_hclust_centers_checkbutton_toggled <- function(button)
{
  #
  # When the hclust centers checkbutton is on, we should not allow any
  # use of the runs spin button because it does not make sens to do
  # both, nor does kmeansruns support starting points for clusters.
  #
  if (button$getActive())
  {
    theWidget("kmeans_runs_spinbutton")$setSensitive(FALSE)
    theWidget("kmeans_runs_label")$setSensitive(FALSE)
  }
  else
  {
    theWidget("kmeans_runs_spinbutton")$setSensitive(TRUE)
    theWidget("kmeans_runs_label")$setSensitive(TRUE)
  }
}    

on_kmeans_iterate_checkbutton_toggled <- function(button)
{
  #
  # When the iterate checkbutton is on, we should not allow any use of
  # the runs spin button because it does not make sense to do both.
  #
  if (button$getActive())
  {
    theWidget("kmeans_runs_spinbutton")$setSensitive(FALSE)
    theWidget("kmeans_runs_label")$setSensitive(FALSE)
  }
  else
  {
    theWidget("kmeans_runs_spinbutton")$setSensitive(TRUE)
    theWidget("kmeans_runs_label")$setSensitive(TRUE)
  }
}    

on_kmeans_seed_button_clicked <- function(button)
{
  rseed <- as.integer(runif(1, 0, 1000000))
  theWidget("kmeans_seed_spinbutton")$setValue(rseed)
}

on_kmeans_stats_button_clicked <- function(button)
{
  displayClusterStatsKMeans()
}

on_kmeans_data_plot_button_clicked <- function(button)
{
  dataPlotKMeans()
}

on_kmeans_discriminant_plot_button_clicked <- function(button)
{
  discriminantPlotKMeans()
}

########################################################################
# Execution

executeClusterKMeans <- function(include)
{
  TV <- "kmeans_textview"
  sampling  <- not.null(crs$sample)

  # Obtain interface information.
  
  nclust <- theWidget("kmeans_clusters_spinbutton")$getValue()
  seed <- theWidget("kmeans_seed_spinbutton")$getValue()
  if (seed == crv$seed) seed <- "crv$seed"
  nruns <- theWidget("kmeans_runs_spinbutton")$getValue()
  usehclust <- theWidget("kmeans_hclust_centers_checkbutton")$getActive()
  useIterate <- theWidget("kmeans_iterate_checkbutton")$getActive()
  
  startLog(commonName(crv$KMEANS))

  # Set the seed so we can repeat.

  seed.cmd <- sprintf('set.seed(%s)', seed)
  appendLog(Rtxt("Reset the random number seed to obtain the same results each time."),
            seed.cmd)
  eval(parse(text=seed.cmd))

  # Determine the dataset to use.

  ds <- sprintf("na.omit(crs$dataset[%s, %s])",
                ifelse(sampling, "crs$sample", ""), include)

  # Check if we should rescale

  if (theWidget("kmeans_rescale_checkbutton")$getActive())
  {
    lib.cmd <- "library(reshape, quietly=TRUE)"
    if (! packageIsAvailable("reshape", Rtxt("rescale for kmeans"))) return(FALSE)
    appendLog(packageProvides('reshape', 'rescaler'), lib.cmd)
    eval(parse(text=lib.cmd))
    
    ds <- sprintf('sapply(%s, rescaler, "range")', ds)
  }
  
  # Calculate the centers

  if (usehclust)
    centers <- sprintf("centers.hclust(%s, crs$hclust, %d)", ds, nclust)
  else
    centers <- nclust
  
  # KMEANS: Log the R command and execute.

  if (! useIterate)
  {
    if (nruns > 1)
    {
      lib.cmd <- "library(fpc, quietly=TRUE)"
      if (! packageIsAvailable("fpc", Rtxt("run kmeans multiple times"))) return(FALSE)
      appendLog(packageProvides('fpc', 'kmeansruns'), lib.cmd)
      eval(parse(text=lib.cmd))

      kmeans.cmd <- sprintf('crs$kmeans <- kmeansruns(%s, %s, runs=%s)',
                            ds, centers, nruns)
    }
    else
    {
      kmeans.cmd <- sprintf('crs$kmeans <- kmeans(%s, %s)', ds, centers)
    }
    
    appendLog(sprintf(Rtxt("Generate a kmeans cluster of size %s%s%s."), nclust,
                      ifelse(nruns>1, Rtxt(" choosing the best from "), ""),
                      ifelse(nruns>1, nruns, "")),
              kmeans.cmd)

    start.time <- Sys.time()

    result <- try(eval(parse(text=kmeans.cmd)), TRUE)
    time.taken <- Sys.time()-start.time

    if (inherits(result, "try-error"))
    {
      if (any(grep("more cluster centers than distinct data points", result))||
          any(grep("cannot take a sample larger than the population", result)))
        errorDialog(Rtxt("The data does not support the number of clusters",
                         "requested. Reduce the number of clusters and try again."))
      else
        errorDialog(errorMessageFun("kmeans", result))
      return(FALSE)
    }

    # Show the resulting model.

    size.cmd <- "paste(crs$kmeans$size, collapse=' ')"
    means.cmd <- sprintf("colMeans(%s)", ds)
    centres.cmd <- "crs$kmeans$centers"
    withinss.cmd <- "crs$kmeans$withinss"
    
    startLog(Rtxt("Report on the cluster characteristics."))
    appendLog(Rtxt("Cluster sizes:"), size.cmd)
    appendLog(Rtxt("Data means:"), means.cmd)
    appendLog(Rtxt("Cluster centers:"), centres.cmd)
    appendLog(Rtxt("Within cluster sum of squares:"), withinss.cmd)

    resetTextview(TV)
    setTextview(TV, Rtxt("Cluster sizes:"), "\n\n",
                collectOutput(size.cmd, TRUE),
                "\n\n", Rtxt("Data means:"), "\n\n",
                collectOutput(means.cmd),
                "\n\n", Rtxt("Cluster centers:"), "\n\n",
                collectOutput(centres.cmd, TRUE),
                "\n\n", Rtxt("Within cluster sum of squares:"), "\n\n",
                collectOutput(withinss.cmd, TRUE),
                "\n")

    # Ensure the kmeans information buttons are now active.

    showModelKMeansExists()

  }
  else # Iterate over the clusters.
  {
    start.time <- Sys.time()
    css <- vector()
    for (i in 2:nclust)
    {
      kmeans.cmd <- sprintf(paste('crs$kmeans <-',
                                  'kmeans(na.omit(crs$dataset[%s, %s]), %s)'),
                            ifelse(sampling, "crs$sample", ""), include, i)
      eval(parse(text=seed.cmd))
      eval(parse(text=kmeans.cmd))
      css[i] <- sum(crs$kmeans$withinss)
    }
    css[1] <- crs$kmeans$totss
    time.taken <- Sys.time()-start.time
    resetTextview(TV)
    setTextview(TV, sprintf(Rtxt("We have iterated over cluster sizes",
                                 "from 2 to %d clusters.\n",
                                 "\nThe plot displays the 'sum(withinss)' for each clustering",
                                 "\nand the change in this value from the previous clustering.\n"),
                            nclust))
    newPlot()
    plot(1:nclust, c(css[1:nclust]), ylim=c(0, max(css[1:nclust])),
         type="b", lty=1, col="blue",
         xlab=Rtxt("Number of Clusters"), ylab=Rtxt("Sum of WithinSS"),
         main=Rtxt("Sum of WithinSS Over Number of Clusters"))
    points(2:nclust, css[1:(nclust-1)]-css[2:nclust],
           type="b", pch=4, lty=2, col="red")
    legend("topright", c(Rtxt("Sum(WithinSS)"), Rtxt("Diff previous Sum(WithinSS)")),
           col=c("blue", "red"), lty=c(1, 2), pch=c(1,4), inset=0.05)
  }

  reportTimeTaken(TV, time.taken, model=commonName(crv$KMEANS))

  return(TRUE)
}

# 100224 Added to support resetting buttons on loading a
# project. Modelled on showModelExistsRPart.

showModelKMeansExists <- function(state=!is.null(crs$kmeans))
{
  # If a kmeans model exists then make available the Stats, Data Plot,
  # and Discriminate Plot buttons on the Descriptive tab.
  
  theWidget("kmeans_stats_button")$setSensitive(state)
  theWidget("kmeans_data_plot_button")$setSensitive(state)
  theWidget("kmeans_discriminant_plot_button")$setSensitive(state)

  # 110911 TODO Might want to move this into a showModelEwkmExists and
  # perhaps even save the textview and swap between kmeans and ewkm in
  # the textview.
  
  theWidget("kmeans_weights_plot_button")$setSensitive(state &&
                                                       theWidget("ewkm_radiobutton")$getActive())
  
}

########################################################################
# Export

exportKMeansTab <- function()
{
  # Make sure we have a model first!

  if (noModelAvailable(crs$kmeans, crv$KMEANS)) return(FALSE)

  startLog(paste(Rtxt("Export"), commonName(crv$KMEANS)))
  
  save.name <- getExportSaveName(crv$KMEANS)
  if (is.null(save.name)) return(FALSE)
  ext <- tolower(get.extension(save.name))

  # Construct the command to produce PMML.

  pmml.cmd <- sprintf("pmml(crs$kmeans%s)",
                      ifelse(length(crs$transforms) > 0,
                             ", transforms=crs$transforms", ""))

  # We can't pass "\" in a filename to the parse command in
  # MS/Windows so we have to run the save/write command separately,
  # i.e., not inside the string that is being parsed.

  if (ext == "xml")
  {
    appendLog(sprintf(Rtxt("Export %s as PMML."), commonName(crv$KMEANS)),
              sprintf('saveXML(%s, "%s")', pmml.cmd, save.name))
    XML::saveXML(eval(parse(text=pmml.cmd)), save.name)
  }
  else if (ext == "c")
  {
    # 090103 gjw Move to a function: saveC(pmml.cmd, save.name,
    # "kmeans")

    # 090223 Why is this tolower being used? Under GNU/Linux it is
    # blatantly wrong. Maybe only needed for MS/Widnows

    if (isWindows()) save.name <- tolower(save.name)
    
    model.name <- sub("\\.c", "", basename(save.name))

    export.cmd <- generateExportPMMLtoC(model.name, save.name, "kmeans_textview")
    
    appendLog(sprintf(Rtxt("Export %s as a C routine."), commonName(crv$KMEANS)),
              sprintf('pmml.cmd <- "%s"\n\n', pmml.cmd),
              export.cmd)

    eval(parse(text=export.cmd))
  }
  
  setStatusBar(sprintf(Rtxt("The model has been exported to '%s'."), save.name))

}

########################################################################
# Score

predict.kmeans <- function(object, data, ...)
{
  # 081228 Initial work on a predict.kmeans function, to allow using a
  # kmeans model to allocate new data to pre-existing clusters using
  # the common model interface function, predict. This makes it easy
  # to use the Rattle modelling code on kmeans. TODO Currently, no
  # support for alternative distance measures. This will be needed
  # eventually
  
  #num.clusters <- nrow(object$centers)
  cluster.names <- rownames(object$centers)
  cluster.vars <- colnames(object$centers)
  #num.rows <- nrow(data)
  #cluster.row.nums <- seq(num.rows+1, num.rows+num.clusters)

  # 081228 Put the data first, to maintain rownames. If there are
  # conflicts in rownames then rbind creates new rownames, and this
  # will be a surprise for the rownames of the returned data. We avoid
  # this by having the original data first, thus its rownames are
  # maintained. This may well change the cluster names, which we then
  # need to change back correctly.

  ## REMOVE 090808
  ## # 081228 Simply calculate the distance between all points - this is
  ## # simpler to code, but perhaps less efficient?
  
  ## d <- as.matrix(dist(rbind(data[cluster.vars], object$centers)))
  ## d <- d[-cluster.row.nums,cluster.row.nums]
  ## colnames(d) <- cluster.names
  
  ## out <- apply(d, 1, which.min)
  ## miss <- attr(na.omit(data[cluster.vars]), "na.action")
  ## out[miss] <- NA
  ## return(out)

  # 090808 Because datasets can be large try to use an efficient
  # approach.

  out <- apply(data[cluster.vars], 1,
               function(d)
               cluster.names[which.min(apply(object$centers, 1,
                                             function(x)
                                             sqrt(sum(abs(d - x)^2))))])
  out <- sapply(out, function(x) ifelse(length(x), x, NA))
  return(out)
}

genPredictKmeans <- function(dataset)
{
  # 081227 Generate a command to obtain the prediction results when
  # applying the model to new data.
  
  return(sprintf("crs$pr <- predict(crs$kmeans, %s)", dataset))
}

genResponseKmeans <- function(dataset)
{
  # 081227 Generate a command to obtain the response when applying the
  # model to new data.
  
  return(genPredictKmeans(dataset))
}

genProbabilityKmeans <- function(dataset)
{
  # 081227 Generate a command to obtain the probability when applying
  # the model to new data. There is probably a prblem with simply
  # using the cluster label as the output, since it won't be a
  # probability or even look like it. Let's do it for now though -
  # should be okay.
  
  return(genPredictKmeans(dataset))
}

########################################################################
# Report on the model.

displayClusterStatsKMeans <- function()
{
  # Make sure there is a cluster first.
  
  if (is.null(crs$kmeans))
  {
    errorDialog("E124: Should not be here.", crv$support.msg)
    return()
  }

  startLog(Rtxt("Generate statistics for the clustering."))

  # LIBRARY: Ensure the appropriate package is available for the
  # plot, and log the R command and execute.
  
  if (!packageIsAvailable("fpc", Rtxt("plot a cluster"))) return()
  lib.cmd <- "library(fpc, quietly=TRUE)"
  appendLog(packageProvides("fpc", "cluster.stats"), lib.cmd)
  eval(parse(text=lib.cmd))

  # Some background information.  Assume we have already built the
  # cluster, and so we don't need to check so many conditions.

  TV <- "kmeans_textview"
  sampling  <- not.null(crs$sample)

  # 091219 Why would we check this here - the cluster is already
  # built?
  
  include <- "crs$numeric" # 20110102 getNumericVariables()
  if (length(include) == 0)
  {
    errorDialog(Rtxt("Clusters are currently calculated only for numeric data.",
                     "No numeric variables were found in the dataset",
                     "from amongst those having an input/target/risk role."))
    return()
  }

  # 091129 For large data, the distance matrix gets very large and
  # calculations take a long time. Under 32 bit we often run out of
  # memory. Under 64 bit we might eventually run out of virtual
  # memory. Pop up a warning.

  large <- length(crs$kmeans$cluster) > crv$cluster.report.max.obs
  
  ## if (large &&
  ##     ! questionDialog("The dataset has a large number of observations.",
  ##                      "This may be too large to calculate the",
  ##                      "required pairwsie distance matrix for this dataset.",
  ##                      "If you continue you may find that the command will",
  ##                      "fail on trying to allocate memory (on a 32 bit",
  ##                      "computer) or else will proceed to fill up all available",
  ##                      "memory. Whlist the calculations are being performed",
  ##                      "the interface will not be responsive and you may not",
  ##                      "be able to interrupt the process.",
  ##                      "Consider saving your project before proceeding.",
  ##                      "\n\nDo you wish to continue anyhow?"))
  ##   return(FALSE)
  
  if (large &&
      ! questionDialog(sprintf(Rtxt("The dataset contains many observations.",
                                    "The statistics are based on a",
                                    "pairwise distance matrix which can be enormous",
                                    "(GBs for 10,000 observations).",
                                    "If you continue, %s will use an",
                                    "auto sampling methodology of size",
                                    "%d to calculate the statistics.",
                                    "\n\nWould you like to continue with auto sampling?"),
                               crv$appname, crv$cluster.report.max.obs)))
      return(FALSE)
  
  # STATS: Log the R command and execute. 080521 TODO Fix a bug by
  # adding the na.omit here (since by default that is done in building
  # the clusters). Not sure if this is generally correct.

  if (large)
  {
    large.sample.cmd <- paste(sprintf("set.seed(%s)", crv$seed),
                              sprintf("smpl <<- sample(length(crs$kmeans$cluster), %d)",
                                      crv$cluster.report.max.obs),    
                              sep="\n")
    appendLog(Rtxt("Select a sample from the dataset to calculate the statistics."),
              sub("<<", "<", large.sample.cmd))
    eval(parse(text=large.sample.cmd))
  }
  
  set.cursor("watch", Rtxt("Determining the cluster statistics...."))
  on.exit(set.cursor("left-ptr"))
  while (RGtk2::gtkEventsPending()) RGtk2::gtkMainIteration()

  stats.cmd <- sprintf(paste("cluster.stats(dist(na.omit(crs$dataset[%s, %s]%s)),",
                             "crs$kmeans$cluster%s)\n"),
                       ifelse(sampling, "crs$sample", ""),
                       include,
                       ifelse(large, "[smpl,]", ""),
                       ifelse(large, "[smpl]", ""))
  appendLog(packageProvides("fpc", "cluster.stats"), stats.cmd)
  result <- try(collectOutput(stats.cmd, use.print=TRUE))
  if (inherits(result, "try-error"))
  {
    if (any(grep("[cC]annot allocate (vector|memory)", result)))
    {
      errorDialog("E144: The call to cluster.stats appears to have failed.",
                  "This is often due, as is the case here,",
                  "to running out of memory.",
                  "A quick solution is to sample the dataset, through the",
                  "Data tab, and rebuild the cluster. On 32bit machines you may be limited to",
                  "less than 4000 observations.")
      setTextview(TV)
    }
    else
      errorDialog(errorMessageFun("cluster.stat", result))
    return(FALSE)
  }

  appendTextview(TV, Rtxt("General cluster statistics:"), "\n\n",
                 result)
  setStatusBar(paste(Rtxt("KMeans cluster statistics have been generated."),
                     Rtxt("You may need to scroll the textview to view them.")))
}

dataPlotKMeans <- function()
{
  
  # Make sure there is a cluster first.

  if (is.null(crs$kmeans))
  {
    errorDialog("E132: No cluster model. Please create one first.", crv$support.msg)
    return()
  }

  startLog(Rtxt("Display a scatterplot matrix for the KMeans clustering."))

  # Some background information.  Assume we have already built the
  # cluster, and so we don't need to check so many conditions.

  sampling  <- not.null(crs$sample)

  include <- "intersect(crs$input, crs$numeric)"
  incvars <- eval(parse(text=include))

  # We can only plot if there is more than a single variable.
  
  if (length(incvars) == 1)
  {
    infoDialog(Rtxt("A data plot of the clusters can not be constructed",
                    "because there is only one numeric variable available",
                    "in the data."))
    return()
  }

  # 091219 Check for very large data, and if so use auto sampling.
  
  large <- length(crs$kmeans$cluster) > crv$cluster.report.max.obs
  manyvars <- length(incvars) > crv$scatter.max.vars

  if (large && manyvars)
  {
    if (! questionDialog(sprintf(Rtxt("The dataset contains many variables and observations.",
                                      "For more than %d variables and %d",
                                      "observations the plot will be cluttered and quite slow.",
                                      "The application will also wait whilst",
                                      "the plot is drawn.",
                                      "If you continue, the first %d variables will be used and",
                                      "%s will use an auto sampling methodology of size %d",
                                      "to reduce the number of observations.",
                                      "\n\nWould you like to continue with the reduced",
                                      "variables and auto sampling?"),
                                 crv$scatter.max.vars, crv$cluster.report.max.obs,
                                 crv$scatter.max.vars, crv$appname, crv$cluster.report.max.obs)))
      return(FALSE)
  }
  else if (large)
  {
    if (! questionDialog(sprintf(Rtxt("The dataset contains many observations.",
                                      "For more than %d",
                                      "observations the plot will be quite slow.",
                                      "The application will also wait whilst",
                                      "the plot is drawn. If you continue, %s",
                                      "will use an auto sampling methodology of size %d",
                                      "to reduce the number of observations.",
                                      "\n\nWould you like to continue with auto sampling?"),
                                 crv$cluster.report.max.obs, crv$appname,
                                 crv$cluster.report.max.obs)))
      return(FALSE)
  }
  else if (manyvars)
  {
    if (! questionDialog(sprintf(Rtxt("The dataset contains many variables.",
                                      "For more than %d",
                                      "variables the plot will be cluttered.",
                                      "The application will also wait whilst",
                                      "the plot is drawn.",
                                      "If you continue, the first %d",
                                      "variables will be used.",
                                      "\n\nWould you like to continue with the reduced",
                                      "variables?"),
                                 crv$scatter.max.vars, crv$scatter.max.vars)))
      return(FALSE)
  }
    
  if (large)
  {
    large.sample.cmd <- paste(sprintf("set.seed(%s)", crv$seed),
                              sprintf("smpl <<- sample(length(crs$kmeans$cluster), %d)",
                                      crv$cluster.report.max.obs),
                              sep="\n")
    appendLog(Rtxt("Select a sample from the dataset to calculate the statistics."),
              sub("<<", "<", large.sample.cmd))
    eval(parse(text=large.sample.cmd))
  }

  if (manyvars)
  {
    top.vars.cmd <- sprintf("vars <<- 1:%d", crv$scatter.max.vars)
    appendLog(sprintf(Rtxt("Keep just the first %d variables for the plot."),
                      crv$scatter.max.vars),
              sub("<<", "<", top.vars.cmd))
    eval(parse(text=top.vars.cmd))
  }
  
  # PLOT: Log the R command and execute. 080521 TODO I've added in
  # na.omit here, since when we cluster the audit data, with missing
  # values for Age we need to ensure the data points correspond to the
  # cluster numbers. Otherwise we get a bad looking plot!!!! But do we
  # always need na.omit. It is not always used on bulding clusters.

  # Alternative plot commands that could be considered:
  #
  #    plot3d (rgl) with type="s"
  #    plotmatrix (ggplot2)
  #    splom (lattice)
  #
  # I think the default plot is quite good. plotmatrix is good, but
  # does not include the scales and takes a long time to render.

  ##  plot.cmd <- sprintf(paste("plot(crs$dataset[%s,%s], ",

  plot.cmd <- sprintf(paste("plot(na.omit(crs$dataset[%s, %s]%s), ",
                            "col=crs$kmeans$cluster)\n%s", sep=""),
                      ifelse(sampling, "crs$sample", ""), include,
                      ifelse(large,
                             ifelse(manyvars, "[smpl, vars]", "[smpl,]"),
                             ifelse(manyvars, "[vars]", "")),
                      genPlotTitleCmd(""))
  appendLog(Rtxt("Generate a data plot."), plot.cmd)

  set.cursor("watch", Rtxt("Rendering the plot. Please wait...."))
  newPlot()
  eval(parse(text=plot.cmd))
  set.cursor("left-ptr", Rtxt("The data plot has been generated."))
}

discriminantPlotKMeans <- function()
{

  # Make sure there is a cluster first.

  if (is.null(crs$kmeans))
  {
    errorDialog("E125: No cluster to plot.",
                "The button should not have been sensitive.",
                crv$support.msg)
    return()
  }

  # Some background information.  Assume we have already built the
  # cluster, and so we don't need to check so many conditions.

  sampling <- not.null(crs$sample)

  include <- "intersect(crs$input, crs$numeric)"
  incvars <- eval(parse(text=include))

  if (length(crs$numeric) == 0 || length(crs$input) == 0)
  {
    errorDialog(Rtxt("Clusters are currently calculated only for numeric data.",
                     "No numeric variables were found in the dataset",
                     "from amongst those having an input/target/risk role."))
    return()
  }

  # We can only plot if there is more than a single variable.
  
  if (length(incvars) == 1)
  {
    infoDialog(Rtxt("A discriminant coordinates plot can not be constructed",
                    "because there is only one numeric variable available."))
    return()
  }

  # PLOT: Log the R command and execute. 080521 Add the na.omit since
  # kmeans is usually built with this. 150821 Move to using
  # cluster::clusplot.

 plot.cmd <- paste(sprintf("cluster::clusplot(na.omit(crs$dataset[%s, %s]), ",
                            ifelse(sampling, "crs$sample", ""), include),
                   "crs$kmeans$cluster, color=TRUE, shade=TRUE, ",
                   "main='Discriminant Coordinates ",
                   crs$dataname, "')\n",
                   sep="")
  appendLog(Rtxt("Generate a discriminant coordinates plot."), plot.cmd)
  newPlot()
  eval(parse(text=plot.cmd))

  setStatusBar(Rtxt("A discriminant coordinates plot has been generated."))
}
