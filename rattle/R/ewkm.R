# Gnome R Data Miner: GNOME interface to R for Data Mining
#
# Time-stamp: <2015-07-09 19:13:01 gjw>
#
# Implement biclust functionality.
#
# Copyright (c) 2011-2013 Togaware Pty Ltd
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

# When a radio button is selected, display the appropriate tab page.

on_ewkm_radiobutton_toggled <- function(button)
{
  if (button$getActive())
    crv$CLUSTER$setCurrentPage(crv$CLUSTER.KMEANS.TAB)
  setStatusBar()
}

on_kmeans_weights_plot_button_clicked <- function(button)
{
  weightsPlotEwkm()
}


########################################################################
# Execution

executeClusterEwkm <- function(include)
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
  
  startLog(commonName(crv$EWKM))

  # Load the required package.
  
  lib.cmd <- "library(wskm, quietly=TRUE)"
  if (! packageIsAvailable("wskm", Rtxt("perform subspace cluster analysis")))
    return(FALSE)
  appendLog(packageProvides('wskm', 'ewkm'), lib.cmd)
  eval(parse(text=lib.cmd))

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
    if (! packageIsAvailable("reshape", Rtxt("rescale for ewkm"))) return(FALSE)
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

  ewkm.cmd <- sprintf('crs$kmeans <- ewkm(%s, %s)', ds, centers)
    
  appendLog(sprintf(Rtxt("Generate a ewkm cluster of size %s."), nclust),
            ewkm.cmd)

  start.time <- Sys.time()

  result <- try(eval(parse(text=ewkm.cmd)), TRUE)
  time.taken <- Sys.time()-start.time

  # Show the resulting model.

  size.cmd <- "paste(crs$kmeans$size, collapse=' ')"
  means.cmd <- sprintf("colMeans(%s)", ds)
  centres.cmd <- "crs$kmeans$centers"
  withinss.cmd <- "crs$kmeans$withinss"
  weights.cmd <- "round(crs$kmeans$weights, 2)"
    
  startLog(Rtxt("Report on the cluster characteristics."))
  appendLog(Rtxt("Cluster sizes:"), size.cmd)
  appendLog(Rtxt("Data means:"), means.cmd)
  appendLog(Rtxt("Cluster centers:"), centres.cmd)
  appendLog(Rtxt("Cluster weights:"), weights.cmd)
  appendLog(Rtxt("Within cluster sum of squares:"), withinss.cmd)

  resetTextview(TV)
  setTextview(TV,
              sprintf(Rtxt("EWKM: %d clusters, %d iterations,",
                           "%d restarts, %d total iterations."),
                      length(crs$kmeans$size),
                      crs$kmeans$iterations,
                      crs$kmeans$restarts,
                      crs$kmeans$total.iterations),
              "\n\n",
              ifelse(crs$kmeans$restarts > 2,
                     Rtxt("With that many restarts you may want to",
                          "reduce the number of clusters.\n\n"), ""),
              Rtxt("Cluster sizes:"), "\n\n",
              collectOutput(size.cmd, TRUE),
              "\n\n", Rtxt("Data means:"), "\n\n",
              collectOutput(means.cmd),
              "\n\n", Rtxt("Cluster centers:"), "\n\n",
              collectOutput(centres.cmd, TRUE),
              "\n\n", Rtxt("Cluster weights:"), "\n\n",
              collectOutput(weights.cmd, TRUE),
              "\n\n", Rtxt("Within cluster sum of squares:"), "\n\n",
              collectOutput(withinss.cmd, TRUE),
              "\n")

  # Ensure the kmeans information buttons are now active.

  showModelKMeansExists()

  reportTimeTaken(TV, time.taken, model=commonName(crv$KMEANS))

  return(TRUE)
}

########################################################################
# Report on the model.

weightsPlotEwkm <- function()
{
  # Make sure there is a cluster first.

  if (is.null(crs$kmeans) || ! "ewkm" %in% class(crs$kmeans))
  {
    errorDialog("E126: No ewkm cluster to plot.",
                "The button should not have been sensitive.",
                crv$support.msg)
    return()
  }

  startLog(sprintf("Plot variable weights from the %s algorithm.", commonName(crv$EWKM)))

  # The wskm package provides the plot and levelplot methods.
  
  if (!packageIsAvailable("wskm", "plot variable weights")) return()
  lib.cmd <- "library(wskm, quietly=TRUE)"
  appendLog(packageProvides("wskm", "plot"), lib.cmd)
  eval(parse(text=lib.cmd))

  advancedPlot <- theWidget("use_ggplot2")$getActive() # Not really ggplot2 but convenient.

  if (advancedPlot)
    plot.cmd <- "plot(levelplot(crs$kmeans))"
  else
    plot.cmd <- "plot(crs$kmeans)"
  
  appendLog(Rtxt("Plot the variable weights."), plot.cmd)
  newPlot()
  eval(parse(text=plot.cmd))

  setStatusBar(Rtxt("The variable weights have been plotted."))
}

########################################################################
# Export

exportEwkmTab <- function()
{
  # Make sure we have a model first!

  if (noModelAvailable(crs$kmeans, crv$EWKM)) return(FALSE)

  startLog(paste(Rtxt("Export"), commonName(crv$EWKM)))
  
  save.name <- getExportSaveName(crv$EWKM)
  if (is.null(save.name)) return(FALSE)
  ext <- tolower(get.extension(save.name))

  # Construct the command to produce PMML.

  pmml.cmd <- sprintf(paste("pmml(crs$kmeans%s, description='%s',",
                            "algorithm.name='EWKM: Liping, Ng, and Huang')"),
                      ifelse(length(crs$transforms) > 0,
                             ", transforms=crs$transforms", ""),
                      commonName(crv$EWKM))

  # We can't pass "\" in a filename to the parse command in
  # MS/Windows so we have to run the save/write command separately,
  # i.e., not inside the string that is being parsed.

  if (ext == "xml")
  {
    appendLog(sprintf(Rtxt("Export %s as PMML."), commonName(crv$EWKM)),
              sprintf('saveXML(%s, "%s")', pmml.cmd, save.name))
    XML::saveXML(eval(parse(text=pmml.cmd)), save.name)
  }
  else if (ext == "c")
  {
    # 090103 gjw Move to a function: saveC(pmml.cmd, save.name,
    # "ewkm")

    # 090223 Why is this tolower being used? Under GNU/Linux it is
    # blatantly wrong. Maybe only needed for MS/Widnows

    if (isWindows()) save.name <- tolower(save.name)
    
    model.name <- sub("\\.c", "", basename(save.name))

    export.cmd <- generateExportPMMLtoC(model.name, save.name, "kmeans_textview")
    
    appendLog(sprintf(Rtxt("Export %s as a C routine."), commonName(crv$EWKM)),
              sprintf('pmml.cmd <- "%s"\n\n', pmml.cmd),
              export.cmd)

    eval(parse(text=export.cmd))
  }
  
  setStatusBar(sprintf(Rtxt("The model has been exported to '%s'."), save.name))

}

