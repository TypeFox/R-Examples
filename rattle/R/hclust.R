# Gnome R Data Miner: GNOME interface to R for Data Mining
#
# Time-stamp: <2015-09-30 06:27:32 gjw>
#
# Implement hclust functionality.
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

# When a radio button is selected, display the appropriate tab page.

on_hclust_radiobutton_toggled <- function(button)
{
  if (button$getActive())
    crv$CLUSTER$setCurrentPage(crv$CLUSTER.HCLUST.TAB)
  setStatusBar()
}

on_hclust_dendrogram_button_clicked <- function(button)
{
  if (theWidget("use_ggplot2")$getActive())
    plotDendrogram2()
  else
    plotDendrogram()
}

on_hclust_stats_button_clicked <- function(button)
{
  set.cursor("watch", Rtxt("Determining the cluster statistics...."))
  try(displayHClustStats())
  set.cursor("left-ptr", Rtxt("Cluster statistics displayed. Scroll to see all."))
}

on_hclust_data_plot_button_clicked <- function(button)
{

  # Make sure there is a cluster first.

  if (is.null(crs$hclust))
  {
    errorDialog(Rtxt("No cluster to plot.",
                     "The button should not have been sensitive."),
                crv$support.msg)
    return()
  }

  # Some background information.  Assume we have already built the
  # cluster, and so we don't need to check so many conditions.

  sampling  <- not.null(crs$sample)
  num.clusters <- theWidget("hclust_clusters_spinbutton")$getValue()
  nums <- seq(1,ncol(crs$dataset))[as.logical(sapply(crs$dataset, is.numeric))]
  if (length(nums) > 0)
  {
    indicies <- getVariableIndicies(crs$input)
    include <- simplifyNumberList(intersect(nums, indicies))
  }

  if (length(nums) == 0 || length(indicies) == 0)
  {
    errorDialog(Rtxt("Clusters are currently calculated only for numeric data.",
                     "No numeric variables were found in the dataset",
                     "from amongst those having an input/target/risk role."))
    return()
  }

  # We can only plot if there is more than a single variable.
  
  if (length(intersect(nums, indicies)) == 1)
  {
    infoDialog(Rtxt("A data plot can not be constructed",
                    "because there is only one numeric variable available",
                    "in the data."))
    return()
  }

  # PLOT: Log the R command and execute.

  set.cursor("watch", Rtxt("Determining the cluster statistics...."))
  plot.cmd <- paste(sprintf(paste("plot(crs$dataset[%s, %s], ",
                                  "col=cutree(crs$hclust, %d))\n",
                                  sep=""),
                            ifelse (sampling, "crs$sample", ""), include,
                            num.clusters),
                    genPlotTitleCmd(""))
  appendLog(Rtxt("Generate a data plot."), plot.cmd)
  newPlot()
  eval(parse(text=plot.cmd))

  set.cursor("left-ptr", Rtxt("Data plot has been generated."))
}

on_hclust_discriminant_plot_button_clicked <- function(button)
{

  # Make sure there is a cluster first.

  if (is.null(crs$hclust))
  {
    errorDialog(Rtxt("No cluster to plot.",
                     "The button should not have been sensitive."),
                crv$support.msg)
    return()
  }

  # Some background information.  Assume we have already built the
  # cluster, and so we don't need to check so many conditions.

  sampling  <- not.null(crs$sample)
  num.clusters <- theWidget("hclust_clusters_spinbutton")$getValue()
  nums <- seq(1,ncol(crs$dataset))[as.logical(sapply(crs$dataset, is.numeric))]
  if (length(nums) > 0)
  {
    indicies <- getVariableIndicies(crs$input)
    include <- simplifyNumberList(intersect(nums, indicies))
  }

  if (length(nums) == 0 || length(indicies) == 0)
  {
    errorDialog(Rtxt("Clusters are currently calculated only for numeric data.",
                     "No numeric variables were found in the dataset",
                     "from amongst those having an input/target/risk role."))
    return()
  }

  # We can only plot if there is more than a single variable.
  
  if (length(intersect(nums, indicies)) == 1)
  {
    infoDialog(Rtxt("A discriminant coordinates plot can not be constructed",
                    "because there is only one numeric variable available",
                    "in the data."))
    return()
  }

  # PLOT: Log the R command and execute.

 plot.cmd <- paste(sprintf(paste("cluster::clusplot(na.omit(crs$dataset[%s, %s]), ",
                                 "cutree(crs$hclust, %d), ",
                                 "color=TRUE, shade=TRUE, ",
                                 "main='Discriminant Coordinates ",
                                 crs$dataname, "')\n", sep=""),
                           ifelse(sampling, "crs$sample", ""), include, num.clusters))
  appendLog(Rtxt("Generate a discriminant coordinates plot."), plot.cmd)
  newPlot()
  eval(parse(text=plot.cmd))

  setStatusBar(Rtxt("Discriminant coordinates plot has been generated."))
}

########################################################################
# EXECUTION

executeClusterHClust <- function(include)
{
  # Initial setup. Ensure the textview is monospace for fixed width
  # output of the centers and other information (so the columns line
  # up).

  TV <- "hclust_textview"
  theWidget(TV)$modifyFont(RGtk2::pangoFontDescriptionFromString(crv$textview.font))
  
  # TODO : If data is large put up a question about wanting to
  # continue?
  
  sampling  <- not.null(crs$sample)

  startLog(Rtxt("Hierarchical Cluster"))

  # The amap library needs to be loaded for hcluster. Also note that
  # hcluster takes about 0.33 seconds, compared to hclust taking 11
  # seconds!

  lib.cmd <- "library(amap, quietly=TRUE)"
  if (packageIsAvailable("amap", Rtxt("perform an efficient hierarchical clustering")))
  {
    amap.available <- TRUE
    appendLog(packageProvides("amap", "hclusterpar"), lib.cmd)
    eval(parse(text=lib.cmd))
  }
  else
    amap.available <- FALSE
  
  # Obtain interface information. 100319 Need to account for
  # transaltions in getting the function options.

  dist <- theWidget("hclust_distance_combobox")$getActiveText()
  Encoding(dist) <- "UTF-8"
  dist.opts <- paste("euclidean", "maximum", "manhattan", "canberra",
                     "binary", "pearson", "correlation", "spearman",
                     sep="\n")
  dist.rtxt <- Rtxt (dist.opts)
  dist.opts <- strsplit(dist.opts, "\n")[[1]]
  dist.rtxt <- strsplit(dist.rtxt, "\n")[[1]]
  dist <- dist.opts[which(dist == dist.rtxt)]
  
  link <- theWidget("hclust_link_combobox")$getActiveText()
  Encoding(link) <- "UTF-8"
  link.opts <- paste("complete", "ward", "single", "average", "mcquitty",
                     "median", "centroid", sep="\n")
  link.rtxt <- Rtxt (link.opts)
  link.opts <- strsplit(link.opts, "\n")[[1]]
  link.rtxt <- strsplit(link.rtxt, "\n")[[1]]
  link <- link.opts[which(link == link.rtxt)]
  
  nbproc <- theWidget("hclust_nbproc_spinbutton")$getValue()

  # Check if user has requested more than a single processor, and if
  # so but amap is not available, inform the user and exit.
  
  if (nbproc != 1 && ! amap.available)
  {
    errorDialog(Rtxt("The 'amap' package is not available and so the efficient",
                     "and parallel hcluster is not available.",
                     "Please set the number of processors to 1 to proceed",
                     "with using the single processor hclust instead.",
                     "Be aware that the amap version is over 10 times faster.",
                     "You may ant to install the 'amap' package."))
    return(FALSE)
  }
  
  # Determine which hclust to use for clustering.

  if (amap.available)

    # Use the more efficient hcluster for clustering.
  
    hclust.cmd <- paste("crs$hclust <- ",
                        sprintf(paste('hclusterpar(na.omit(crs$dataset[%s, %s]),',
                                      '\n    method="%s", link="%s",',
                                      'nbproc=%d)'),
                                ifelse(sampling, "crs$sample", ""),
                                include, dist, link, nbproc),
                        sep="")
  else

    # Use the standard hclust for clustering.
    
    hclust.cmd <- paste("crs$hclust <- ",
                        sprintf(paste('hclust(dist(crs$dataset[%s, %s],',
                                      'method="%s"),',
                                      'method="%s")'),
                                ifelse(sampling, "crs$sample", ""),
                                include, dist, link),
                        sep="")

  # Log the R command.

  appendLog(Rtxt("Generate a hierarchical cluster of the data."),
          hclust.cmd)
  
  # Perform the commands.

  start.time <- Sys.time()
  result <- try(eval(parse(text=hclust.cmd)), silent=TRUE)
  time.taken <- Sys.time()-start.time
  if (inherits(result, "try-error"))
  {
    if (any(grep("[cC]annot allocate (vector|memory)", result)))
    {
      errorDialog(Rtxt("The call to hclust appears to have failed.",
                       "This is often due, as is the case here,",
                       "to running out of memory",
                       "as hclust is rather memory hungry.",
                       "A quick solution is to sample the dataset, through the",
                       "Data tab. On 32bit machines you may be limited to",
                       "less than 2000 observations."))
      setTextview(TV)
    }
    else
      errorDialog(errorMessageFun("hclust", result))
    return(FALSE)
  }

  setTextview(TV, Rtxt("Hierachical Cluster"), "\n", collectOutput("crs$hclust", TRUE))

  showModelHClustExists()
  
  reportTimeTaken(TV, time.taken,
                  msg=Rtxt("A hierarchical cluster has been generated."))
  
  return(TRUE)
}

centers.hclust <- function(x, h, nclust=10, use.median=FALSE)
{
  if (!inherits(h, "hclust")) stop(Rtxt("Not a legitimate hclust object"))
    
  if (class(x) != "matrix") x <- as.matrix(x)
  if (use.median)
    centres <- round(tapply(x, list(rep(cutree(h, nclust), ncol(x)),
                                    col(x)), median))
  else
    centres <- tapply(x, list(rep(cutree(h, nclust), ncol(x)),
                              col(x)), mean)
  dimnames(centres) <- list(NULL, dimnames(x)[[2]])
  return(centres)
}

plotDendrogram2 <- function()
{

  # Make sure there is a hclust object first.

  if (is.null(crs$hclust))
  {
    errorDialog(Rtxt("No Hierarchical Cluster to plot."), crv$support.msg)
    return()
  }

  startLog(Rtxt("Dendrogram Plot"))
  
  # Load the required package into the library.

  lib.cmd <- "library(ggplot2, quietly=TRUE)"
  if (! packageIsAvailable("ggplot2", Rtxt("plot a dendrogram"))) return(FALSE)
  appendLog(packageProvides("ggplot2", "ggplot"), lib.cmd)
  eval(parse(text=lib.cmd))

  lib.cmd <- "library(ggdendro, quietly=TRUE)"
  if (! packageIsAvailable("ggdendro", Rtxt("plot a dendrogram"))) return(FALSE)
  appendLog(packageProvides("ggdendro", "dendro_data"), lib.cmd)
  eval(parse(text=lib.cmd))

  # Show a busy cursor whilst drawing the plot.

  set.cursor("watch", Rtxt("Rendering the hierarchical cluster dendrogram...."))
  on.exit(set.cursor("left-ptr", ""))
  
  ttl <- genPlotTitleCmd(Rtxt("Cluster Dendrogram"), crs$dataname, vector=TRUE)
  plot.cmd <- paste('ddata <- dendro_data(crs$hclust, type="rectangle")',
                    'g <- ggplot(segment(ddata))',
                    'g <- g + geom_segment(aes(x = y, y = x, xend = yend, yend = xend))',
                    'g <- g + scale_y_discrete(labels = ddata$label$label)',
                    'g <- g + labs(x="Height", y="Observation")',
                    paste("g <- g +",
                          sprintf('ggtitle(expression(atop("%s", atop(italic("%s")))))',
                                  ttl[1], ttl[2])),
                    "print(g)",
                    sep="\n")

  # Log the R command and execute.
  
  appendLog(Rtxt("Generate the dendrogram plot."), plot.cmd)
  newPlot()
  eval(parse(text=plot.cmd))

  # TODO 130311 How to identify the clusters in the plot, if
  # specified.

  ## nclust <- theWidget("hclust_clusters_spinbutton")$getValue()
  ## if (nclust > 1 && nclust <= length(crs$hclust$height))
  ## {
  ##   rect.cmd <- sprintf("rect.hclust(crs$hclust, k=%d)", nclust)
  ##   appendLog(Rtxt("Add in rectangles to show the clusters."), rect.cmd)
  ##   eval(parse(text=rect.cmd))
  ## }
}

plotDendrogram <- function()
{

  # Make sure there is a hclust object first.

  if (is.null(crs$hclust))
  {
    errorDialog(Rtxt("There is no Hierarchical Cluster yet we are",
                     "trying to plot it."), crv$support.msg)
    return()
  }

  # Load the required package into the library.  The library, cba,
  # should already be loaded. But check anyhow.

  lib.cmd <- "library(cba, quietly=TRUE)"
  if (! packageIsAvailable("cba", Rtxt("plot a dendrogram"))) return(FALSE)
  appendLog(packageProvides("cba", "plot"), lib.cmd)
  eval(parse(text=lib.cmd))

  # Show a busy cursor whilst drawing the plot.

  set.cursor("watch", Rtxt("Rendering the hierarchical cluster dendrogram...."))
  on.exit(set.cursor("left-ptr", ""))
  
  # Generate the plot command to not print the xaxis labels if there
  # are too many observations.

  if (length(crs$hclust$order) > 100)
    limit <- ", labels=FALSE, hang=0"
  else
    limit <- ""
  plot.cmd <- paste(sprintf('plot(crs$hclust, main="", sub="", xlab=""%s)\n',
                            limit),
                    genPlotTitleCmd(Rtxt("Cluster Dendrogram"), crs$dataname),
                    sep="")

  # Log the R command and execute.
  
  appendLog(Rtxt("Generate a dendrogram plot."), plot.cmd)
  newPlot()
  eval(parse(text=plot.cmd))

  # Identify the clusters in the plot, if specified.

  nclust <- theWidget("hclust_clusters_spinbutton")$getValue()
  if (nclust > 1 && nclust <= length(crs$hclust$height))
  {
    rect.cmd <- sprintf("rect.hclust(crs$hclust, k=%d)", nclust)
    appendLog(Rtxt("Add in rectangles to show the clusters."), rect.cmd)
    eval(parse(text=rect.cmd))
  }
}

displayHClustStats <- function()
{
  # Initial setup.
  
  TV <- "hclust_textview"

  # Make sure there is a cluster first.
  
  if (is.null(crs$hclust))
  {
    errorDialog(Rtxt("No cluster to plot.",
                     "The button should not have been sensitive."),
                crv$support.msg)
    return()
  }

  # The fpc package provides is available for cluster.stats function.
  
  if (!packageIsAvailable("fpc", Rtxt("calculate cluster statistics"))) return()
  lib.cmd <- "library(fpc, quietly=TRUE)"
  appendLog(packageProvides("fpc", "cluster.stats"), lib.cmd)
  eval(parse(text=lib.cmd))

  # 090323 Don't reset the textview since we want to reatin the build
  # information.

  # 090323 REMOVE resetTextview(TV)

  # Some background information.  Assume we have already built the
  # cluster, and so we don't need to check so many conditions.

  nclust <- theWidget("hclust_clusters_spinbutton")$getValue()
  sampling  <- not.null(crs$sample)
#  nums <- seq(1,ncol(crs$dataset))[as.logical(sapply(crs$dataset, is.numeric))]
#  if (length(nums) > 0)
#  {
#    indicies <- getVariableIndicies(crs$input)
#    include <- simplifyNumberList(intersect(nums, indicies))
#  }
#
#  if (length(nums) == 0 || length(indicies) == 0)
#  {
#    errorDialog("Clusters are currently calculated only for numeric data.",
#                "No numeric variables were found in the dataset",
#                "from amongst those having an input/target/risk role.")
#    return()
#  }

  include <- "crs$numeric" # 20110102 getNumericVariables()
  
  # Cluster centers.

  centers.cmd <- sprintf("centers.hclust(na.omit(crs$dataset[%s, %s]), crs$hclust, %d)",
                       ifelse(sampling, "crs$sample", ""), include, nclust)
  appendLog(Rtxt("List the suggested cluster centers for each cluster"), centers.cmd)
  appendTextview(TV, Rtxt("Cluster means:"), "\n\n",
                 collectOutput(centers.cmd, use.print=TRUE))
  
  # STATS: Log the R command and execute.

  stats.cmd <- sprintf(paste("cluster.stats(dist(na.omit(crs$dataset[%s, %s])),",
                             "cutree(crs$hclust, %d))\n"),
                       ifelse(sampling, "crs$sample", ""), include,
                       nclust)
  appendLog(Rtxt("Generate cluster statistics using the fpc package."), stats.cmd)
  appendTextview(TV, Rtxt("General cluster statistics:"), "\n\n",
                 collectOutput(stats.cmd, use.print=TRUE))

  setStatusBar(Rtxt("HClust cluster statistics have been generated."))
}

## THIS IS NOT EVEN RELATED TO hclust!!!! USES PAM

## on_hclust_seriation_button_clicked <- function(button)
## {

##   ## Make sure there is a hclust object first.

##   if (is.null(crs$hclust))
##   {
##     errorDialog("SHOULD NOT BE HERE.", crv$support.msg)
##     return()
##   }

##   ## The library, cba, should already be loaded. But check anyhow. I
##   ## think this is required for the seriation. Need to check.

##   lib.cmd <- "library(cba, quietly=TRUE)"
##   if (! packageIsAvailable("cba", "generate a seriation plot")) return()
##   appendLog(packageProvides("cba", "Seriation"), lib.cmd)
##   eval(parse(text=lib.cmd))
  
##   ## Some background information.

##   sampling  <- not.null(crs$sample)
##   nums <- seq(1,ncol(crs$dataset))[as.logical(sapply(crs$dataset, is.numeric))]
##   if (length(nums) > 0)
##   {
##     indicies <- getVariableIndicies(crs$input)
##     include <- simplifyNumberList(intersect(nums, indicies))
##   }

##   plot.cmd <- paste("d <- dist(as.matrix(crs$dataset",
##                     sprintf("[%s, %s]",
##                             ifelse(sampling, "crs$sample", ""),
##                             include),
##                     "))\n",
##                     "l <- pam(d, 10, cluster.only = TRUE)\n",
##                     "res <- cluproxplot(d, l, method = ",
##                     'c("Optimal", "Optimal"), plot = FALSE)\n',
##                     'plot(res, plotOptions = list(main = "PAM + ',
##                     'Seriation (Optimal Leaf ordering)", ',
##                     'col = terrain.colors(64)))', sep="")

##   appendLog("Generate a seriation plot.", plot.cmd)
##   newPlot()
##   eval(parse(text=plot.cmd))
  
##   setStatusBar("Seriation plot completed.")
## }

# 100424 Support resetting of widgets on loading a project.

showModelHClustExists <- function(state=!is.null(crs$hclust))
{
  # If a model exists then make available the appropriate widgets.
  
  theWidget("hclust_dendrogram_button")$setSensitive(TRUE)
  theWidget("hclust_clusters_label")$setSensitive(TRUE)
  theWidget("hclust_clusters_spinbutton")$setSensitive(TRUE)
  theWidget("hclust_stats_button")$setSensitive(TRUE)
  theWidget("hclust_data_plot_button")$setSensitive(TRUE)
  theWidget("hclust_discriminant_plot_button")$setSensitive(TRUE)

}

########################################################################
# EXPORT

exportHClustTab <- function()
{
  # Make sure we have a model first!
  
  if (noModelAvailable(crs$hclust, crv$HCLUST)) return(FALSE)

  # Get some required information

  sampling  <- not.null(crs$sample)
  nclust <- theWidget("hclust_clusters_spinbutton")$getValue()
  include <- "crs$numeric" # 20110102 getNumericVariables()
  
  startLog(paste(Rtxt("Export"), commonName(crv$HCLUST)))
  
  save.name <- getExportSaveName(crv$HCLUST)
  if (is.null(save.name)) return(FALSE)
  ext <- tolower(get.extension(save.name))

  # Construct the command to produce PMML.

  pmml.cmd <- sprintf(paste("pmml(crs$hclust, centers=centers.hclust(",
                           "na.omit(crs$dataset[%s, %s]), crs$hclust, %d)%s)",
                            sep=""),
                      ifelse(sampling, "crs$sample", ""), include, nclust,
                      ifelse(length(crs$transforms) > 0,
                             ", transforms=crs$transforms", ""))

  # We can't pass "\" in a filename to the parse command in
  # MS/Windows so we have to run the save/write command separately,
  # i.e., not inside the string that is being parsed.

  if (ext == "xml")
  {
    appendLog(Rtxt("Export hierarchical cluster as PMML."),
              sprintf('saveXML(%s, "%s")', pmml.cmd, save.name))
    XML::saveXML(eval(parse(text=pmml.cmd)), save.name)
  }
  else if (ext == "c")
  {
    if (isWindows()) save.name <- tolower(save.name)
    
    model.name <- sub("\\.c", "", basename(save.name))

    export.cmd <- generateExportPMMLtoC(model.name, save.name, "hclust_textview")
    
    appendLog(sprintf(Rtxt("Export %s as a C routine."), commonName(crv$HCLUST)),
              sprintf('pmml.cmd <- "%s"\n\n', pmml.cmd),
              export.cmd)

    eval(parse(text=export.cmd))

  }
  
  setStatusBar(sprintf(Rtxt("The model has been exported to '%s'."), save.name))

}

########################################################################
# SCORE

predict.hclust <- function(object, data, x, nclust=10, ...)
{
  # 090126 Initial work on a predict.hclust function, to allow using a
  # hclust model to allocate new DATA to pre-existing clusters that
  # are built from another dataset X. This uses the common model
  # interface function, predict. This makes it easy to use the Rattle
  # modelling code on kmeans. We use a kmeans encoding to generate the
  # clusters. This is only an approximation. Gets pretty close for
  # ward link and euclidean distance.

  object$centers <- centers.hclust(x, object, nclust=nclust, use.median=FALSE)
  rownames(object$centers) <- seq_len(nclust)
  return(predict.kmeans(object, data))
}

genPredictHclust <- function(dataset)
{
  # 081227 Generate a command to obtain the prediction results when
  # applying the model to new data.

  nclust <- theWidget("hclust_clusters_spinbutton")$getValue()
  sampling  <- not.null(crs$sample)
  include <- "crs$numeric" # 20110102 getNumericVariables()

  return(sprintf("crs$pr <- predict(crs$hclust, %s, na.omit(crs$dataset[%s, %s]), %s)",
                 dataset, ifelse(sampling, "crs$sample", ""), include, nclust))
}

genResponseHclust <- function(dataset)
{
  # 081227 Generate a command to obtain the response when applying the
  # model to new data.
  
  return(genPredictHclust(dataset))
}

genProbabilityHclust <- function(dataset)
{
  # 081227 Generate a command to obtain the probability when applying
  # the model to new data. There is probably a prblem with simply
  # using the cluster label as the output, since it won't be a
  # probability or even look like it. Let's do it for now though -
  # should be okay.
  
  return(genPredictHclust(dataset))
}
