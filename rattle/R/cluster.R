# Gnome R Data Miner: GNOME interface to R for Data Mining
#
# Time-stamp: <2014-07-18 15:08:01 gjw>
#
# Implement cluster functionality.
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
# EXECUTION

executeClusterTab <- function()
{
  # Can not cluster without a dataset.

  if (noDatasetLoaded()) return()

  # If it looks like the VARIABLES page has not been executed, complain..

  if (variablesHaveChanged(Rtxt("building clusters"))) return()

  # Check if sampling needs executing.

  if (sampleNeedsExecute()) return()

  # 091216 Automatically handle any selected categorics by converting
  # them to numeric, so they then become included variables. This
  # works, but it risks suprising the user with hte addition of new
  # variables outside their control. So let's leave it to the user to
  # do the transforms, or use clara.

  # factors <- crs$input[sapply(crs$input, function(x)
  #                             is.factor(crs$dataset[[x]]))]
  # sapply(factors, executeTransformRemapPerform, action="indicator",
  #        remap.prefix="TIN")
  
  # Kmeans and hclust only work for numeric data, so identify
  # variables to include.  Only work with the INPUT/TARGET/RISK
  # variables. That is, only exclude the IGNORE and IDENT variables.

  include <- "crs$numeric" # 20110102 getNumericVariables()
  if (! length(include))
  {
    errorDialog(Rtxt("Clusters are currently calculated only for numeric data.",
                     "No numeric variables were found in the dataset",
                     "from amongst those having an input/target/risk role."))
    return()
  }

  # Dispatch.

  if (theWidget("kmeans_radiobutton")$getActive())
  {
    if (executeClusterKMeans(include))
      theWidget("evaluate_kmeans_checkbutton")$setActive(TRUE)
  }
  else if (theWidget("ewkm_radiobutton")$getActive())
  {
    if (executeClusterEwkm(include))
      theWidget("evaluate_kmeans_checkbutton")$setActive(TRUE)
  }
  ## else if (theWidget("clara_radiobutton")$getActive())
  ## {
  ##   infoDialog(Rtxt("Not yet implemented."))
  ##   if (executeClusterClara(include))
  ##     theWidget("evaluate_clara_checkbutton")$setActive(TRUE)
  ## }
  ## else if (theWidget("pam_radiobutton")$getActive())
  ## {
  ##   infoDialog(Rtxt("Not yet implemented."))
  ##   if (executeClusterPam(include))
  ##     theWidget("evaluate_pam_checkbutton")$setActive(TRUE)
  ## }
  else if (theWidget("hclust_radiobutton")$getActive())
  {
    if (executeClusterHClust(include))
      theWidget("evaluate_hclust_checkbutton")$setActive(TRUE)
  }
  else if (theWidget("biclust_radiobutton")$getActive())
  {
    executeClusterBiclust(include)
#      theWidget("evaluate_biclust_checkbutton")$setActive(TRUE)
  }
}

########################################################################
# EXPORT

exportClusterTab <- function()
{
  
  if (noDatasetLoaded()) return()

  if (theWidget("kmeans_radiobutton")$getActive())
  {
    exportKMeansTab()
  }
  else if (theWidget("ewkm_radiobutton")$getActive())
  {
    exportEwkmTab()
  }
  else if (theWidget("hclust_radiobutton")$getActive())
  {
    exportHClustTab()
  }
  else
  {
    errorDialog(Rtxt("PMML export for this model is not yet implemented."))
    return()
  }
}
