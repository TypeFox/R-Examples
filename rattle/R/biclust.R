# Gnome R Data Miner: GNOME interface to R for Data Mining
#
# Time-stamp: <2015-05-17 08:54:15 gjw>
#
# Implement biclust functionality.
#
# Copyright (c) 2010 Togaware Pty Ltd
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
# ToDo 100121
#
# Graphical display of output.
# Allow choice of methods and options

########################################################################
# Callbacks

# When a radio button is selected, display the appropriate tab page.

on_biclust_radiobutton_toggled <- function(button)
{
  if (button$getActive())
    crv$CLUSTER$setCurrentPage(crv$CLUSTER.BICLUST.TAB)
  setStatusBar()
}

########################################################################
# Execution

executeClusterBiclust <- function(include)
{
  TV <- "biclust_textview"
  sampling  <- not.null(crs$sample)

  # Obtain interface information.

  method <- "BCCC"
  seed <- "crv$seed"

  # Start the log.
  
  startLog(commonName(crv$BICLUST))

  # Load the required package.
  
  lib.cmd <- "library(biclust, quietly=TRUE)"
  if (! packageIsAvailable("biclust", Rtxt("perform bicluster analysis"))) return(FALSE)
  appendLog(packageProvides('biclust', 'biclust'), lib.cmd)
  eval(parse(text=lib.cmd))

  # Set the seed so we can repeat.

  seed.cmd <- sprintf('set.seed(%s)', seed)
  appendLog(Rtxt("Reset the random number seed to obtain the same results each time."),
            seed.cmd)
  eval(parse(text=seed.cmd))

  # Build the model.
  
  biclust.cmd <- sprintf(paste('crs$biclust <- biclust(',
                               'as.matrix(na.omit(crs$dataset[%s, %s])),',
                               'method=%s)', sep=""),
                         ifelse(sampling, "crs$sample", ""),
                         include, method)

  appendLog(sprintf(Rtxt("Generate %s using method '%s'."),
                    commonName(crv$BICLUST), method),
            biclust.cmd)

  start.time <- Sys.time()

  result <- try(eval(parse(text=biclust.cmd)), TRUE)
  time.taken <- Sys.time()-start.time

  # Check for errors.
  
  if (inherits(result, "try-error"))
  {
    errorDialog(errorMessageFun("biclust", result))
    return(FALSE)
  }

  # Show the results.

  print.cmd <- "print(crs$biclust)"
  
  appendLog(sprintf(Rtxt("Generate a textual view of the %s model."),
                    commonName(crv$BICLUST)),
            print.cmd)
  
  resetTextview(TV)
  setTextview(TV,
              sprintf(Rtxt("Summary of the %s model (built using '%s'):"),
                      commonName(crv$BICLUST), "biclust"),
              "\n",
              collectOutput(print.cmd))

  reportTimeTaken(TV, time.taken, model=commonName(crv$BICLUST))

  return(TRUE)
}

