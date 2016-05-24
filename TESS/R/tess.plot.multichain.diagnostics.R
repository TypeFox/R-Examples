################################################################################
#
# tess.plot.multichain.diagnostics.R
#
# Copyright (c) 2012- Michael R May
#
# This file is part of TESS.
# See the NOTICE file distributed with this work for additional
# information regarding copyright ownership and licensing.
#
# TESS is free software; you can redistribute it and/or modify
# it under the terms of the GNU Lesser General Public License as
# published by the Free Software Foundation; either version 2
# of the License, or (at your option) any later version.
#
#  TESS is distributed in the hope that it will be useful,
#  but WITHOUT ANY WARRANTY; without even the implied warranty of
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#  GNU Lesser General Public License for more details.
#
# You should have received a copy of the GNU Lesser General Public
# License along with TESS; if not, write to the
# Free Software Foundation, Inc., 51 Franklin St, Fifth Floor,
# Boston, MA  02110-1301  USA
#
################################################################################

################################################################################
#
# @brief Plotting the multi-chain mcmc diagnostics of a episodic diversification rate analysis with mass-extinction events.
#
# @date Last modified: 2014-10-05
# @author Michael R May
# @version 2.0
# @since 2014-10-04, version 2.0.0
#
# @param    output            list          The processed output for plotting.
# @param    parameters        character     Which parameters to diagnose. See details for a complete description.
# @param    diagnostics       character     Which diagnostics to use. Options are "ESS" and "geweke".
# @param    gelman.crit       numeric       The critical value above which a Rubin-Gelman statistic is considered a failure. Default is 1.5
# @param    xlab              character     The label of the x-axis. By default, millions of years.
# @param    col               character     Colors used for printing. Must be of same length as fig.types.
# @param    xaxt              character     The type of x-axis to plot. By default, no x-axis is plotted (recommended).
# @param    yaxt              character     The type of y-axis to plot.
# @param    pch               integer       The type of points to draw (if points are drawn).
# @param    ...                             Parameters delegated to various plotting functions.
#
#
################################################################################

tess.plot.multichain.diagnostics = function(outputs,
                                            parameters=c("speciation rates",
                                                         "speciation shift times",
                                                         "extinction rates",
                                                         "extinction shift times",
                                                         "net-diversification rates",
                                                         "relative-extinction rates",
                                                         "mass extinction times"),
                                            diagnostics="Gelman-Rubin",
                                            gelman.crit=1.05,
                                            xlab="million years ago",
                                            col=NULL,
                                            xaxt="n",
                                            yaxt="s",
                                            pch=19,
                                            ...) {

  # Check that parameter type is valid
  validFigTypes <- c("speciation rates",
                     "speciation shift times",
                     "extinction rates",
                     "extinction shift times",
                     "net-diversification rates",
                     "relative-extinction rates",
                     "mass extinction times")

  invalidFigTypes <- parameters[!parameters %in% validFigTypes]

  if ( length( invalidFigTypes ) > 0 ) {
    stop("\nThe following figure types are invalid: ",paste(invalidFigTypes,collapse=", "),".",
         "\nValid options are: ",paste(validFigTypes,collapse=", "),".")
  }

  # Check that the diagnostics are valid
  validDiagnostics <- c("Gelman-Rubin")

  invalidDiagnostics <- diagnostics[!diagnostics %in% validDiagnostics]

  if ( length( invalidFigTypes ) > 0 ) {
    stop("\nThe following diagnostics are invalid: ",paste(invalidDiagnostics,collapse=", "),".",
         "\nValid options are: ",paste(validDiagnostics,collapse=", "),".")
  }

  # Make color vector
  if ( is.null(col) ) {
    col <- c("#E41A1C","#FF7F00","#377EB8")
  } else if ( length(col) != 3 ){
    stop("\nYou must supply 3 input colors.")
  }

  # Compute the axes
  treeAge <- max(branching.times(outputs[[1]]$tree))
  numIntervals <- length(outputs[[1]]$intervals)-1
  plotAt <- 0:numIntervals 
  intervalSize <- treeAge/numIntervals
  labels <- pretty(c(0,treeAge))
  labelsAt <- numIntervals - (labels / intervalSize)

  for( type in parameters ) {
    for( diag in diagnostics ) {
      if ( diag == "Gelman-Rubin" ){

        # Get the outputs for this parameter
        this_output_list <- lapply(outputs,function(x) x[[type]]  )

        # Get the chain with the least number of samples and trim all of the chains to be of this length
        min_num_samples <- min(sapply(this_output_list,nrow))
        this_output_list_trimmed <- lapply(this_output_list,function(output){
          as.mcmc(output[(nrow(output) - min_num_samples + 1):nrow(output),])
        })

        # Find the parameters with variance 0
        this_output_variance <- do.call(rbind,lapply(this_output_list_trimmed,function(output){
          apply(output,2,var)
        }))
        outputs_with_no_variance <- apply(this_output_variance == 0,2,any)

        # Make the container for the PSRF
        psrf <- numeric(length(outputs_with_no_variance))

        # Compute the PSRF for each parameter
        for(i in 1:length(outputs_with_no_variance)){
          if( !outputs_with_no_variance[i] ) {
            this_list <- as.mcmc.list(lapply(this_output_list_trimmed,function(x) x[,i]))
            psrf[i] <- gelman.diag(this_list)$psrf[1,1]
          } else {
            psrf[i] <- NA
          }
        }

        if ( sum( is.finite(psrf) ) > 0 ) {
          # Compute the colors
          failed <- abs(1 - psrf) > abs(1 - gelman.crit)
          thisCol <- ifelse(failed,col[1],col[3])

          # Plot the PSRF
          plot(psrf,col=thisCol,type="p",xaxt=xaxt,ylab="Rubin-Gelman statistic",main=type,xlab=xlab,xlim=range(plotAt),pch=pch,...)
          abline(h=c(1,gelman.crit),lty=2,...)
          axis(1,at=labelsAt,labels=labels)
        }
      }
    }
  }
}