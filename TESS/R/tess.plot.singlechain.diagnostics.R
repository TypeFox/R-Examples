################################################################################
#
# tess.plot.singlechain.diagnostics.R
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
# @brief Plotting the single-chain mcmc diagnostics of a episodic diversification rate analysis with mass-extinction events.
#
# @date Last modified: 2014-10-05
# @author Michael R May
# @version 2.0
# @since 2014-10-04, version 2.0.0
#
# @param    output            list          The processed output for plotting.
# @param    parameters        character     Which parameters to diagnose. See details for a complete description.
# @param    diagnostics       character     Which diagnostics to use. Options are "ESS" and "geweke".
# @param    ess.crit          numeric       Two values which correspond to low ESS threshold and acceptable ESS threshold. Default values are 100 and 200.
# @param    geweke.crit       numeric       The p-value cutoff for Geweke's diagnostic. Default is the canonical 0.05.
# @param    correction        character     What type of multiple-correction method to use. Options are "bonferroni" and "sidak".
# @param    xlab              character     The label of the x-axis. By default, millions of years.
# @param    col               character     Colors used for printing. Must be of same length as fig.types.
# @param    xaxt              character     The type of x-axis to plot. By default, no x-axis is plotted (recommended).
# @param    yaxt              character     The type of y-axis to plot.
# @param    pch               integer       The type of points to draw (if points are drawn).
# @param    ...                             Parameters delegated to various plotting functions.
#
#
################################################################################

tess.plot.singlechain.diagnostics = function(output,
                                      parameters=c("speciation rates",
                                                   "speciation shift times",
                                                   "extinction rates",
                                                   "extinction shift times",
                                                   "net-diversification rates",
                                                   "relative-extinction rates",
                                                   "mass extinction times"),
                                      diagnostics=c("ESS","geweke"),
                                      ess.crit=c(100,200),
                                      geweke.crit=0.05,
                                      correction="bonferroni",
                                      xlab="million years ago",
                                      col=NULL,
                                      xaxt="n",
                                      yaxt="s",
                                      pch=19,
                                      ...){

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
  validDiagnostics <- c("ESS","geweke")

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
  treeAge <- max(branching.times(output$tree))
  numIntervals <- length(output$intervals)-1
  plotAt <- 0:numIntervals
  intervalSize <- treeAge/numIntervals
  labels <- pretty(c(0,treeAge))
  labelsAt <- numIntervals - (labels / intervalSize)

  for ( type in parameters ) {
    for ( diag in diagnostics ) {

      if ( diag == "ESS" ) {

        thisOutput <- output[[type]]
        thisESS    <- effectiveSize(thisOutput)
        thisCol    <- col[findInterval(thisESS,ess.crit)+1]
        ylim       <- range(pretty(c(0,max(thisESS))))

        # For samples that had 0 ESS (always had the same parameter value),
        # reset the value to MAX and color the bar grey.
        thisCol[thisESS == 0] <- 'grey90'
        thisESS[thisESS == 0] <- max(ylim)

        barplot(thisESS,space=0,xaxt=xaxt,col=thisCol,border=NA,main=type,ylab="effective sample size",xlab=xlab,ylim=ylim,...)
        abline(h=ess.crit,lty=2,...)
        axis(1,at=labelsAt,labels=labels)
        box()

      } else if ( diag == "geweke" ) {

        thisOutput <- output[[type]]
        thisGeweke <- geweke.diag(thisOutput)$z

        # Compute the p-values
        if ( !is.null(correction) ) {
          if ( correction == "bonferroni" ) {
            crit <- geweke.crit / numIntervals
          }
          if ( correction == "sidak" ) {
            crit <- 1 - (1 - geweke.crit)^(1/numIntervals)
          }
        } else {
          crit <- geweke.crit
        }

        thisPvalue <- pnorm(thisGeweke)
        thisPvalue[is.na(thisPvalue)] <- 0
        failed <- thisPvalue < crit | thisPvalue > 1 - (crit)
        thisCol <- ifelse(failed,col[1],col[3])
        criticalGewekeValues <- qnorm(c(crit,1-crit))
        ylim <- range(pretty(c(thisGeweke,criticalGewekeValues)),finite=TRUE)

        plot(thisGeweke,col=thisCol,type="p",xaxt=xaxt,ylab="Geweke statistic",main=type,xlab=xlab,ylim=ylim,xlim=range(plotAt),pch=pch,...)
        abline(h=criticalGewekeValues,lty=2,...)
        axis(1,at=labelsAt,labels=labels)

      }

    }

  }

}
