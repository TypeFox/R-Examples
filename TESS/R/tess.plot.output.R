################################################################################
#
# tess.plot.output.R
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
# @brief Plotting the output of a episodic diversification rate analysis with mass-extinction events.
#
# @date Last modified: 2014-10-05
# @author Michael R May
# @version 2.0
# @since 2014-10-04, version 2.0.0
#
# @param    output            list          The processed output for plotting.
# @param    fig.types         character     Which aspects of the model to visualize. See details for a complete description.
# @param    xlab              character     The label of the x-axis. By default, millions of years.
# @param    col               character     Colors used for printing. Must be of same length as fig.types.
# @param    col.alpha         numeric       Alpha channel parameter for credible intervals.
# @param    xaxt              character     The type of x-axis to plot. By default, no x-axis is plotted (recommended).
# @param    yaxt              character     The type of y-axis to plot.
# @param    pch               integer       The type of points to draw (if points are drawn).
# @param    ...                             Parameters delegated to various plotting functions.
#
#
################################################################################

tess.plot.output = function(output,fig.types=c("speciation rates","speciation shift times","speciation Bayes factors",
                                               "extinction rates","extinction shift times","extinction Bayes factors",
                                               "net-diversification rates","relative-extinction rates",
                                               "mass extinction times","mass extinction Bayes factors"),
                            xlab="million years ago",col=NULL,col.alpha=50,
                            xaxt="n",yaxt="s",pch=19,plot.tree=FALSE,
                            ...){

  # Check that fig type is valid
  validFigTypes <- c("speciation rates","speciation shift times","speciation Bayes factors",
                     "extinction rates","extinction shift times","extinction Bayes factors",
                     "net-diversification rates","relative-extinction rates",
                     "mass extinction times","mass extinction Bayes factors")
  invalidFigTypes <- fig.types[!fig.types %in% validFigTypes]

  if ( length( invalidFigTypes ) > 0 ) {
    stop("\nThe following figure types are invalid: ",paste(invalidFigTypes,collapse=", "),".",
         "\nValid options are: ",paste(validFigTypes,collapse=", "),".")
  }

  # Make color vector
  if ( is.null(col) ) {
    col <- c("speciation rates"="#984EA3",
             "speciation shift times"="#984EA3",
             "speciation Bayes factors"="#984EA3",
             "extinction rates"="#E41A1C",
             "extinction shift times"="#E41A1C",
             "extinction Bayes factors"="#E41A1C",
             "net-diversification rates"="#377EB8",
             "relative-extinction rates"="#FF7F00",
             "mass extinction times"="#4DAF4A",
             "mass extinction Bayes factors"="#4DAF4A")
  } else {
    names(col) <- fig.types
  }

  # Compute the axes
  treeAge <- max(branching.times(output$tree))
  numIntervals <- length(output$intervals)-1
  plotAt <- 0:numIntervals
  intervalSize <- treeAge/numIntervals
  labels <- pretty(c(0,treeAge))
  labelsAt <- numIntervals - (labels / intervalSize)

  for( type in fig.types ) {

    if ( grepl("times",type) ) {

      thisOutput <- output[[type]]
      meanThisOutput <- colMeans(thisOutput)
      criticalPP <- output[[grep(strsplit(type," ")[[1]][1],grep("CriticalPosteriorProbabilities",names(output),value=TRUE),value=TRUE)]]

      if(plot.tree){
        plot(output$tree,show.tip.label=FALSE,edge.col=rgb(0,0,0,0.10),x.lim=c(0,treeAge))
        par(new=TRUE)
      }

      barplot(meanThisOutput,space=0,xaxt=xaxt,col=col[type],border=col[type],main=type,ylab="posterior probability",xlab=xlab,ylim=c(0,1),...)
      abline(h=criticalPP,lty=2,...)
      axis(4,at=criticalPP,labels=2*log(output$criticalBayesFactors),las=1,tick=FALSE,line=-0.5)
      axis(1,at=labelsAt,labels=labels)
      box()

    } else if ( grepl("Bayes factors",type) ) {

      thisOutput <- output[[type]]
      ylim <- range(c(thisOutput,-10,10),finite=TRUE)

      if(plot.tree){
        plot(output$tree,show.tip.label=FALSE,edge.col=rgb(0,0,0,0.10),x.lim=c(0,treeAge))
        par(new=TRUE)
      }
      plot(x=plotAt[-1]-diff(plotAt[1:2])/2,y=thisOutput,type="p",xaxt=xaxt,col=col[type],ylab="Bayes factors",main=type,xlab=xlab,ylim=ylim,xlim=range(plotAt),pch=pch,...)
      abline(h=2 * log(output$criticalBayesFactors),lty=2,...)
      axis(4,at=2 * log(output$criticalBayesFactors),las=1,tick=FALSE,line=-0.5)
      axis(1,at=labelsAt,labels=labels)

    } else {

      thisOutput <- output[[type]]
      meanThisOutput <- colMeans(thisOutput)
      quantilesThisOutput <- apply(thisOutput,2,quantile,prob=c(0.025,0.975))
      if( type %in% c("speciation rates","extinction rates")){
        quantilesSpeciation <- apply(output[["speciation rates"]],2,quantile,prob=c(0.025,0.975))
        quantilesExtinction <- apply(output[["extinction rates"]],2,quantile,prob=c(0.025,0.975))
        ylim <- c(0,max(quantilesSpeciation,quantilesExtinction))
      } else {
        ylim <- c(0,max(quantilesThisOutput))
      }

      if(plot.tree){
        plot(output$tree,show.tip.label=FALSE,edge.col=rgb(0,0,0,0.10),x.lim=c(0,treeAge))
        par(new=TRUE)
      }
      plot(x=plotAt,y=c(meanThisOutput[1],meanThisOutput),type="l",ylim=ylim,xaxt=xaxt,col=col[type],ylab="rate",main=type,xlab=xlab,...)
      polygon(x=c(0:ncol(quantilesThisOutput),ncol(quantilesThisOutput):0),y=c(c(quantilesThisOutput[1,1],quantilesThisOutput[1,]),rev(c(quantilesThisOutput[2,1],quantilesThisOutput[2,]))),border=NA,col=paste(col[type],col.alpha,sep=""))
      axis(1,at=labelsAt,labels=labels)

    }

  }

}