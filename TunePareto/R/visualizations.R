#################################################################
# Visualization routines
#################################################################

# Retrieves a list of Pareto fronts from 
# a domination matrix <mat>.
#
# Returns a list of vectors and a rank vector
calculateParetoFronts <- function(mat)
{
  ranks <- rep(NA,nrow(mat))
  paretoFronts <- list()
  indexList <- 1:nrow(mat)
  
  paretoFront <- 1
  # calculate the Pareto fronts
  while (length(indexList) != 0)
  {
    nonDominated <- which(!apply(mat,2,any))
    
    ranks[indexList[nonDominated]] <- paretoFront
    paretoFronts[[paretoFront]] <- indexList[nonDominated]
   
    # remove the Pareto front from the list of points to consider 
    indexList <- indexList[-nonDominated]
    mat <- mat[-nonDominated, -nonDominated, drop=FALSE]
    paretoFront <- paretoFront + 1
  }
  return(list(paretoFronts=paretoFronts, ranks=ranks))
}

# Plots a graph in which each column of nodes represents one Pareto front
# and edges represent a domination relation.
# <tuneParetoResult> is a TuneParetoResult object to be plotted.
# If <transitiveReduction> is true, transitive edges are removed from the graph.
# If <drawDominatedObjectives> is true, color bars that indicate the objectives in
# which a configuration is the best in its Pareto front are drawn.
# <drawLabels> specifies whether the configuration descriptions should be drawn.
# <drawLegend> specifies whether a legend for the color bars should be plotted.
# <col.indicator> is a vector of colors for the objectives in the indicators
# <pch.indicator> is a vector of symbols for the objectives in the indicators
# <cex.indicator> is a value or vector specifying the sizes of the symbols in the indicators
# <x.legend> specifies the position of this legend.
# <cex.legend> specifies the character size of the legend
#
# Invisibly returns the igraph object.
plotDominationGraph <- function(tuneParetoResult, 
                                transitiveReduction=TRUE, 
                                drawDominatedObjectives=TRUE,
                                drawLabels=TRUE,
                                drawLegend=TRUE,
                                x.legend="topleft", 
                                cex.legend=0.7,
                                col.indicator, 
                                pch.indicator=15, 
                                cex.indicator=0.8,
                                 ...)
{
  if(!inherits(tuneParetoResult, "TuneParetoResult"))
    stop("\"tuneParetoResult\" must be a TuneParetoResult object!")
  
  require(igraph)
  
  if (missing(col.indicator))
      colorSet <- c("blue","green","red","darkgoldenrod","gold","brown","cyan",
        "purple","orange","seagreen","tomato","darkgray","chocolate",
        "maroon","darkgreen","gray12","blue4","cadetblue","darkgoldenrod4",
        "burlywood2")
  else
  {
    if (length(col.indicator) != length(tuneParetoResult$minimizeObjectives))
      stop("Please specify one color for each objective function in col.indicator!")
    colorSet <- col.indicator
  }
  mat <- tuneParetoResult$dominationMatrix
  edges <- tuneParetoResult$dominationMatrix
  
  if (transitiveReduction)
  # remove transitive edges
  {
    transitiveGraph <- edges
    # first, build the complete transitive graph
    for (i in 1:nrow(edges))
    {
      for (j in 1:nrow(edges))
      {
        for (k in 1:nrow(edges))
        {
          if (edges[i,j] && edges[j,k])
            transitiveGraph[i,k] <- TRUE
        }
      }
    }

    # now check if there are transitive dependencies
    # and remove the corresponding edges
    for (i in 1:nrow(transitiveGraph))
    {
      for (j in 1:nrow(transitiveGraph))
      {
        for (k in 1:nrow(transitiveGraph))
        {
          if (transitiveGraph[i,j] && transitiveGraph[j,k])
            edges[i,k] <- FALSE
        }
      }
    }

  }
  
  paretoFronts <- calculateParetoFronts(mat)$paretoFronts
  
  # create igraph object
  g <- graph.adjacency(t(edges), mode="directed")
  
  dim_y <- 10
  
  positions <- matrix(ncol=2, nrow=nrow(edges))
  
  # calculate graph layout
  for (i in seq_along(paretoFronts))
  {

    indices <- paretoFronts[[i]]
    
    if (i %% 2 == 0)
    {
      distance <- dim_y/length(indices)
      offset <- 0
    }
    else
    {
      distance <- (dim_y - 1)/length(indices)
      offset <- 0.5
    }  
    positions[indices, ] <- t(sapply(seq_along(indices),function(j)
                              {
                                c(i + 1, offset + (j-1)*distance+distance/2)
                              }))
  }
  
  args <- list(...)
  
  # check for certain graphical parameters in ... 
  # that have different default values in this plot
  if (is.null(args$vertex.size))
    args$vertex.size <- 2
    
  if (is.null(args$vertex.color))
    args$vertex.color <- "grey"
    
  if (is.null(args$vertex.label.cex))
    args$vertex.label.cex <- 0.7
    
  if (is.null(args$vertex.label.dist))
    args$vertex.label.dist <- 0.2
    
  if (is.null(args$edge.arrow.size))
    args$edge.arrow.size <- 0.3
    
  if (is.null(args$edge.curved))
    args$edge.curved <- !transitiveReduction
    
  if(!is.null(args$legend.x))
  {
  # backward compatibility: 
  # interpret old parameter legend.x which has been renamed to x.legend
    x.legend <- args$legend.x
    args <- args[names(args) != "legend.x"]
  }
    
  if (drawLabels)
    args$vertex.label <- rownames(edges)
  else
    args$vertex.label <- NA
  
  do.call(plot,c(list(x=g),args,list(layout=positions[1:nrow(edges),,drop=FALSE])))
    
  # normalize layout
  positions <- layout.norm(positions,-1, 1, -1, 1)

  cols <- c()
  pchs <- c()
  cexs <- c()
  if (drawDominatedObjectives)
  {
  
    cin <- par("cin")
    Cex <- cex.indicator * par("cex") * 0.7

    pchWidth <- Cex * xinch(cin[1L], warn.log = FALSE)
  
    # calculate the dominated objectives for the color indicators
    dominatedObjectives <- lapply(paretoFronts,function(front)
    {
     r <- lapply(1:ncol(tuneParetoResult$testedObjectiveValues),function(j)
            {
              vec <- tuneParetoResult$testedObjectiveValues[front,j]
              if (tuneParetoResult$minimizeObjectives[j])
                front[which(vec == min(vec))]
              else
                front[which(vec == max(vec))]
            })
      matrix(sapply(front,function(x)
            {
              sapply(r,function(y)
              {
                x %in% y
              })                
            }), ncol=length(front))
    })
        
    # calculate the positions and colors of the color indicators      
    for (i in seq_along(paretoFronts))
    {
      for (c in ncol(dominatedObjectives[[i]]):1)
      {
        #sumDom <- 1          
        sumDom <- xinch(cin[1L], warn.log=FALSE)
        for (r in 1:nrow(dominatedObjectives[[i]]))
        {       
          if (dominatedObjectives[[i]][r,c])
          {
            r_inv <- nrow(dominatedObjectives[[i]]) - r
            positions <- rbind(positions,
                               c(positions[paretoFronts[[i]][c],1] - sumDom,#(sumDom+1)*pchWidth[(r - 1) %% length(pchWidth) + 1], 
                                 positions[paretoFronts[[i]][c],2]))
            cols <- c(cols, colorSet[(r_inv) %% length(colorSet) + 1])
            pchs <- c(pchs, pch.indicator[(r_inv) %% length(pch.indicator) + 1])
            cexs <- c(cexs, cex.indicator[(r_inv) %% length(cex.indicator) + 1])
            #sumDom <- sumDom + 1
            sumDom <- sumDom + pchWidth[(r_inv) %% length(pchWidth) + 1]
          }         
        }
      }
    }
    
    # draw the indicators     
    pts <- positions[-(1:nrow(edges)),,drop=FALSE]
    points(pts[,1], pts[,2], pch=pchs, col=cols,
           cex=cexs)
           
    if (drawLegend)
      # draw the legend       
      legend(x=x.legend, pch=pch.indicator, ncol=1, cex=cex.legend,
             col=colorSet[(1:ncol(tuneParetoResult$testedObjectiveValues) - 1) %% length(colorSet) + 1],
             legend = paste("Dominates front in",colnames(tuneParetoResult$testedObjectiveValues)),
             pt.cex=cex.indicator[(1:ncol(tuneParetoResult$testedObjectiveValues) - 1) %% length(cex.indicator) + 1], box.lty=0)
  } 
     
  # return the graph
  return(invisible(g))
}

# Internal function called from plotParetoFronts2D and plotObjectivePairs
# to plot the Pareto fronts of two objectives.
# <objectiveVals> is a matrix of objective values.
# <minimizeObjectives> is a Boolean vector specifying which objectives are minimized.
# <boundaries> is a vector of boundaries for the objectives, or NULL.
# If <drawLabels> is TRUE, the parameter configurations are printed.
# If <drawBoundaries is true, the boundaries of the objectives are plotted as lines.
# If <plotNew> is true, a new plot is started, otherwise lines are added to existing plots.
# If <fitLabels> is TRUE, overlaps of the parameter configuration labels are removed.
# <xlim>, <ylim>, <xlab>, <ylab>, <xaxs>, <yaxs>,
# <type> and <lwd> correspond to the parameters in the generic plot routine.
# <labelPos> specifies the position of the labels with respect to the points
# <cex.conf> is the text size of the configuration labels.
# <lty.fronts>, <pch.fronts> amd <col.fronts> specify the line types, the characters,
# and the colors for the Pareto fronts.
internal.plotParetoFronts2D <- function(objectiveVals, minimizeObjectives, boundaries,
                                        drawLabels=TRUE, drawBoundaries=TRUE, 
                                        plotNew=TRUE, 
                                        fitLabels=TRUE,
                                        xlim, ylim, xlab, ylab, xaxs="r", yaxs="r",
                                        labelPos=4, 
                                        type="o", lwd=1, cex.conf=0.5,
                                        lty.fronts=1, pch.fronts=8, col.fronts, ...)
{
  
  # calculate domination matrix of the 2 chosen objectives
  domination <- calculateDominationMatrix(objectiveVals, minimizeObjectives)
  
  # calculate the Pareto fronts
  paretoFronts <- calculateParetoFronts(domination)
  
  
  if (missing(xlim))
    xlim <- c(min(objectiveVals[,1]), max(objectiveVals[,1]))
  
  if (missing(ylim))
    ylim <- c(min(objectiveVals[,2]), max(objectiveVals[,2]))
    
  if (missing(xlab))
    xlab <- colnames(objectiveVals)[1]
  
  if (missing(ylab))
    ylab <- colnames(objectiveVals)[2]
    
  if (missing(col.fronts) || is.null(col.fronts))
    colorSet <- c("blue","green","red","darkgoldenrod","gold","brown","cyan",
      "purple","orange","seagreen","tomato","darkgray","chocolate",
      "maroon","darkgreen","gray12","blue4","cadetblue","darkgoldenrod4",
      "burlywood2")
  else
    colorSet <- col.fronts
  
    
  if (plotNew)
  {
        plot(1, type="n", 
           xlim=xlim,
           ylim=ylim,
           xlab=xlab,
           ylab=ylab,
           xaxs=xaxs,
           yaxs=yaxs, 
           lwd=lwd, ...)
  }
  
  if (drawBoundaries && !is.null(boundaries))
  {
    if (!is.na(boundaries[1]) && !is.null(boundaries[1]))
      abline(v=boundaries[1], col="darkgrey", lty=2)
    if (!is.na(boundaries[2]) && !is.null(boundaries[2]))
      abline(h=boundaries[2], col="darkgrey", lty=2)
  }
  
  pl <- c()
  cl <- c()
  for (i in seq_along(paretoFronts$paretoFronts))
  {
    # extract points and order them by their x value
    pts <- objectiveVals[paretoFronts$paretoFronts[[i]],,drop=FALSE]
    pts <- pts[order(pts[,1]),,drop=FALSE]
    
    pl <- rbind(pl,pts)
    cl <- c(cl, rep(colorSet[(i - 1) %% length(colorSet) + 1],nrow(pts)))
    
    # add new Pareto front
    lines(pts[,1], pts[,2], col=colorSet[(i - 1) %% length(colorSet) + 1], 
          pch=pch.fronts[(i - 1) %% length(pch.fronts) + 1], 
          type=type, lty=lty.fronts[(i - 1) %% length(lty.fronts) + 1], lwd=lwd)
 
  }
  
  if (drawLabels)
  {
    # add the configuration descriptions
    if (!fitLabels)
      # draw all labels
      remaining <- 1:nrow(pl)
    else
      # remove overlapping labels
      remaining <- removeOverlaps(labels=rownames(pl), 
                                  positions=pl,
                                  labelPos=labelPos, 
                                  cex=cex.conf, 
                                  xlim=xlim, 
                                  ylim=ylim,
                                  xaxs=xaxs,
                                  yaxs=yaxs)
  
    text(pl[remaining,1], pl[remaining,2], rownames(pl)[remaining], 
         cex=cex.conf, pos=labelPos, col=cl[remaining])
  }
}

removeOverlaps <- function(labels, positions, labelPos=4, cex=0.5, offset=0.5, xlim, ylim, xaxs, yaxs)
{
  cxy <- par("cxy")
  height <- strheight(labels,cex=cex)
  width <- strwidth(labels,cex=cex)
  
  if (xaxs == "r")
  # extend x range by 4%
  {
    xlimRangeExt <- (xlim[2] - xlim[1]) * 0.02
    xlim[1] <- xlim[1] - xlimRangeExt
    xlim[2] <- xlim[2] + xlimRangeExt
  }
  
  if (yaxs == "r")
  # extend y range by 4%
  {
    ylimRangeExt <- (ylim[2] - ylim[1]) * 0.02
    ylim[1] <- ylim[1] - ylimRangeExt
    ylim[2] <- ylim[2] + ylimRangeExt
  }
  
  # calculate positions of labels
  rects <- switch(labelPos,
    "1" = {
            cbind(positions[,1] - width/2,
                  positions[,2] - height - cxy[2]*offset,
                  positions[,1] + width/2,
                  positions[,2] - cxy[2]*offset)
          },
    "2" = {
            cbind(positions[,1] - width - cxy[1]*offset,
                  positions[,2] - height/2,
                  positions[,1] - cxy[1]*offset,
                  positions[,2] + height/2)          
          },      
    "3" = {
            cbind(positions[,1] - width/2,
                  positions[,2] + cxy[2]*offset,
                  positions[,1] + width/2,
                  positions[,2] + height + cxy[2]*offset)
        },
    "4" = {
          cbind(positions[,1] + cxy[1]*offset,
                positions[,2] - height/2,
                positions[,1] + width + cxy[1]*offset,
                positions[,2] + height/2)   
          })
  
  # determine which labels do not exceed the plotting region
  remaining <- which(apply(rects,1,function(pos)
                    {
                      pos[1] >= xlim[1] && pos[2] >= ylim[1] &&
                      pos[3] <= xlim[2] && pos[4] <= ylim[2]
                    }))
  
  # determine which labels overlap
  overlap <- matrix(apply(rects[remaining,,drop=FALSE],1,function(pos1)
                     apply(rects[remaining,,drop=FALSE],1,function(pos2)
                     {
                        pos1[1] <= pos2[3] && pos1[3] >= pos2[1] &&
                        pos1[2] <= pos2[4] && pos1[4] >= pos2[4]
                     })), ncol=length(remaining))
               
  diag(overlap) <- FALSE
  
  # calculate the number of overlaps for each label
  sums <- apply(overlap,2,sum)
  
  # apply a greedy strategy by iteratively
  # removing the labels with the most overlaps
  while(any(sums > 0))
  {
    # choose between equal overlap count randomly
    removeIdx <- sample(which(sums==max(sums)),size=1)
    
    # reduce overlap matrix and recalculate sum
    overlap <- overlap[-removeIdx,-removeIdx,drop=FALSE]
    sums <- apply(overlap,2,sum)
    
    # remove the label
    remaining <- remaining[-removeIdx]
  }
  
  return(remaining)               
}

# Plots a 2D plot of two objectives in an optimization.
# <tuneParetoResult> is an object of class TuneParetoResult.
# <objectives> is a vector of indices or names of the objectives to plot.
# If <drawLabels> is true, the descriptions of the configurations are written next to the plot
# If <drawBoundaries is true, the boundaries of the objectives are plotted as lines.
# <labelPos> specifies the position of the labels with respect to the points
# If <fitLabels> is TRUE, overlaps of the parameter configuration labels are removed.
# <cex.conf> is the text size of the configuration labels.
# <labelPos> is the position of the text relative to the points.
# <lty.fronts>, <pch.fronts> amd <col.fronts> specify the line types, the characters,
# and the colors for the Pareto fronts.
plotParetoFronts2D <- function(tuneParetoResult, objectives, drawLabels=TRUE,
                               drawBoundaries=TRUE, labelPos=4, fitLabels=TRUE,
                               cex.conf=0.5, lty.fronts=1, pch.fronts=8, col.fronts, ...)
{
  if(!inherits(tuneParetoResult, "TuneParetoResult"))
  {
    stop("\"tuneParetoResult\" must be a TuneParetoResult object!")
  }
  
  if (missing(objectives))
  {
    if (ncol(tuneParetoResult$testedObjectiveValues) == 2)
      objectives = c(1,2)
    else
      stop("Please supply exactly 2 objectives!")
  }
  
  if (missing(col.fronts))
    col.fronts <- NULL
    
  if (length(objectives) != 2)
  {
    stop("Please supply exactly 2 objectives!")
  }
  
  colorSet <- c("blue","green","red","darkgoldenrod","gold","brown","cyan",
      "purple","orange","seagreen","tomato","darkgray","chocolate",
      "maroon","darkgreen","gray12","blue4","cadetblue","darkgoldenrod4",
      "burlywood2")
  
  objectiveVals <- tuneParetoResult$testedObjectiveValues[,objectives,drop=FALSE]
  minimizeObjectives <- tuneParetoResult$minimizeObjectives[objectives]
  boundaries <- tuneParetoResult$objectiveBoundaries[objectives]
  
  internal.plotParetoFronts2D(objectiveVals, minimizeObjectives, boundaries, 
                              drawLabels=drawLabels, ,
                              drawBoundaries=drawBoundaries,
                              labelPos=labelPos,
                              fitLabels=fitLabels,                              
                              cex.conf=cex.conf,
                              lty.fronts=lty.fronts, 
                              pch.fronts=pch.fronts,
                              col.fronts=col.fronts,
                              ...)
}

# Plots the Pareto fronts of pairs of objectives in <tuneParetoResult>
# in a matrix-like plot. If <drawLabels> is true, the parameter configurations
# are printed in the plot.
# If <drawBoundaries is true, the boundaries of the objectives are plotted as lines.
# <cex.conf> is the text size of the configuration labels.
# <labelPos> is the position of the text relative to the points.
# If <fitLabels> is TRUE, overlaps of the parameter configuration labels are removed.
# <lty.fronts>, <pch.fronts> amd <col.fronts> specify the line types, the characters,
# and the colors for the Pareto fronts.
plotObjectivePairs <- function(tuneParetoResult, drawLabels=TRUE,
                               drawBoundaries=TRUE, labelPos=4, fitLabels=TRUE,
                               cex.conf=0.5, lty.fronts=1, pch.fronts=8, col.fronts,...)
{

  if (missing(col.fronts))
    col.fronts <- NULL

  xpos <- 0
  ypos <- 1
  
  pairs(data.frame(tuneParetoResult$testedObjectiveValues), 
        panel=function(x,y)
              {
                # calculate the current x and y position in the plot
                xpos <<- xpos + 1
                if (xpos == ypos)
                  xpos <<- xpos + 1
                  
                if (xpos > length(tuneParetoResult$minimizeObjectives))
                {
                  ypos <<- ypos + 1
                  xpos <<- 1
                }
                
                # build output matrix
                data <- data.frame(x,y)
                colnames(data) <- names(tuneParetoResult$minimizeObjectives)[c(xpos,ypos)]
                rownames(data) <- rownames(tuneParetoResult$testedObjectiveValues)
                
                # plot Pareto fronts
                internal.plotParetoFronts2D(data, 
                                            tuneParetoResult$minimizeObjectives[c(xpos,ypos)],
                                            tuneParetoResult$objectiveBoundaries[c(xpos,ypos)], 
                                            drawLabels=drawLabels, 
                                            drawBoundaries=drawBoundaries,
                                            labelPos=labelPos,         
                                            fitLabels=fitLabels,
                                            cex.conf=cex.conf,
                                            lty.fronts=lty.fronts, 
                                            pch.fronts=pch.fronts,
                                            col.fronts=col.fronts,
                                            plotNew=FALSE, ...)
              })
}

