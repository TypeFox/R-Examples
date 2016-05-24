#******************************************************************************* 
#
# Bayesian Regression and Adaptive Sampling with Gaussian Process Trees
# Copyright (C) 2005, University of California
# 
# This library is free software; you can redistribute it and/or
# modify it under the terms of the GNU Lesser General Public
# License as published by the Free Software Foundation; either
# version 2.1 of the License, or (at your option) any later version.
# 
# This library is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
# Lesser General Public License for more details.
# 
# You should have received a copy of the GNU Lesser General Public
# License along with this library; if not, write to the Free Software
# Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301  USA
#
# Questions? Contact Robert B. Gramacy (rbgramacy@ams.ucsc.edu)
#
#*******************************************************************************


## tgp.trees:
##
## plot the MAP tree found at each tree in the Markov chain
## for the tgp-class object (or, possibly constrin the plotting
## to certain heights -- requires the maptree library for plotting

"tgp.trees" <-
function(out, heights=NULL, main=NULL, ...)
{
  ## get the full set of heights if none specified, and length
  if(is.null(heights)) heights <- out$posts$height
  else if(heights == "map") { ## only plot the MAP
    heights <- out$post$height[which.max(out$posts$lpost)]
  }
  howmany <- length(heights)

  ## calculate how many sub-windows to make with par
  if(howmany > 1) {
    h <- howmany
    if(sum(out$posts$height == 1) >= 1) { h <- h - 1; }
    rows <- floor(sqrt(h)); cols <- floor(h / rows)
    while(rows * cols < h) cols <- cols + 1
    par(mfrow=c(rows, cols), bty="n")
  } else par(mfrow=c(1,1), bty="n")

  ## create a vector of names for the main text section of each plot
  names <- names(out$X)
  if(is.null(names)) {
    for(i in 1:out$d) { names <- c(names, paste("x", i, sep="")) }
  }

  ## plot each tree
  for(j in 1:howmany) { 
    if(is.null(out$trees[[heights[j]]])) next;
    p <- (1:length(out$posts$height))[out$posts$height == heights[j]]
    tgp.plot.tree(out$trees[[heights[j]]], names, out$posts[p,],
                  main=main, ...); 
  }
}


## tgp.plot.tree:
##
## actually use maptree to plot each tree specified in the
## tree frame with specified name and posterior probability

"tgp.plot.tree" <-
function(frame, names, posts, main=NULL, ...)
{
  ## don't plot (null) trees of height one
  if(dim(frame)[1] == 1) {
    cat(paste("NOTICE: skipped plotting tree of height 1, with lpost =", 
              posts$lpost, "\n"))
    return()
  }

  ## concatenate the log-posterior probability to the main text
  main <- paste(main, " height=", posts$height, ", log(p)=",
                posts$lpost, sep="")
  
  ## create a frame vector that maptree understands
  frame[,2] <- as.character(frame[,2])
  n.i <- frame[,2] != "<leaf>"
  frame[n.i,2] <- names[as.numeric(frame[n.i,2])+1]
  frame[,2] <- factor(frame[,2])
  splits <- as.matrix(data.frame(cutleft=as.character(frame[,6]),
                                 cutright=as.character(frame[,7])))
  new.frame <- data.frame(frame[,2:5], splits=I(splits), row.names=frame[,1])
  tree <- list(frame=new.frame)

  ## draw the tree and add a title
  draw.tree(tree, ...)
  title(main)
}

