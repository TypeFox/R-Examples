# LDheatmap - Plots measures of pairwise linkage disequilibria for SNPs
# Copyright (C) 2004  J.Shin, S. Blay, N. Lewin-Koh, J.Graham, B.McNeney

# To cite LDheatmap in publications use:
# Shin J-H, Blay S, McNeney B and Graham J (2006). LDheatmap: An R
# Function for Graphical Display of Pairwise Linkage Disequilibria
# Between Single Nucleotide Polymorphisms. J Stat Soft, 16 Code Snippet 3

# This program is free software; you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation; either version 2 of the License, or
# (at your option) any later version.

# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.

# You should have received a copy of the GNU General Public License
# along with this program; if not, write to the Free Software
# Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA

###########################################################################


#______________________________________Add Scatter Plot________________________________##
LDheatmap.addScatterplot <- function(LDheatmap, P, height=0.2, ylab=NULL, ylim=NULL, type="points") {
  if (dim(LDheatmap$LDmatrix)[1] != length(P) ) {
     print ("Length of vector not equal number of SNPs in LDheatmap")
     return()
  }

  flip <- !is.null(LDheatmap$flipVP)
  vp <- constructVP(LDheatmap$LDheatmapGrob, 0.03, flip)
  vp$height <- unit(height, "npc")
  vp$name <- "associationVP"
  if(is.null(ylim)) ylim <- c(floor(min(P)), ceiling(max(P)))
  vp$yscale <- ylim
  vp$xscale <- c(min(LDheatmap$genetic.distances), max(LDheatmap$genetic.distances))
  xaxis <- linesGrob(x=vp$xscale, y=0, default.units="native", name="xaxis")
  yaxis <- linesGrob(x=min(LDheatmap$genetic.distances), y=vp$yscale, default.units="native",
		name="yaxis") 
  yaxisT <- yaxisGrob(name="yaxis_ticks", gp=gpar(fontsize=7))
  ylab <- textGrob(ylab, rot=90, gp=gpar(fontsize=9), name="yaxis_title",
	x=unit(min(LDheatmap$genetic.distances), "native")- unit(10, "millimeters"))
  vpstack <- vp; if(flip) vpstack <- vpStack(LDheatmap$flipVP,vp)
  association <- gTree(children=gList(xaxis, yaxis, yaxisT, ylab), name="association", vp=vpstack)
  if (type=="points" || type == "both") {
     graph_points <- pointsGrob(LDheatmap$genetic.distances, P, size=unit(2, "millimeters"),
 		name="points")
     association <- addGrob(association, graph_points)
  }
  if (type=="lines" || type == "both") {
     graph_lines <- linesGrob(LDheatmap$genetic.distances, P, default.units="native", name="lines")
     association <- addGrob(association, graph_lines)
  }
  LDheatmap$LDheatmapGrob <- addGrob(LDheatmap$LDheatmapGrob, association)
  LDheatmap$LDheatmapGrob <- moveTitles(LDheatmap$LDheatmapGrob, vp)
  return(LDheatmap)
}


