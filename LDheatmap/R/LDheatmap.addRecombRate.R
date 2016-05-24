# ldheatmap - Plots measures of pairwise linkage disequilibria for SNPs
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

#__________Add Recomb Rate track from UCSC genome Browser to an LDheatmap_________##
LDheatmap.addRecombRate <- function(LDheatmap, chromosome, genome=NULL, recombRateLocation=0.02, view="dense") {
  minRange <- min(LDheatmap$genetic.distances)
  maxRange <- max(LDheatmap$genetic.distances)
#  minRange <- 129000000
#  maxRange <- 140000000

  flip <- !is.null(LDheatmap$flipVP)
  vp <- constructVP(LDheatmap$LDheatmapGrob, recombRateLocation, flip)
  vp$height <- unit(0.06, "npc")
  vp$name <- "recombRateVP"
  pushViewport(LDheatmap$LDheatmapGrob$vp) # the vps are necessary for the calculation of grobWidth
  if (flip) pushViewport(LDheatmap$flipVP) # in the recombRate function
  recombRate <- recombRate(minRange, maxRange, chromosome, genome, vp=vp, view)
  vp <- recombRate$vp
  vpstack <- vp; if(flip) vpstack <- vpStack(LDheatmap$flipVP,vp)
  grobT <- editGrob(recombRate, vp=vpstack)
  LDheatmap$LDheatmapGrob <-addGrob(LDheatmap$LDheatmapGrob, grobT)
  LDheatmap$LDheatmapGrob <- moveTitles(LDheatmap$LDheatmapGrob, vp)
  return(LDheatmap)
}

