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

#__________Construct a viewport for additional graphs and gene tracks__________________##
constructVP <- function(LDheatmapGrob, location=0, flip) {

  x0 <- convertX(getGrob(LDheatmapGrob, "diagonal")[[1]][1], "npc", valueOnly=TRUE)
  x1 <- convertX(getGrob(LDheatmapGrob, "diagonal")[[1]][2], "npc", valueOnly=TRUE)
  y0 <- convertX(getGrob(LDheatmapGrob, "diagonal")[[2]][1], "npc", valueOnly=TRUE)
  y1 <- convertX(getGrob(LDheatmapGrob, "diagonal")[[2]][2], "npc", valueOnly=TRUE)
  map_len = sqrt((x1-x0)^2 + (y1-y0)^2)		# genetic map length in npc units

  g_height <- g_x0 <- g_y0 <- 0
  if(!is.null(getGrob(LDheatmapGrob, "transcripts"))) {	# if gene track has been plotted
    transcriptsVP <- getGrob(LDheatmapGrob, "transcripts")$vp
    if (flip) transcriptsVP <- transcriptsVP[[2]]
    g_x0 <- convertX(transcriptsVP$x, "npc", valueOnly=TRUE)
    g_y0 <- convertX(transcriptsVP$y, "npc", valueOnly=TRUE)
    g_height <- convertX(transcriptsVP$height, "npc", valueOnly=TRUE)
  }

  r_height <- r_x0 <- r_y0 <- 0  
  if(!is.null(getGrob(LDheatmapGrob, "recombRate"))) {	# if recombRate track has been plotted
    recombRateVP <- getGrob(LDheatmapGrob, "recombRate")$vp
    if (flip) recombRateVP <- recombRateVP[[2]]
    r_x0 <- convertX(recombRateVP$x, "npc", valueOnly=TRUE)
    r_y0 <- convertX(recombRateVP$y, "npc", valueOnly=TRUE)
    r_height <- convertX(recombRateVP$height, "npc", valueOnly=TRUE)
  }

  m_height <- m_x0 <- m_y0 <- 0  
  if(!is.null(getGrob(LDheatmapGrob, "association"))) {	# if association scatterplot has been plotted
    assocVP <- getGrob(LDheatmapGrob, "association")$vp
    if (flip) assocVP <- assocVP[[2]]
    m_x0 <- convertX(assocVP$x, "npc", valueOnly=TRUE)
    m_y0 <- convertX(assocVP$y, "npc", valueOnly=TRUE)
    m_height <- convertX(assocVP$height, "npc", valueOnly=TRUE)
  }

  # Set the viewport
  if (!flip) {				# flip = FALSE
     angle <- 45
     genome_vp_just <- c("left", "top")
  }
  else {				# flip = TRUE
     angle <- 45
     genome_vp_just <- c("left", "bottom")
  }
  vp <- viewport(angle=angle, just=genome_vp_just, width=map_len,
	x=min(x0,  g_x0 - g_height*0.8, r_x0 - r_height*0.8, m_x0 - m_height*0.8) - location, 
	y=max(y0,  g_y0 + g_height*0.8, r_y0 + r_height*0.8, m_y0 + m_height*0.8) + location)
  return (vp)
}

#_______________________Move titles if covered up by additional plots____________________##
moveTitles <- function(LDheatmapGrob, vp) {

  genemap_title_y <- convertX(getGrob(LDheatmapGrob,"geneMap::title")$y,
                        "npc", valueOnly=TRUE)
  genemap_title_x <- convertX(getGrob(LDheatmapGrob,"geneMap::title")$x,
                        "npc", valueOnly=TRUE)
  flipVP <- getGrob(LDheatmapGrob,"geneMap::diagonal")$vp

  if (is.null(flipVP)) {					# flip = FALSE
     if (genemap_title_y == 0.3 & genemap_title_x == 0.5)	# user used default setting
	LDheatmapGrob <- editGrob(LDheatmapGrob,"geneMap::title", y=unit(0.1, "npc"), 
			x=unit(1.1, "npc"), just="right")
     grid.newpage()
     grid.draw(LDheatmapGrob)
     return(LDheatmapGrob)
  }

  # Get top of viewport coordinates in inches on the device
  pushViewport(LDheatmapGrob$vp)
  pushViewport(flipVP)
  vp_trans <- current.transform()
  temp <- c(
convertX(vp$x, "inches", valueOnly=TRUE) - convertX(vp$height, "inches",valueOnly=TRUE)/sqrt(2), 
convertX(vp$y, "inches", valueOnly=TRUE) + convertX(vp$height, "inches", valueOnly=TRUE)/sqrt(2), 1)
  upViewport()
  tr <- temp %*% vp_trans  # (x, y, 1) on device

  # Get genemap title coordinates in inches on the device
  vp_trans1 <- current.transform()
  genemap_title_y_inch <- convertY(getGrob(LDheatmapGrob,"geneMap::title")$y, "inches", valueOnly=TRUE)
  temp1 <- c(0, genemap_title_y_inch, 1)
  tr1 <- temp1 %*% vp_trans1


  # Move gene map title # if necessary
  new_genemap_title_y_inch <- t(solve(t(vp_trans1), t(tr)))[1,2]
  new_genemap_title_y_npc <- convertY(unit(new_genemap_title_y_inch, "inches"), "npc", valueOnly=TRUE) + 0.05

  LDheatmapGrob <- editGrob(LDheatmapGrob, "geneMap::title", 
	y = unit(new_genemap_title_y_npc, "npc"), just=c("left","bottom"))


  # Move heat map title if necessary
  heatmap_title_y <- convertY(getGrob(LDheatmapGrob,"heatMap::title")$y, "npc", valueOnly=TRUE)
  genemap_title_height <- convertHeight(grobHeight(getGrob(LDheatmapGrob,"geneMap::title")),
 		"npc", valueOnly=TRUE)

  if (heatmap_title_y < new_genemap_title_y_npc + genemap_title_height*3) {
     new_heatmap_title_y <- new_genemap_title_y_npc + genemap_title_height*3
	 LDheatmapGrob <- editGrob(LDheatmapGrob, "heatMap::title", 
		y = unit(new_heatmap_title_y, "npc"))
  }

  drawLDheatmapGrob(LDheatmapGrob)
  return(LDheatmapGrob)
}

drawLDheatmapGrob <- function(LDheatmapGrob) {
  heatmap_title_y <- convertY(getGrob(LDheatmapGrob,"heatMap::title")$y, "npc", valueOnly=TRUE)
  vp = viewport(height=1/heatmap_title_y, width=1, y=0.05, just="bottom",
        gp=gpar(cex=1/heatmap_title_y), name="container")
  grid.newpage()
  pushViewport(vp)
  grid.draw(LDheatmapGrob)
  popViewport()
}


