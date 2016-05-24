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



#_______________________Add marks to heatmap__________________________________##
# Adds a symbol to the i,jth cell of the heatmap. 
# The default is to add a symbol to the diagonal (j=i). 
# i and j can be vectors.
LDheatmap.marks <- function(LDheatmap, i, j=NULL, pch=20, gp=gpar(...), ...){
    nSNP <- dim(LDheatmap$LDmatrix)[1]
    if(is.null(j)) j<-i
    ind <- i>j
    if(any(ind)){
      ind <- as.numeric(ind)
      ind[ind>0] <- i[ind>0]
      i[ind>0] <- j[ind>0]
      j[ind>0] <- ind[ind>0]
    }
    pts<-list(x=(i-0.5)*1/nSNP,y=(j-0.5)*1/nSNP)
    if(!is.null(LDheatmap$flipVP)) pts<-list(x=(j-0.5)*1/nSNP,y=(i-0.5)*1/nSNP)
    heatmap.vp <- LDheatmap$heatmapVP$name
    #If heatmap.vp is on the grid display list, i.e., it is included in the 
    #returned value of current.vpTree(), a[1] <- 1 else a[1] <- NA
    a <- grep(paste("[", heatmap.vp, "]", sep=""), as.character(current.vpTree()), fixed=TRUE)
    if(!is.na(a[1]))   seekViewport(heatmap.vp)
    else               pushViewport(LDheatmap$heatmapVP)
    if (!is.null(LDheatmap$flipVP)) pushViewport(LDheatmap$flipVP)
    Symbols <- pointsGrob(pts$x, pts$y, pch=pch, gp=gp, name="symbols")
    symbols <- gTree(children=gList(Symbols), name="Symbols", cl="symbols")
    grid.draw(symbols)
    if(!is.na(a[1]))  upViewport(0)  #back to the root viewport
    else              popViewport() 
    invisible(pts)
}


