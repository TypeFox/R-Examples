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


#_______________________Highlight a region in the heatmap____________________________##
LDheatmap.highlight <- function(LDheatmap, i, j, fill="NA", col="black", lwd=1, lty=1){

  # Highlights the perimeter of selected cells in the heatmap as a block
  backbone <- function(i,j,nSNP){
     x <- c(i-1,i-1,j-1)/nSNP
     y <- c(i,j,j)/nSNP
     cbind(x,y)
  }
  backboneFlip <- function(i,j,nSNP){
     x <- c(i,j,j)/nSNP
     y <- c(i-1,i-1,j-1)/nSNP
     cbind(x,y)
  }

  zigzag <- function(i,j,nSNP){
    c1 <- j-i
    nvert <- (2*c1)-1
    x <-c(j-1,rep((j-2):(j-c1),each=2))
    y <- c(rep((j-1):(j-(c1-1)),each=2),j-c1)
    cbind(x,y)/nSNP 
  }
  zigzagFlip <- function(i,j,nSNP){
    c1 <- j-i
    nvert <- (2*c1)-1
    y <-c(j-1,rep((j-2):(j-c1),each=2))
    x <- c(rep((j-1):(j-(c1-1)),each=2),j-c1)
    cbind(x,y)/nSNP 
  }
                       
  nSNP <- dim(LDheatmap$LDmatrix)[1]
  if(length(i)>1 | length (j) > 1) stop("i and j must be scalar indices")
  if((i<1 | i>nSNP) |(j<1 | j>nSNP) )
    stop(paste("index out of bounds, i and j must be in (1,",nSNP,")",sep=""))
  if(i==j) stop("i cannot be equal to j")
  if(i>j){
     h<-i
     i <- j
     j <- h
  }
  pgon <- data.frame(rbind(backbone(i,j,nSNP), zigzag(i,j,nSNP)))
  if(!is.null(LDheatmap$flipVP)) pgon <- data.frame(rbind(backboneFlip(i,j,nSNP), zigzagFlip(i,j,nSNP)))
  ## Square or almost square interior Blocks
  names(pgon) <- c("x","y")
  heatmap.vp <- LDheatmap$heatmapVP$name
  #If heatmap.vp is on the grid display list, i.e., it is included in the 
  #returned value of current.vpTree(), a[1]=1 else a[1]=NA:
  a <- grep(paste("[", heatmap.vp, "]", sep=""), as.character(current.vpTree()), fixed=TRUE)
  if(!is.na(a[1]))   seekViewport(heatmap.vp)
  else               pushViewport(LDheatmap$heatmapVP)
  if (!is.null(LDheatmap$flipVP)) pushViewport(LDheatmap$flipVP)
  highlight <- polygonGrob(x=pgon$x, y=pgon$y, 
     gp=gpar(col=col, fill=fill, lwd=lwd, lty=lty), name="highlight")
  grid.draw(highlight)
  if(!is.na(a[1]))  upViewport(0)  #back to the root viewport
  else              popViewport() 
  invisible(pgon)
}


