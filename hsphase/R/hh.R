# Copyright (C) 2014 Mohammad H. Ferdosi
#
# HSPhase is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.

# HSPhase program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.

# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <http:#www.gnu.org/licenses/>.

hh <- function (oh, inferredPedigree, realPedigree, pedOnly = TRUE) 
{
  if (is.null(oh)) 
    stop("Invalid input!")
  if (!is.matrix(oh)) 
    stop("Opposing Homozygote matrix should be a MATRIX")
  if (nrow(oh) != ncol(oh)) 
    stop("Opposing Homozygote matrix must be a square matrix ")
  if (!missing(inferredPedigree) && !missing(realPedigree)) {
#     if (nrow(inferredPedigree) != nrow(realPedigree)) 
#       stop("Pedigrees must have the same number of rows")
    colour <- c("#E41A1C", "#0000FF", "#00FF00", "#FFA500", 
                "#E7298A", "#000080", "#00FFFF", "#008000", "#999999", 
                "#800080", "#008080", "#000000", "#FF4500", "#FFFF00", 
                "#4DAF4A", "#800000", "#A65628", "#F781BF", "#1B9E77", 
                "#808000", "#377EB8")
    maxrows <- max(nrow(inferredPedigree),nrow(realPedigree))
    colour <- rep(colour,ceiling(maxrows/length(colour)))
    
    if(pedOnly == TRUE)
       oh = oh[rownames(oh)%in%realPedigree[,1],colnames(oh)%in%realPedigree[,1]]
    
    index1 = match(colnames(oh), realPedigree[, 1])
    index2 = match(colnames(oh), inferredPedigree[, 1])
  
    sires <- unique(c(names(table(as.character(realPedigree[,2]))),names(table(as.character(inferredPedigree[,2])))))
    
    realPedigree <- as.character(realPedigree[index1,2])
    inferredPedigree <- as.character(inferredPedigree[index2,2])
    
    for(i in 1:length(sires))
    {
      realPedigree[which(realPedigree==sires[i])] <- colour[i]
      inferredPedigree[which(inferredPedigree==sires[i])] <- colour[i]
    }
    
    heatmap(oh, symm = T, col = gray.colors(16, start = 0, 
                                            end = 1), RowSideColors = inferredPedigree , ColSideColors = realPedigree )
  }
  else heatmap(oh, symm = T, col = gray.colors(16, start = 0, 
                                               end = 1))
}