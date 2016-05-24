###########################################################################
# Copyright 2009 Nobody                                                   #
#                                                                         #
# This file is part of hergm.                                             #
#                                                                         # 
#    hergm is free software: you can redistribute it and/or modify        #
#    it under the terms of the GNU General Public License as published by #
#    the Free Software Foundation, either version 3 of the License, or    #
#    (at your option) any later version.                                  #
#                                                                         # 
#    hergm is distributed in the hope that it will be useful,             #
#    but WITHOUT ANY WARRANTY; without even the implied warranty of       #
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the        #
#    GNU General Public License for more details.                         #
#                                                                         #
#    You should have received a copy of the GNU General Public License    #
#    along with hergm.  If not, see <http://www.gnu.org/licenses/>.       #
#                                                                         # 
###########################################################################

hergm.plot <- function(sample = NULL,
                       ...)
# input: network, postprocess output
# output: plot of network block membership probabilities
{
  # Extract
  network <- sample$network 

  # Plot
  if (is.directed(network)) gmode <- "digraph"
  else gmode <- "graph"
  p <- gplot(network, gmode=gmode, mode="fruchtermanreingold", vertex.cex=1, vertex.col=0, vertex.border=0, displaylabels=TRUE, label.cex=0.8)
  for(i in 1:nrow(p))
    {
    ergmm.drawpie(center=p[i,], radius=0.25, probs=sample$p_i_k[i,])
    }

}

