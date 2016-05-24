### network.make.dot  (2008-08-18)
###
###   Generate Dot File For Graphviz Network Plot
###
### Copyright 2006-08 Rainer Opgen-Rhein and Korbinian Strimmer
###
###
### This file is part of the `GeneNet' library for R and related languages.
### It is made available under the terms of the GNU General Public
### License, version 3, or at your option, any later version,
### incorporated herein by reference.
### 
### This program is distributed in the hope that it will be
### useful, but WITHOUT ANY WARRANTY; without even the implied
### warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
### PURPOSE.  See the GNU General Public License for more
### details.
### 
### You should have received a copy of the GNU General Public
### License along with this program; if not, write to the Free
### Software Foundation, Inc., 59 Temple Place - Suite 330, Boston,
### MA 02111-1307, USA


network.make.dot = function(filename, edge.list, node.labels, 
          main=NULL, show.edge.labels=FALSE)
{
   f = file(filename, "w")  # open an output file connection

   w = edge.list[,1]  # vector of weights
   el = as.character(round(edge.list[,1], digits=4) )
   n1 = edge.list[,2]
   n2 = edge.list[,3]  

   # thresholds for line width and coloring
   cutoff = quantile(abs(w), c(0.2, 0.8)) 

   cat("Digraph G {\n", file=f)
   cat("K=\"0.3\";\n", file=f)
   cat("ratio=\"0.75\";\n", file=f)

   if(!is.null(main))
   {
     cat("label = \"\\n\\n", main,"\";\n", file=f)
     cat("fontsize=20;\n", file=f)
   }

   # check whether we might have directed edges
   if (ncol(edge.list) == 11) 
     undirected=FALSE
   else 
     undirected=TRUE


   for(i in 1:length(w))   
   {
      cat("\"", node.labels[n1[i]], "\" -> \"", 
        node.labels[n2[i]], "\" [", sep="", file = f)
     
      if (undirected)
        cat("dir=\"none\", ", sep="", file=f)
      else
      {
        if(edge.list[i,11]=="1to2")
        {
          cat("dir=\"forward\",", file=f)
        }
        else if (edge.list[i,11]=="2to1")
        {
          cat("dir=\"back\",", file=f)
        }
        else
        {
          cat("dir=\"none\",", file=f)
        }
      } 

      if(show.edge.labels)
      {
         cat("label=\"", el[i], "\", ", sep="", file=f)
      }
     
      # line thickness and color depends on relative strengh	
      if (abs(w[i]) < cutoff[1]) # lower 20% quantile
      {
         if(w[i] < 0) 
           cat("color=grey, style=\"dashed,setlinewidth(1)\"", file=f)
         else 
           cat("color=grey, style=\"solid,setlinewidth(1)\"", file=f)
      }
      else if (abs(w[i]) < cutoff[2]) # from 20%-80%
      {
         if(w[i] < 0) 
           cat("color=black, style=\"dashed,setlinewidth(1)\"", file=f)
         else 
           cat("color=black, style=\"solid,setlinewidth(1)\"", file=f)
       }
      else # top 80%-100%
      {
          if(w[i] < 0) 
           cat("color=black, style=\"dashed,setlinewidth(2)\"", file=f)
         else 
           cat("color=black, style=\"solid,setlinewidth(2)\"", file=f)
      }

      cat("];\n", file=f)
   }
   cat("}\n", file=f)
   close(f)
}
