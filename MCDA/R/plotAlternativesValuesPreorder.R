#############################################################################
#
# Copyright Patrick Meyer and Sébastien Bigaret, 2013
#
# Contributors:
#   Patrick Meyer <patrick.meyer@telecom-bretagne.eu>
#   Sébastien Bigaret <sebastien.bigaret@telecom-bretagne.eu>
#		
# This software, MCDA, is a package for the R statistical software which 
# allows to use MCDA algorithms and methods. 
# 
# This software is governed by the CeCILL license (v2) under French law
# and abiding by the rules of distribution of free software. You can
# use, modify and/ or redistribute the software under the terms of the
# CeCILL license as circulated by CEA, CNRS and INRIA at the following
# URL "http://www.cecill.info".
# 
# As a counterpart to the access to the source code and rights to copy,
# modify and redistribute granted by the license, users are provided only
# with a limited warranty and the software's author, the holder of the
# economic rights, and the successive licensors have only limited
# liability.
#		
# In this respect, the user's attention is drawn to the risks associated
# with loading, using, modifying and/or developing or reproducing the
# software by the user in light of its specific status of free software,
# that may mean that it is complicated to manipulate, and that also
# therefore means that it is reserved for developers and experienced
# professionals having in-depth computer knowledge. Users are therefore
# encouraged to load and test the software's suitability as regards their
# requirements in conditions enabling the security of their systems and/or
# data to be ensured and, more generally, to use and operate it in the
# same conditions as regards security.
#		
# The fact that you are presently reading this means that you have had
# knowledge of the CeCILL license and that you accept its terms.
#
##############################################################################

plotAlternativesValuesPreorder <- function(alternativesValues, decreasing = TRUE, alternativesIDs = NULL){
  
  # if (!require(Rgraphviz)) stop("Rgraphviz package could not be loaded")
  
  if (!requireNamespace("Rgraphviz", quietly = TRUE)) stop("Rgraphviz package could not be loaded")
  
  ## check the input data

  if (!(is.vector(alternativesValues)))
    stop("alternativesValues should be a vector")
  
  if (!(is.null(alternativesIDs) || is.vector(alternativesIDs)))
    stop("alternativesIDs should be a vector")
  
  ## filter the data according to the given alternatives and criteria
  
  if (!is.null(alternativesIDs)){
    alternativesValues <- alternativesValues[alternativesIDs]
  } 
  
  # data is filtered, check for some data consistency
  
  # if there are less than 2 criteria or 2 alternatives, there is no MCDA problem
  
  if (length(alternativesValues)<2) 
    stop("less than 2 criteria or 2 alternatives")
  
  # -------------------------------------------------------
  
  numAlt <- length(alternativesValues)
  
  
  
  graph<-NULL
  
  # we first order the ranks from best to worst
  
  orderedRanks<-alternativesValues[order(alternativesValues,decreasing=decreasing)]
  
  # we now construct the nodes of the graph (alternatives with equal rank are put together into one node)
  
  prevRank<-orderedRanks[1]
  curLab<-names(orderedRanks)[1]
  labels<-curLab
  
  for (i in 2:length(orderedRanks))
  {
    
    if (orderedRanks[i] != orderedRanks[i-1]){
      
      labels<-c(labels, names(orderedRanks)[i])
    }
    else
    {
      labels[length(labels)] <- paste(labels[length(labels)],names(orderedRanks)[i],sep=",")
    }
  }
  
  # we now construct the edges of the graph
  
  edg<-vector("list",length=length(labels))
  names(edg)<-labels
  if (length(labels)==1){
    edg[[1]]<-list(edges=character(0))
  }
  else{
    for (i in 1:(length(labels)-1))
    {
      edg[[i]]<-list(edges=labels[i+1])
    }
    edg[[length(labels)]]<-list(edges=character(0))
  }
  
  # finally we construct the graph
  
  graph <- new("graphNEL", nodes=labels, edgeL=edg, edgemode="directed")
  
  Rgraphviz::plot(graph, attrs = list(node = list(shape = "box", fixedsize = FALSE)))  
}
