#This file belongs to
#c3net: C3NET, <https://r-forge.r-project.org/projects/c3net/>
#This R package allows inferring regulatory networks from expression data using C3NET.
#The inferred network consists of only direct physical interactions.
## Copyright (C) January 2011 Gokmen Altay <altayscience@gmail.com>
## This program is a free software for only academic useage but not for commercial useage; you can redistribute it and/or
## modify it under the terms of the GNU GENERAL PUBLIC LICENSE
## either version 3 of the License, or any later version.
##
## This program is distributed WITHOUT ANY WARRANTY; 
## You can get a copy of the GNU GENERAL PUBLIC LICENSE
## from
## http://www.gnu.org/licenses/gpl.html
## See the licence information for the dependent package from
## igraph package itself.


checknet <- function(finalrelationmatrix, realrelationmatrix)
{
finalrelationmatrix[finalrelationmatrix != 0] <- 1 
TP <- sum( finalrelationmatrix * realrelationmatrix )
CompRealrelationmatrix <- abs(realrelationmatrix-1)
FP <- sum( finalrelationmatrix * CompRealrelationmatrix )
CompFinalrelationmatrix <- abs(finalrelationmatrix-1)
FN <- sum( realrelationmatrix * CompFinalrelationmatrix )
TP <-TP/2
FP <-FP/2
FN <-FN/2
precision <- TP/(TP+FP); recall <- TP/(TP+FN)
Fscore <- 2*precision*recall/(precision+recall)
output <- c(precision, Fscore, recall, TP, FP, FN)
namesv <- c("precision", "F-score", "recall", "TP", "FP", "FN")
names(output)<-namesv

output
}



