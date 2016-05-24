#    Function: Determines different response patterns in a sample.
#    Copyright (C) 2011  David Preinerstorfer
#    david.preinerstorfer@univie.ac.at
#
#    This program is free software; you can redistribute it and/or modify
#    it under the terms of the GNU General Public License as published by
#    the Free Software Foundation; either version 2 of the License, or
#    (at your option) any later version.
#
#    This program is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#    GNU General Public License for more details. A copy may be obtained at
#    http://www.r-project.org/Licenses/

resp.patterns <-function(data.matrix){

data.matrix <- na.omit(as.matrix(data.matrix))

#recode answer vectors to strings
fu <- function(x){paste(x, collapse = "")}
pat <- apply(data.matrix, 1, fu)

#count patterns
tab <- table(pat)

#return pattern-strings into answer vectors
patterns <- matrix(as.numeric(unlist(strsplit(names(tab), ""))), ncol = dim(data.matrix)[2], byrow = TRUE)

#compute scores
scores <- apply(patterns, 1, sum)

#return different patterns, their counts and scores
return(cbind(patterns, "count" = as.numeric(tab), scores)[order(scores),])
}

