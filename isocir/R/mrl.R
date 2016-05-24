"mrl"<-function(data){
# DATE WRITTEN: 4 Ago 2010          LAST REVISED:  05 Jan 2012
# AUTHOR: Sandra Barragan based on the SAS routines written by Miguel A. Fernandez.
# DESCRIPTION: data is a vector or a matrix of data where:
#                the rows are the items and 
#                the colums are the replications.
# REFERENCE: Mardia, K. and Jupp, P. (2000). Directional Statistics.
# SEE ALSO: CTi, CIREi.

if(all(c(!is.matrix(data),!is.vector(data),!is.data.frame(data)))){stop("data must be a matrix or a vector")}
if(!is.circular(data)){data <- suppressWarnings(as.circular(data))}

if(is.matrix(data)){ r <- suppressWarnings(apply(data,1,rho.circular))}
if(is.null(dim(data))){r <- suppressWarnings(rho.circular(data))}
return(r)
}