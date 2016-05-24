
# ===============================
# private function
# NEAR NEIGHBOUR FINDER based on ANN C++ library and RANN functions
# This function is a modified version of the nn2.R function of the RANN package,
# copyrigth under Samuel Kemp 2005-9 and Gregory Jefferis 2009-2013.
# Date of last edit: 17-09-2013
# ===============================
nn.search <- function(data, query, k=min(10,nrow(data)),treetype=c("kd","bd"),
                      searchtype=c("standard","priority","radius"),radius=0.0,eps=0.0)
{
  dimension  <- ncol(data)
  ND  	    <- nrow(data)
  NQ		    <- nrow(query)
  
  # Check that both datasets have same dimensionality
  if(ncol(data) != ncol(query) )
    stop("Query and data tables must have same dimensions")	
  
  if(k>ND)
    stop("Cannot find more nearest neighbours than there are points")
  
  searchtypeInt=pmatch(searchtype[1],c("standard","priority","radius"))
  if(is.na(searchtypeInt)) stop(paste("Unknown search type",searchtype))
  treetype=match.arg(treetype,c("kd","bd"))
  
  # Coerce to matrix form
  if(is.data.frame(data))
    data <- unlist(data,use.names=FALSE)
  
  # Coerce to matrix form
  if(!is.matrix(query))
    query <- unlist(query,use.names=FALSE)
  
  # void get_NN_2Set(double *data, double *query, int *D, int *ND, int *NQ, int *K, double *EPS,
  # int *nn_index, double *distances)
  
  results <- .C("get_NN_2Set",
                as.double(data),
                as.double(query),
                as.integer(dimension),
                as.integer(ND),
                as.integer(NQ),
                as.integer(k),
                as.double(eps),
                as.integer(searchtypeInt), 
                as.integer(treetype=="bd"), 
                as.double(radius*radius),
                nn.idx   = integer(k*NQ),
                nn       = double(k*NQ), PACKAGE="nonlinearTseries")
  
  # now put the returned vectors into (nq x k) arrays
  nn.indexes=matrix(results$nn.idx,ncol=k,byrow=TRUE)
  nn.dist=matrix(results$nn,ncol=k,byrow=TRUE)
  
  return(list(nn.idx=nn.indexes, nn.dists=nn.dist^2))
}
