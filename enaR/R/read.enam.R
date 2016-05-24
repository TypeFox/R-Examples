#' R function to read in a matrix formatted as Mdloti (Ursula Sharler)
#' Borrett | Sept. 12, 2012, MKL July 2013
#' ------------------------

read.enam<- function(file="file path and name"){
                                        #I have assumed the file is formatted as an excel speadsheet.
                                        #The data must be on the first sheet in the workbook.
  x <- as.matrix(read.xls(file,sheet=1,header=FALSE))
  mname <- as.character(x[1,1]); # Get Model ID
  n <- as.numeric(as.character(x[2,2])) # number of nodes
  liv <- as.numeric(as.character(x[3,2])) # number of nodes
  a <- n+6+1 # ending row of flows matrix -- assumes Flows start on row 6 and Imports and Biomasses are at the end
  b <- n+2+2 # ending column of flows matrix -- assumes exports and respirations are at the end
  m <- x[6:a,3:b] # Matrix of Flows
  m <- apply(m,2,as.numeric)
  rownames(m) <- colnames(m) <- as.character(x[6:a,2]) # node names
  m[is.na(m)] <- 0 # replace NAs with zeros
  Flow <- m[1:n,1:n] # flow matrix
  imports <- m[(n+1),1:n]
  biomass <- as.numeric(unlist(m[(n+2),1:n]))
  exports <- as.numeric(unlist(m[1:n,(n+1)]))
  respiration <- as.numeric(unlist(m[1:n,(n+2)]))
  y <- network(Flow,directed=TRUE)
                                        # packing up the attributes into the network object (y)
  set.vertex.attribute(y,'input',imports)
  set.vertex.attribute(y,'export',exports)
  set.vertex.attribute(y,'respiration',respiration)
  set.vertex.attribute(y,'storage',biomass)
  set.vertex.attribute(y,'output',exports+respiration)
  y%v%'vertex.names' <- rownames(Flow)
  y%v%'living'=c(rep(TRUE,liv),rep(FALSE,n-liv))
  set.network.attribute(y, 'flow', Flow[Flow>0])
  return(y)
}
