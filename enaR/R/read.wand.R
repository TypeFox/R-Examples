#' read.wand --- WAND formatted file into R
#' INPUT = file path
#' OUTPUT = network object
#' S. Borrett | May 2012
#' ------------------------------------

read.wand <- function(file='file name with path'){
                                        # file is the full excel file name
                                        # asssumes that first sheet is "Main" and second sheet is "Flows".
  x <- as.matrix(read.xls(file,sheet="Main"))
  d1 <- x[1:8,1] #model info
  n <- as.numeric(as.character(d1[3])) #Number of compartments
  dat.main <- x[8:(n+9),2:6] #isolate the stocks,imports,exports,respirations
  dat.main[is.na(dat.main)] <- 0 #zero NA
  vn <- dat.main[1:n,1] #vertex names
  dat.main <- apply(dat.main[,2:5],2,as.numeric)
                                        # get flows
  Flow <- read.xls(file,sheet="Flows")
  flow.mat <- as.matrix(Flow[1:(n),2:(n+1)])
  flow.mat[is.na(flow.mat)] <- 0
  flow.mat <- apply(flow.mat,2,as.numeric)
  rownames(flow.mat) <- colnames(flow.mat) <- vn
                                        #pack for export
  x <- list('flow'=flow.mat,
    'input'=dat.main[,2],
    'exports'=dat.main[,3],
    'respiration'=dat.main[,4],
    'storage'=dat.main[,1])
  ## --- Create Network Object From Data ---
  y <- network(x[[1]],directed=TRUE)
                                        # packing up the attributes into the network object (y)
  set.vertex.attribute(y,'input',x[[2]])
  set.vertex.attribute(y,'export',x[[3]])
  set.vertex.attribute(y,'respiration',x[[4]])
  set.vertex.attribute(y,'storage',x[[5]])
  set.vertex.attribute(y,'output',x[[3]]+x[[4]])
  y%v%'vertex.names' <- vn
  set.edge.attribute(y,'flow', flow.mat[flow.mat>0])
  return(y)
}
