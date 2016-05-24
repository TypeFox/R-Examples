indexDatacube = function(
##title<< Create logical index matrices for multidimensional datacubes
  datacube = c() ##<< array: datacube from which to extract the sub-parts
                 ##   datacube and dims.datacube should be supplied.
  , logical.ind  ##<< logical array: TRUE/FALSE index matrix for a subset of the 
                 ##   dimensions of the datacube. The size of logical.ind`s dimensions 
                 ##   has to match the sizes of the corresponding dimensions in 
                 ##   datacube.
  , dims='auto'  ##<< integer vector or 'auto' : indices of the dimensions in 
                 ##   datacube corresponding to the dimensions of logical.ind. 
                 ##   If set to 'auto' this matching is tried to be accomplished 
                 ##   by comparing the sizes of the dimensions of the two objects.
  , dims.datacube = dim(datacube) ##<< integer vector: dimensions of the datacube. Only one of
                 ##   dims.datacube or datacube should be supplied!
             
)
  ##description<< This function facilitates supplying logical index array for only some but not all 
  ##              of the dimensions of a data array. This mimics Matlabs indexing scheme.
  ##              The indexing mechanisms of R only allow supplying logical indices for all
  ##              dimensions.
{
  if (sum(logical.ind) == 0) 
    stop('No TRUE value in index matrix!')
  if (dims[1] != 'auto' && length(dims) != length(dim(logical.ind)))
    stop('Argument dims needs to have one entry for each dimension in logical.ind!')
  if (dims[1] != 'auto' && max(dims) > length(dims.datacube))
    stop('Index in dims exceed dimensions of target datacube!')
  if (length(datacube) != 0 && !identical(dim(datacube), dims.datacube))
    stop('dims.datacube does not match dim(datacube)! Supplying values for both does not make sense!')
  
  
  if (dims[1] == 'auto') {
    if (is.null(dim(logical.ind)[1])) {
      size.ind     <- length(logical.ind)
      logical.ind  <- matrix(logical.ind, ncol=1)
    } else {
      size.ind     <- dim(logical.ind)
    }
    dims             <- match(size.ind, dims.datacube)
    if (sum(duplicated(size.ind)) > 0 || sum(duplicated(dims)) > 0 )
      stop('dimensions do not match unambiguously. Supply dims manually!')
  }
  if (is.element(class(logical.ind), c('matrix', 'array'))) {
    dims.size <- dim(logical.ind)
  } else {
    dims.size <- length(logical.ind)
  } 
  logical.ind             <- array(logical.ind, dim = dims.size)
  dims.all                <- setdiff(1:length(dims.datacube), dims)
  ind.matrix.choice       <- array(which(logical.ind, arr.ind = TRUE), 
                                   dim = c(sum(logical.ind), length(dims.size)))
  
  dims.all.expand         <- list()
  for (i in 1:length(dims.all))
    dims.all.expand[[i]] <- 1:dims.datacube[dims.all[i]]
  dims.all.grid           <-  as.matrix(do.call(expand.grid, dims.all.expand))
  
  dims.all.grid.exp       <- apply(dims.all.grid, 2, 
                                   function(x){rep(x, times = dim( ind.matrix.choice)[1])}) 
  ind.matrix.choice.exp   <- apply(ind.matrix.choice, 2, 
                                   function(x){rep(x, each = dim( dims.all.grid)[1])}) 

  
  ind.matrix.all          <- cbind(dims.all.grid.exp, ind.matrix.choice.exp)
  ind.matrix.ord          <- ind.matrix.all[, order(c(dims.all, dims))]
  
  colnames(ind.matrix.ord)    <- paste('dim', 1:length(dims.datacube), sep='')
  ##value<< integer index matrix which can be used to index datacube
  ind.matrix.ord
}

