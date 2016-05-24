### The function "mlData" reads in the orginal examinee x items
### matrix, and generate data for maximum likehood estimateion of Log
### Multiplicative Association Models
###
### it is a general function that can be used for multiple latent
### variables
###
### Arguments:
###   data : the orginal examinee response matrix (nexaminee x nitem)
###   item.mtx : adjacency matrix between items and latent traits
###   trait.mtx: adjacency matrix among latent traits
###
### Value:
###   A matrix of the stacked data is returned. The columns are
###   'count': number of examinees that has the item response pattern
###   'item..': item response pattern 
###   'phi..': corresponding to the parameters in LMA
###
### Notes:
###   for simplification, we use 'raw' scores for the repsonpse
### data, score matrix and data maxtrix are the same( 'data' in
### arguments)
###
### Author: Zhushan "Mandy" Li
### Last Updated: August 20, 2006

mlData <- function(data, item.mtx, trait.mtx){
  ## make sure the arguments are matrices
  if(! is.matrix(data)) data <- as.matrix(data)
  if(! is.matrix(item.mtx)) item.mtx <- as.matrix(item.mtx)
  if(! is.matrix(trait.mtx)) trait.mtx <- as.matrix(trait.mtx)

  
  nexaminee <- nrow(data)  # number of examinees
  nitem <- ncol(data)  # number of items
  ntrait <- ncol(item.mtx)  # number of latent trait
  
  if(nitem != nrow(item.mtx)){
    stop("Error: the number of items in data matrix is not consistent with
the number of rows in item matrix" )
  }

  ## lower trianglular of the latent traits adjacency matrix
  trait.low <- trait.mtx & lower.tri(trait.mtx,diag=TRUE)
  ## total number of parameters for latent variables
  nphi <- sum(trait.low)

  which.low <- which(trait.low, arr.ind=TRUE)
  ## the names of the parameters
  phi.names <- paste("phi",apply(which.low, 1, paste, collapse=''), sep='')

  ## construct matricies used for calculating phi
  phi.array <- array(0, dim=c(ntrait, ntrait, nphi))
  phi.array[cbind(which.low, 1:nphi)] <- 1
  phi.array[cbind(which.low[,2:1], 1:nphi)] <- 1

  jjmtx <- array(dim=c(nitem, nitem, nphi))
  for (i in 1:nphi){
    jjmtx[,,i] <- item.mtx %*% phi.array[,,i] %*% t(item.mtx)
    diag(jjmtx[,,i]) <- 0
  }

  
  ctdata <- countData(data)
  
  ## caluculate phi
  phi <- matrix(NA, nrow=nrow(ctdata), ncol=nphi)
  colnames(phi) <- phi.names
  for(i in 1:nphi){
    jjmtx.1.idx <- which( (jjmtx[,,i] != 0) & lower.tri(jjmtx[,,i]), arr.ind=TRUE)
    phi[,i] <- apply(ctdata[,-1], 1,
        function(x)crossprod(x[jjmtx.1.idx[,1]], x[jjmtx.1.idx[,2]]) )
  }

  newdata <- cbind(ctdata, phi)
  return(newdata);
}  

### The function 'countData' enumerates all possible response patterns and then counts
### the number of occurance of each response pattern in the data
###
### Syntax
###   countData(data, ncat, nitem)
###
### Arguments
###   'data': the orginal examinee response matrix (nexaminee x nitem)
###   'ncat': number of choices of response for one item
###   'nitem': number of items in the data
###
### Value
###   The returned value is a data frame with first column denoting the counts,
### and the rest columns are response patterns.
###
### Notes
###   The number of all possible response patterns is ncat^nitem and this number
### becomes huge if the number of items is large. Computer may not have enough memory
### to run.
###
### Author: Zhushan "Mandy" Li
### Last Updated: Aug 23, 2006
countData <- function(data, ncat=length(unique(as.matrix(data[!is.na(data)]))), nitem=ncol(data)){
  ### Transform the data format
  ## itemcode is a table of the all possible combinations of item repsonses.
  ##   it uses function integer.base.b, which changes a interger into
  ## a vector of base b representation 
  itemcode <- data.frame(integer.base.b(0:(ncat^nitem-1), b=ncat, ndigits=nitem))
  if(!is.null(colnames(data))){
    colnames(itemcode) <- colnames(data)
  }
     
  count <- rep(0, ncat^nitem);
  nexaminee <- nrow(data);
  for(i in 1:nexaminee){
    ## the function dec.base.b changes a vector of base b representation of a integer into its decimal integer form.
    pattern.idx <- dec.base.b(data[i,], base=ncat) + 1;
    count[pattern.idx] <- count[pattern.idx] + 1;
  }
  
  return(data.frame(count=count, itemcode))
}

### The function 'integer.base.b'
### changes a interger into a vector of base b representation 
### this function is called by 'countData'
integer.base.b <- function(x, b=2, ndigits=NULL){
  xi <- as.integer(x)
  if(any(is.na(xi) | ((x-xi)!=0)))
    print(list(ERROR="x not integer", x=x))
  N <- length(x)
  xMax <- max(x)
  if(is.null(ndigits)){
    ndigits <- (floor(logb(xMax, base=b))+1)
  }
  Base.b <- array(NA, dim=c(N, ndigits))
  for(i in 1:ndigits){#i <- 1
    Base.b[, ndigits-i+1] <- (x %% b)
    x <- (x %/% b)
  }
  if(N ==1) Base.b[1, ] else Base.b
}

### the function 'dec.base.b' changes a vector of base b representation of
### a integer into its decimal integer form.
### called by countData
dec.base.b <- function( x, base=2 ){
  basis <- base^((length(x)-1):0);
  return( x%*%basis )
}

### Instead of enumerate all possible response patterns, 'countData2' only extract
### the unique response patterns in the data and gives the number of counts of each
### pattern in the data.
###
### Notes
###   This function has not been called by any other functions yet
###
### Author: Zhushan "Mandy" Li
### Last Updated: Aug 23, 2006
countData2 <- function(data){
  dup.ind <- duplicated(data)
  uniquePattern <- as.data.frame(data[!dup.ind,])
  count <- rep(1, nrow(uniquePattern))
  dup.match <- match(do.call(paste, c(as.data.frame(data[dup.ind,]), sep="\r")),
                     do.call(paste, c(uniquePattern, sep="\r")))
  for(i in dup.match){
    count[i] <- count[i] + 1;
  }
  return(data.frame(count=count, uniquePattern))
}

