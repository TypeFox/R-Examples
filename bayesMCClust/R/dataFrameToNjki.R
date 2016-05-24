dataFrameToNjki <-
function(dataFrame) {     # dataFrame ... matrix or data.frame with time series of equal length in rows

    # class(dataFrame) has to be "matrix" or "data.frame" 
    stopifnot( class(dataFrame)=="matrix" | class(dataFrame)=="data.frame" )
    
    N <- dim(dataFrame)[1]       # number of individuals/units/objects 
    T <- dim(dataFrame)[2]       # length of time series --> equal for all individuals!! and corresponds to (T-1) transitions!!  
    
    minCat <- min(dataFrame)     # store indicator for lowest/first category (consecutive numbering necessary!)
    
    K <- max(dataFrame) - minCat  # K+1 categ. (consecutive numbers): min(dataFrame), ..., max(dataFrame) to be transformed to 0,...,K (see below)
    
    # initialise 3-dim array, in which each individual transition matrix (absolute frequencies) of each unit i is to be stored
    Njk.i <- array(0, c(K + 1, K + 1, N), dimnames = list( as.numeric(names(table(dataFrame))), as.numeric(names(table(dataFrame))), NULL) ) 
    
    # define function to update transition matrix Njk.i[,,i] with each single transition of individual i (see below)
    counttr.i <- function(tr1, tr2) {
        Njk.i[tr1 + 1, tr2 + 1, i] <<- Njk.i[tr1 + 1, tr2 + 1, i] + 1  # ( Njk.i has to exist outside this function ) 
    }
    
    # if required transform numbers of categories to 0,...,K!!
    dataM <- dataFrame # as.matrix(dataFrame) - minCat
    
    # count and store absolute frequencies of (T-1) transitions of each single individual into Njk.i
    for( i in 1:N ) {  
        mapply( counttr.i, dataM[i, 1:(T-1)], dataM[i, 2:T] ) # count and store the abs. frequencies of transitions of individual i into Njk.i
    }
    
    return( invisible( Njk.i ) )
}
