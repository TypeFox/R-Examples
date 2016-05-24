dataListToNjki <-
function(dataList) {     # dataList ... list of time series with individual length

    # class(dataList) has to be "list"
    stopifnot( class(dataList)=="list" )
    
    N <- length(dataList)             # number of individuals/units/objects  
    
    minCat <- range(dataList)[1]      # store indicator for lowest/first category  (consecutive numbering necessary!)
    
    K <- range(dataList)[2] - minCat  # K+1 categ. (consecutive numbers): min(dataList), ..., max(dataList)
    
    # initialise 3-dim array, in which each individual transition matrix (absolute frequencies) of each unit i is to be stored
    Njk.i <- array(0, c(K + 1, K + 1, N), dimnames = list( as.numeric( names( table( unlist(dataList) ) ) ), as.numeric( names( table( unlist(dataList) ) ) ), NULL) ) 
    
    # define function to update transition matrix Njk.i[,,i] with each single transition of individual i (see below)
    counttr.i <- function(tr1, tr2) {
        Njk.i[tr1 + 1, tr2 + 1, i] <<- Njk.i[tr1 + 1, tr2 + 1, i] + 1  # ( Njk.i has to exist outside this function ) 
    }
    
    # if required transform indicators of categories to 0,...,K!!
    transDataList <- dataList # lapply(dataList, function(x) x - minCat)   
    
    lengthTS <- sapply(dataList, length)    # store time series length
    
    # count and store absolute frequencies of (T-1) transitions of each single individual into Njk.i
    for( i in 1:N ) {     # count and store the abs. frequencies of transitions of individual i into Njk.i
        mapply( counttr.i, transDataList[[i]][ 1:(lengthTS[i] - 1 ) ], transDataList[[i]][ 2:lengthTS[i] ] )
    }
    
    return( invisible( Njk.i ) )
}
