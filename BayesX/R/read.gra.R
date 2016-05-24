read.gra <- function(file, sorted=FALSE)
{
    ## disable warnings and revert to old settings afterwards
    oldOptions <- options(warn = -1)
    on.exit(options(oldOptions))    
    
    ## read information from the graph file
    mapinfo <- scan(file)

    ## extract the number of districts
    S <- mapinfo[1]
    mapinfo <- mapinfo[-1]
    cat("Note: map consists of", S,"regions\n")

    cat("Reading map ...")

    ## create a vector to contain the names of the districts and the adjacency matrix
    districts <- rep(0,S)
    pmat <- matrix(0,S,S)

    ## loop over the districts
    for(i in 1:S){
        ## extract the name of the i-th district
        districts[i] <- paste(mapinfo[1])
        mapinfo <- mapinfo[-1]

        ## extract the number of neighbors for the i-tz district
        nrn <- mapinfo[1]
        mapinfo <- mapinfo[-1]
        pmat[i,i] <- nrn

        ## off-diagonal entries in pmat
        ## Note: the indices are zero-based!
        if(nrn > 0){
            for(j in 1:nrn){
                pmat[i,mapinfo[j]+1] <- pmat[mapinfo[j]+1,i] <- -1
            }
            mapinfo <- mapinfo[-(1:nrn)]
        }
    }

    cat(" finished\n")
    
    if(sorted){
        ## sort districts and adjacency matrix
        if(sum(is.na(as.numeric(districts))) == 0){
            districts <- as.numeric(districts)
            cat("Note: regions sorted by number\n") 
        }
        else 
            cat("Note: regions sorted by name\n")
        ord <- order(districts)
        districts <- sort(districts)
        
        pmat.sort <- matrix(0,S,S)

        for(i in 1:S){
            ordi <- ord[i]
            for(j in 1:S){
                ordj <- ord[j]
                pmat.sort[i,j] <- pmat[ordi, ordj]
            } 
        }
        
        pmat <- pmat.sort
    }   

    rownames(pmat) <- districts
    colnames(pmat) <- districts

    class(pmat) <- "gra"
    return(pmat)
}

