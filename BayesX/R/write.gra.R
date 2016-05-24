write.gra <- function(map, file, replace=FALSE)
{
    if(! inherits(map,"gra"))
        stop("Argument 'map' is not an object of class 'gra'!")

    ## check whether the file exists
    if(replace & file.exists(file))
        test <- file.remove(file)
    if(!replace & file.exists(file))
        stop("Specified file already exists!")

    ## names of districts
    districts <- as.integer(rownames(map))

    ## no. of regions
    S <- length(districts)
    write(S, file)
    
    ## loop over the regions
    for(i in 1:S){
        ## write name of the district
        write(districts[i], file, append=TRUE)
        ## write no. of neighbors
        write(map[i,i], file, append=TRUE)
        
        ## derive and write neighbors
        ind <- which(map[i,]==-1)-1
        write(ind, file, ncolumns=length(ind), append=TRUE)
    }
    
    return(invisible())
}
