mbind2 <- function(...) { 
 isnull   <- sapply(list(...),is.null)
 ismagpie <- sapply(list(...),is.magpie)
 if(any(!ismagpie & !isnull)) stop("Input(s) ",paste(which(!ismagpie & !isnull),collapse=", ")," is/are no MAgPIE objects(s)!")
 if(length(list(...))==1 & is.list(list(...)[[1]])) stop("Input is a list!") #will be stopped earlier because of magpie-check ... inputs <- inputs[[1]]

 if(sum(ismagpie)==1) return(list(...)[[which(ismagpie)]]) #exactly one MAgPIE object
 if(sum(ismagpie)==0) return(NULL) #only NULL values
 #compare size of MAgPIE objects
 diffdata  <- !(length(unique(lapply(list(...)[ismagpie],getNames)))==1)
 difftemp  <- !(length(unique(lapply(list(...)[ismagpie],getYears)))==1)
 .tmp <- function(x) return(getCells(clean_magpie(x)))
 diffcells <- !(length(unique(lapply(list(...)[ismagpie],.tmp)))==1)
 if(diffcells) stop("Cannot handle different spatial dimensions!")
 if(diffdata & difftemp) stop("Cannot handle objects! Data as well as temporal dimensions differ!")      
 if(difftemp) {
   years <- unlist(lapply(list(...)[ismagpie],getYears))
   nyears <- sapply(list(...)[ismagpie],nyears)
   if(sum(nyears)!=length(years)) stop("Combining MAgPIE objects with and without years is not possible!")
   if(length(years)!=length(unique(years))) stop("Some years occur more than once! Cannot handle this case!")
    output <- new.magpie(getCells(list(...)[ismagpie][[1]]),years,getNames(list(...)[ismagpie][[1]])) 
    for(i in 1:length(list(...)[ismagpie])) output[,getYears(list(...)[ismagpie][[i]]),] <- list(...)[ismagpie][[i]]
    return(output)
  } else {
    elems <- unlist(lapply(list(...)[ismagpie],getNames))
    nelem <- sapply(list(...)[ismagpie],ndata)
    if(sum(nelem)!=length(elems)) stop("Combining MAgPIE objects with and without data names is not possible!")
    if(length(elems)!=length(unique(elems))) {
      #duplicates exist -> make dimnames unique
      lelems <- lapply(list(...)[ismagpie],getNames)
      for(i in 1:length(lelems)) lelems[[i]] <- paste(lelems[[i]],i,sep=".x")
      elems <- unlist(lelems)
      output <- new.magpie(getCells(list(...)[ismagpie][[1]]),getYears(list(...)[ismagpie][[1]]),elems)
      for(i in 1:length(list(...)[ismagpie])) {
        getNames <- paste(getNames(list(...)[ismagpie][[i]]),i,sep=".x")
        output[,,getNames] <- setNames(list(...)[ismagpie][[i]],getNames)
      }
    } else {
      output <- new.magpie(getCells(list(...)[ismagpie][[1]]),getYears(list(...)[ismagpie][[1]]),elems)
      names(dimnames(output)) <- names(dimnames(list(...)[ismagpie][[1]]))
      for(i in 1:length(list(...)[ismagpie])) output[,,getNames(list(...)[ismagpie][[i]])] <- list(...)[ismagpie][[i]]      
    }
    return(output)
  }
}  
