mbind <- function(...) {  
  inputs <- list(...)
  #Remove NULL elements from list
  for(i in length(inputs):1) if(is.null(inputs[[i]])) inputs[[i]] <- NULL
  
  if(length(inputs)==1 & is.list(inputs[[1]])) inputs <- inputs[[1]]
  regio <- NULL
  cells <- NULL
  elems <- NULL
  years <- NULL
  diffspat <- FALSE
  difftemp <- FALSE
  diffdata <- FALSE
  for(i in 1:length(inputs)) {
    if(!is.magpie(inputs[[i]])) stop("Inputs must all be MAgPIE-objects")
    if(is.null(dimnames(inputs[[i]])[[3]])) dimnames(inputs[[i]])[[3]] <- paste("dummydimname",1:ndata(inputs[[i]]),sep="")
    #Check which dimensions differ
    if(suppressWarnings(any(sort(dimnames(inputs[[1]])[[1]])!=sort(dimnames(inputs[[i]])[[1]])))) diffspat <- TRUE
    if(suppressWarnings(any(sort(dimnames(inputs[[1]])[[2]])!=sort(dimnames(inputs[[i]])[[2]])))) difftemp <- TRUE
    if(suppressWarnings(any(sort(dimnames(inputs[[1]])[[3]])!=sort(dimnames(inputs[[i]])[[3]])))) diffdata <- TRUE
    years <- c(years, getYears(inputs[[i]]))
    elems <- c(elems, getNames(inputs[[i]]))
    cells <- c(cells, getCells(inputs[[i]]))   
    if(!diffspat & ncells(inputs[[1]])>1) inputs[[i]] <- inputs[[i]][getCells(inputs[[1]]),,]
    if(!difftemp & nyears(inputs[[1]])>1) inputs[[i]] <- inputs[[i]][,getYears(inputs[[1]]),]
    if(!diffdata &  ndata(inputs[[1]])>1) inputs[[i]] <- inputs[[i]][,,getNames(inputs[[1]])]
  }
  
  if(!(length(grep(".",cells,fixed=TRUE)) %in% c(0,length(cells)))) stop("Mixture of regional (no cell numbers) and cellular (with cell numbers) data objects! Cannot handle this case!")
  
  if(diffspat & difftemp) stop("Cannot handle objects! Spatial as well as temporal dimensions differ!")      
  if(difftemp & diffdata) stop("Cannot handle objects! Data as well as temporal dimensions differ!")      
  if(diffdata & diffspat) stop("Cannot handle objects! Data as well as spatial dimensions differ!") 
  if(!diffspat) {
    
  }
  if(difftemp) {
    if(length(years)!=length(unique(years))) stop("Some years occur more than once! Cannot handle this case!")
    output <- new("magpie",abind::abind(inputs,along=2))
  } else if(diffspat){
    if(length(cells) != length(unique(cells))) stop("Some regions/cells occur more than once! Cannot handle this case!")
    output <- new("magpie",abind::abind(inputs,along=1))
  } else {
    output <- new("magpie",abind::abind(inputs,along=3))
  }
  if(length(grep("dummydimname",getNames(output),fixed=TRUE))==ndata(output)) dimnames(output)[[3]] <- NULL 
  names(dimnames(output)) <- names(dimnames(inputs[[1]]))
  return(output)
}  
