# Check map for NAs, NaNs, etc
checkMap<-function(map, filter.map){
  if(!grepl('proj=utm|proj=longlat',proj4string(map))) stop("Projection needs to be in utm or longlat")

  if(grepl('+proj=utm.*',proj4string(map)))
    if(!grepl('+units=m',proj4string(map)))
      message('Assuming UTM +units=m')

  if(any(is.na(getValues(map)))) stop('NAs in habitat map. Replace with 0s')

  if(any(is.nan(getValues(map)))) stop('NaNs in habitat map')

  if(!is.null(filter.map)){
    if(proj4string(filter.map) != proj4string(map))
      stop('map and filter.map must have the same projection')
    if(any(is.nan(getValues(filter.map))))
      stop('NaNs in filter.map')
    if(any(is.na(getValues(map))))
      stop('NAs in filter.map. Replace with 0s')
    if(extent(filter.map)!=extent(map))
      stop('map and filter.map must have the same extent')
    if(any(res(filter.map)!=res(map)))
      stop('map and filter.map must have the same resolution')
    }

  return(map)
}


# Check parameter list for missing values, update names, etc
checkParameters<-function(pList,argList){

  if(is.null(pList$maxDistQ))
    pList$maxDistQ <- rep(1, length(pList$MFratio))

  if(is.null(pList$wghts))
    pList$wghts <- T

  if(is.null(pList$filter.cutoff)){
    if(!is.null(argList$filter.map)){
      if(all(getValues(argList$filter.map)<=1)){
        pList$filter.cutoff <- 0.95
      }
    }
  }

  if('trunk' %in% names(pList))
    stop("Truncation parameter has been reworked using probabilities instead of SDs.
      Use 'maxDistQ' instead of 'trunk'")

  oldnames <- c('howfar', 'howmuch', 'trunk', 'HRcenter.cutoff')
  newnames <- c('moveDist', 'moveDistQ', 'maxDistQ', 'habitat.cutoff')

  if(any(oldnames %in% names(pList))){
    old<-match(oldnames, names(pList))
    names(pList)[old] <- newnames[!is.na(old)]
    warning(paste('The following parameter names are deprecated:', 
                    oldnames[!is.na(old)],
                  '\n  Replace with: ', newnames[!is.na(old)]))
  }

  return(pList)
}