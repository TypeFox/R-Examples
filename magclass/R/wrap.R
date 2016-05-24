wrap <- function (x, map = list(NA), sep = ".") {
    if (!is.array(x) && !is.matrix(x)) 
        stop("Argument 'x' is not an array or a matrix: ", class(x)[1])
    if (!is.list(map)) 
        stop("Argument 'map' is not a list: ", class(map)[1])
    umap <- unlist(map)
    if (any(duplicated(umap))) {
        stop("Argument 'map' contains duplicated dimension indices: ", 
            paste(umap[duplicated(umap)], collapse = ", "))
    }
    dim <- dim(x)
    ndims <- length(dim)
    missingDims <- setdiff(1:ndims, umap)
    if (length(missingDims) > 0) {
        wildcard <- is.na(map)
        if (any(wildcard)) {
            map[[which(wildcard)]] <- missingDims
            umap <- unlist(map)
        }
        else {
            stop("Argument 'map' miss some dimensions: ", paste(missingDims, 
                collapse = ", "))
        }
    }
    falseDims <- setdiff(umap, 1:ndims)
    if (length(falseDims) > 0) {
        stop("Argument 'map' contains non-existing dimensions: ", 
            paste(falseDims, collapse = ", "))
    }
    if (any(diff(umap) < 0)) {
        perm <- umap
        x <- aperm(x, perm = perm)
        map <- lapply(map, FUN = function(ii) match(ii, perm))
    }
    dim <- dim(x)
    dim2 <- lapply(map, FUN = function(ii) prod(dim[ii]))
    dimnames <- dimnames(x)
    
    tmp_dn<-function(map,dimnames) {
      dimnames2 <- list()
      nn <- NULL
      for(dim in 1:length(map)){
        names<-NULL
        for (ii in map[[dim]]) {
            if (is.null(names)) {
                names <- dimnames[[ii]]
                name_names<-names(dimnames)[ii]
            }
            else {
                names <- paste(names, rep(dimnames[[ii]], each = length(names)), 
                  sep = sep)
                name_names<-paste(name_names,names(dimnames)[ii],sep=sep)
            }
        }
        dimnames2[[dim]]<-names
        nn <- c(nn, name_names)
      }
      #Trick to set names even for NULL entries
      dimnames2[[dim+1]] <- "fake"
      names(dimnames2) <- c(nn,"fake")
      dimnames2[[dim+1]] <- NULL
      return(dimnames2)        
    }
    
    dim(x) <- dim2
    dimnames <- tmp_dn(map,dimnames)
    if(any(dim(x)==0)) {
      dimnames[dim(x)==0] <- NULL
    }
    dimnames(x) <- dimnames
    return(x)
}