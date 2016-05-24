.filter.channel.names <- function(spikes, ids) {
  ## Filter out some channel names.
  ## Keep only the channels mentioned in IDS.
  ## If the elements of IDS are numeric, they are assumed to be the
  ## indexes of the spike trains; otherwise, they are assumed to be the 
  ## names of cells.
  ## e.g.
  ## spikes2 <- .filter.channel.names(spikes, c('-', 'g4a', 'a6a'))
  ## spikes2 <- .filter.channel.names(spikes, c('g4a', 'a6a'))
  ## spikes2 <- .filter.channel.names(spikes, c(5, 3, 1))
  ## first call throws away two channels; second call keeps just two channels.
  ## third just keeps the three trains mentioned.
  
  if (any(is.character(ids)))
    ids = .names.to.indexes(names(spikes), ids)
  
  spikes[ids]
}

.remove.empty.channels <- function(spikes) {
  ## Remove any spike trains that are empty, i.e. zero spikes in them.
  ## This can happen if the beg, end range is too narrow, or if a
  ## datafile is empty, which happens sometime for the Feller data.
  ## TODO: currentlly only used by the Feller reader, perhaps it could
  ## also be used for other routines too?
  
  nspikes <- sapply(spikes, length)
  empty <- which(nspikes==0)
  if ( any(empty) ) {
    spikes <- spikes[-empty]
  }
  spikes 
}

.names.to.indexes <- function(names, elems, allow.na=FALSE, allow.regex=TRUE) {
  ## Return the indexes of where each element of ELEMS is within NAMES.
  ## If the first element of ELEMS is '-', then return all indexes except
  ## those matching ELEMS.  If ELEMS is NULL, then 1:n is returned, where n is
  ## the length of NAMES.
  ## Example:
  ## names = c('a', 'b', 'c', 'd', 'e')
  ## .names.to.indexes(names, c('d', 'b', 'a'))  ## 4 2 1
  ## .names.to.indexes(names, c( '-', 'c', 'a')) ## 2 4 5
  ## .names.to.indexes(names, NULL)
  
  ## to check if first element is "-", we have to use this more
  ## complex expression, as elems[1] == "-" is an error if the first element
  ## by chance is NA.
  if (is.null(elems)) {
    return (1:length(names))
  }
  if ( isTRUE(all.equal("-", elems[1])) ) {
    invert = TRUE
    elems = elems[-1]
  } else {
    invert = FALSE
    
  }
  
  indexes = match(elems, names)
  
  
  if (allow.regex) {
    ## see if any elements returned NA, in which case try them individually
    ## as regular expressions.
    which.na <- which(is.na(indexes))
    if (any(which.na)) {
      regex.elems <- elems[which.na]
      new.indexes <- lapply(regex.elems, function(r) {grep(r, names)})
      new.indexes <- unique(unlist(new.indexes))
      indexes <- indexes[-which.na]
      indexes <- c(indexes, new.indexes) #TODO, preserve order?
    }
    allow.na <- TRUE                    #allow NAs through now.
  }
  
  if (!allow.na) {
    if (any(is.na(indexes)))
      stop('some indexes not found.')
  }
  
  if (invert)
    indexes = setdiff(1:(length(names)), indexes)
  
  indexes
  
}

get.array.info <- function(data) {
  ## Array-specific information; maybe this could go in a file, rather
  ## than be read-in separately.  Useful for the HDF5 functions.
  
  pos <- data$epos;  rownames(pos) <- data$names
  array <- data$array

  if ( any(grep('^Axion', array))) {
    ## e.g. Neurotox ongoing project.
    xlim <- c(0, 8000)
    ylim <- c(0, 6000)
    spacing <- 200
    corr.breaks <-  0                   #TODO; by default, do no breaks!
  }

  ## HACK for HDF5: the arrayname (stored in array) is a string, which
  ## can either be stored as a character vector or as an array.  R was
  ## getting confused about this, and all.equal() was failing...
  ## Need a clearer example...
  
  array <- as.character(array)
  layout <- list(xlim = xlim, ylim = ylim, spacing = spacing, 
                 pos = pos, array=array)
  class(layout) <- "mealayout"
  list(layout=layout, corr.breaks=corr.breaks)
}
