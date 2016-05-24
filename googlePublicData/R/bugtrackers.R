.checkSlices <- function(dims, by) {
################################################################################
# Checks if there's any pair of slices with the same dimentions (error)
################################################################################
  # Builds a list of dims by slice
  slices.dims <- do.call(list,by(dims, by, function(x) unlist(x[order(x)])))
  
  # Generates a list of cohorts (unique) that shouldn't be store more than once
  unique.dims <- unique(slices.dims)
  
  for (i in 1:length(unique.dims)) {
    # Sets the counters
    counter <- 0
    fields <- NULL
    for (j in names(slices.dims)) {
      test <- identical(slices.dims[[j]], unique.dims[[i]])
      
      if (test) {counter <- counter + 1; fields <- paste(fields, j, sep=',')}
      if (counter > 1) {
        stop("Dimention(s) ", unlist(unique.dims[[i]]),
             " apear more than once in the collection at ", fields,
             ". Variables in those tables should be grouped in only one table.")
      }
    }
  }
}

.checkDuplConcepts <- function(concepts) {
################################################################################
# Checks if there's any concepts duplicated as a result of multiple data types
# In the case of beeing all numeric, DSPL assumes the minumum common (float) and
# fixes the error. Output = Warning
################################################################################  
  
  # Partial fix (should work, need to see
  # http://stackoverflow.com/questions/23475309/in-r-is-it-possible-to-suppress-note-no-visible-binding-for-global-variable)
  # type <- is.dim.tab <- id <- label <- freq <- NULL
  
  concepts2 <- unique(
    # subset(concepts,subset=type != 'date' & is.dim.tab==F, select=c(id, label, type))
    concepts[
      with(concepts, type!='date' & is.dim.tab == FALSE),
      colnames(concepts) %in% c('id','label','type'),FALSE]
    )
  
  # Frequency table
  freq.tab <- as.data.frame(table(concepts2$id), stringsAsFactors=F)
  
  colnames(freq.tab) <- c('id','freq')  
  dpl.concepts <- freq.tab[freq.tab[["freq"]] > 1,,FALSE]
  
  # Number of duplicated concepts
  ndpl.concepts <- nrow(dpl.concepts)
  
  # If there are any dpl concepts
  if (ndpl.concepts > 0) {
    
    # Loop for each and one of the duplicated concepts
    for (dpl in dpl.concepts$id) {
      
      # Testing if all the data types of the dpl concepts is numeric
      test <- all(concepts$type[concepts$id == dpl] %in% c('float', 'integer'))
      
      # Fixing the concept type
      if (test) {
        touse <- which(concepts$id == dpl)
        concepts$type[touse] <- "float"
        warning(dpl,' concept was fixed at slices: \n - ',
                paste(unique(concepts$slice[touse]), collapse='\n - '))
      }
      else {
        stop('Duplicated concepts cannot be homogenized\n',dpl,
             concepts[concepts$id == dpl, c('id','type')])
      }
    }    
  }
  return(concepts)
}

.checkPath <- function(x, type='output') {
################################################################################
# DESCRIPTION:
# Checks if the output path exists, otherwise stops the routine.
################################################################################  
  if (!is.na(x)) {
    ER <- try(file.exists(x), silent=T)
    if (class(ER) == 'try-error') {
      stop('Incorrect ', type,' path:\n\t\t\t', x, '\n\t\t couldn\'t be found')
    } 
  }
}