############################################################################################################################

# Function: CreateTableSampleSize .
# Argument: data.strucure and label (optional).
# Description: Generate a summary table of sample size for the report.

CreateTableSampleSize = function(data.structure, label = NULL) {

  # Number of sample ID
  n.id <- length(data.structure$id)
  id.label = c(unlist(lapply(lapply(data.structure$id, unlist), paste0, collapse = ", ")))

  if (!any(is.na(data.structure$sample.size.set))){
    # Number of sample size
    n.sample.size = nrow(data.structure$sample.size.set)

    # Label
    if (is.null(label)) label = paste0("Sample size ", 1:n.sample.size)
    else label = unlist(label)
    if (length(label) != n.sample.size)
      stop("Summary: Number of the sample size labels must be equal to the number of sample size sets.")

    # Summary table
    sample.size.table <- matrix(nrow = n.id*n.sample.size, ncol = 4)
    ind <-1
    for (i in 1:n.sample.size) {
      for (j in 1:n.id) {
        sample.size.table[ind, 1] = i
        sample.size.table[ind, 2] = label[i]
        sample.size.table[ind, 3] = id.label[j]
        sample.size.table[ind, 4] = data.structure$sample.size.set[i,j]
        ind <- ind+1
      }
    }

    sample.size.table = as.data.frame(sample.size.table)
    colnames(sample.size.table) = c("sample.size","Sample size set", "Sample", "Size")
  }
  else if (!any(is.na(data.structure$event.set))){
    # Number of sample size
    n.events = nrow(data.structure$event.set)

    # Label
    if (is.null(label)) label = paste0("Event ", 1:n.events)
    else label = unlist(label)
    if (length(label) != n.events)
      stop("Summary: Number of the events labels must be equal to the number of events sets.")

    # Summary table
    sample.size.table <- matrix(nrow = n.events, ncol = 3)
    ind <-1
    for (i in 1:n.events) {
      sample.size.table[i, 1] = i
      sample.size.table[i, 2] = label[i]
      sample.size.table[i, 3] = data.structure$event[i,1]
    }

    sample.size.table = as.data.frame(sample.size.table)
    colnames(sample.size.table) = c("sample.size","Event set", "Total number of events")
  }


  return(sample.size.table)
}
# End of CreateTableSampleSize