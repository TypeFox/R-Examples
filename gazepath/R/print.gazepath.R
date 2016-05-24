## Make a print method for class gazepath
print.gazepath <- function(x, ...){
  miss <- length(which(is.na(unlist(x[[2]])))) / length(unlist(x[[2]])) * 100
  fix <- length(summary(x)[,1])
  raw <- length(unlist(x[[1]]))
  cat(' Fixations were estimated using the', x[[4]], 'method',
      '\n', '\n',
      round(miss), 'percent of the data points is missing',
      '\n',
      raw, 'raw data points were classifies into', fix, 'fixations',
      '\n',
      'The plot(gazepath-object, i = trial-number) function can be used to visualize the results',
      '\n',
      'The summary(gazepath-object) can be used get all fixations')
}
