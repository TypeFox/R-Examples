"writeGPDascii" <-
function(x, counts = NULL, depths = NULL,
                          sample.names = NULL,
                          absolute.ages = NULL,
                          taxa = NULL,
                          short.names = NULL,
                          tax.cat = NULL,
                          metadata = NULL,
                          output.file = NULL, ...){

if(is.strat.column(x)){
  counts <- x$counts
  depths <- x$depths
  absolute.ages <- x$absolute.ages
  taxa <- x$taxa
  short.names <- x$short.names
  tax.cat <- x$tax.cat
  metadata <- x$metadata
}

levels <- length(depths)
ntaxa <- ncol(counts)

tax.vector <- 1:ncol(counts)

tax.vector <- paste(format(tax.vector, width = 5), ' ',
                    format(short.names, width = 8), ' ',
                    tax.cat, ' ', taxa, sep = '')

GPDascii.vector <- paste('#', metadata, '\n', collapse = '', sep = '')
GPDascii.vector <- c(GPDascii.vector,
                     paste('# Created on ', date(), ' by the function writeGPDascii()\n# in the R package stratigraph written by W. A. Green (walton.green@yale.edu)\n', sep = ''))
GPDascii.vector <- c(GPDascii.vector, paste(ntaxa, ', ', levels, '\n',
                                            sep = ''))
GPDascii.vector <- c(GPDascii.vector, paste(tax.vector, '\n', sep = ''))
for(i in 1:levels){
  GPDascii.vector <- c(GPDascii.vector,
                       paste(depths[i],
                             ', ',
                             sample.names[i],
                             ', ',
                             absolute.ages[i],
                             '\n', sep = ''))
  # warnings are turned off for the matrix call, because the behavior
  #  desired elicits a lot of irritating warnings
  old.warn.opt <- options()$warn
  level.matrix <- matrix(counts[i,], ncol = 10, byrow = TRUE)
  options(warn = old.warn.opt)
  level.matrix[nrow(level.matrix),(11 - ((nrow(level.matrix) *
                                          ncol(level.matrix)) -
                                          ntaxa)):10] <- NA
  for(j in 1:nrow(level.matrix)){
     GPDascii.vector <- c(GPDascii.vector,
                          paste(level.matrix[j,], sep = '',
                                collapse = '\t'))
     GPDascii.vector <- c(GPDascii.vector, '\n')
  }
  GPDascii.vector <- c(GPDascii.vector, '\n')
}

if(is.null(output.file)){
  cat(GPDascii.vector)
} else {
  sink(output.file)
  cat(GPDascii.vector)
  sink()
}

} # End of function

