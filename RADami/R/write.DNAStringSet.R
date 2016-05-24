write.DNAStringSet <-
function(x, format= 'phylip', padding = 30, filename = 'DNAStringSetOut.phy') {
  # writes a sequence matrix to phylip format
  x.width <- width(x)[1]
  x <- as.character(x)
  for(i in 1:length(x)) x[i] <- paste(names(x)[i], paste(rep(" ", (padding - nchar(names(x)[i]))), collapse = ''), x[i], sep = '')
  writeLines(c(paste(length(x), x.width), x), filename)
  return(0)
  }
