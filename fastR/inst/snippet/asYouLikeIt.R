apropos("^as\\.")[1:10]      # just a small sample
# convert numbers to strings (this drops attributes)
as.character(xm)             
# convert matrix to vector
as.vector(xm)
as.logical(xm)
alpha <- c("a","1","b","0.5")    
mode(alpha)
as.numeric(alpha)      # can't do the coersion, so NAs are introduced
as.integer(alpha)      # notice coersion of 0.5 to 0
