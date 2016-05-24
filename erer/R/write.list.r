write.list <- function(z, file, t.name = NULL, row.names = FALSE, ...) {
  if (!inherits(z, "list")) stop("\nNeed an 'list' object.\n")

  # create table names
  if (is.null(t.name)) {
    if (is.null(names(z))) {
      add.name <- paste("result.", 1:length(z), sep='')
    } else {
      add.name <- names(z)
    }
  } else {
    if (length(t.name) != length(z)) {
      stop("\n 't.name' and 'z' should have the same length.\n")
    } else {    
      add.name <- t.name
    }
  }
  
  # write.table for each list
  options(warn = - 1) # suppressWarnings    
  for (k in 1:length(z)) {
    dat <- as.data.frame(z[[k]])
    if (row.names) {
      h2 <- as.data.frame(cbind(Result = add.name[k], 
        row.name = rownames(dat), dat))
    } else {
      h2 <- as.data.frame(cbind(Result = add.name[k], dat))    
    }
    h3 <- rbind(apply(h2, 2, as.character), "") # add a blank row to each table
    ap <- ifelse(k==1, FALSE, TRUE)
    write.table(x=h3, file=file, sep=",", append=ap, row.names=FALSE, ...)
  }
  options(warn = 0) 
}