export.fasta <- function (x, outfile = "alignment.fa", ncol = 60, open = "w") {

  testClass<-0
  if (inherits(x,"list")) {
    testClass<-1
    nbClus <- length(x)
    id <- names(x)
    if (!inherits(x[[1]], "align"))
      testClass<-0
  }
  if (inherits(x, "align")) {
    testClass<-2
    nbClus <- 1
  }
  if (testClass == 0)
    stop("object of class 'align' or named list of 'class' expected")
  for(i in 1:nbClus) {
    if (testClass == 2)
      clus <- x
    else
      clus <- x[[i]]
    if (length(outfile) == nbClus)
      out <- outfile[i]
    else {
      if (nbClus == 1)
	out <- outfile[1]
      else
	out <- paste(c(id[i],"_",outfile[1]),collapse="",sep="")
    }
      
    file <- file(description = out, open = open)
    nb.seq <- length(clus)
    name <- names(clus)
    ln.seq <- length(clus[[1]])
    j <- round(ln.seq/ncol)
    supJ <- ln.seq - (j * ncol)
    for (i in 1:nb.seq) {
	writeLines(paste(">", name[i], sep = ""), file)
	sapply(c(1:j),function(y) writeLines(paste(clus[[i]][(((y-1)*ncol)+1):(y*ncol)],sep="",collapse=""),file))
	if (supJ >= 1) 
	  writeLines(paste(clus[[i]][((j*ncol)+1):ln.seq],sep="",collapse=""),file)
    }
    close(file)
  }
}

