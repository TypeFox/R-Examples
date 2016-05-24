mafft.merge <- function(subMSA, mafft.exe, quiet = TRUE){
  
  quiet <- ifelse(quiet, " --quiet ", " ")
  
  ## create subMSAtable
  n <- sapply(subMSA, nrow)
  subMSAtable <- vector(length = length(n))
  init <- 0
  for ( i in seq_along(n) ){
    nn <- 1:n[i] + init
    init <- max(nn)
    subMSAtable[i] <- paste(nn, collapse = " ")
  }
  write(subMSAtable, "subMSAtable")
  
  ## prepare input file
  
  subMSA <- lapply(subMSA, as.list)
  subMSA <- do.call(c, subMSA)
  names(subMSA) <- gsub("^.+[.]", "", names(subMSA))
  fns <- c("mafftin.fas", "mafftout.fas")
  write.fas(subMSA, fns[1])
  
  ## assemble call to MAFFT
  call.mafft <- paste(mafft.exe, " --merge subMSAtable", 
                      quiet, fns[1], " > ", fns[2], sep = "")
#   if ( os == "unix" ){
    system(call.mafft, intern = FALSE, ignore.stdout = FALSE)
    res <- length(scan(fns[2], what = "c", quiet = TRUE))
    if ( res != 0 ) res <- read.fas(fns[2])
#   }
  unlink(fns[file.exists(c("subMSAtable", fns))])
  return(res)
}
