import.fasta <- function (file, aa.to.upper = TRUE, gap.to.dash = TRUE) {

  if(missing(file)) {
    stop("file is missing")
  }

  #read as a vector of lines
  lines <- readLines(file)

  #localize sequence identifiers
  #and check fasta format
  loc <- grep(">", lines)
  if (length(loc) == 0)
      stop("file is not in fasta format")
 
  #get sequence identifiers
  id <- sub("^>(\\S+).*$","\\1", lines[loc])
  nb.seq <- length(id)

  #localize sequence pieces for each identifier
  start <- loc + 1
  end <- loc - 1  
  end <- c(end[-1], length(lines))

  seq <- sapply(seq_len(nb.seq), function(i) {paste(lines[start[i]:end[i]], collapse = "")})
  seq <- gsub("\\s", "", seq)

  #turn aa into upper case
  if (aa.to.upper)
    seq <- toupper(seq)

  #give a list of split sequences
  seq <- strsplit(seq, split = "")
  names(seq) <- id

  #turn gap into dash character
  if (gap.to.dash)
    seq <- lapply(seq, function (i) {i[is.gap(i)] <- "-"; return(i)})
  class (seq) <- c("align")
  return(seq)
}
