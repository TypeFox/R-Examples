import.msf <- function (file, aa.to.upper = TRUE, gap.to.dash = TRUE) {
  
  if(missing(file)) {
    stop("file is missing")
  }

  #read as a vector of lines
  lines <- readLines(file)

  #check msf format
  check1 <- grep("MSF:.*Type:.*Check:.*", lines)
  check2 <- grep("Name:.*Len:.*Check:.*Weight:.*", lines)
  limit <- grep("^//", lines)
  if (length(check1) == 0 || length(check2) == 0 || length(limit) == 0)
      stop("file is not in msf format")

  #get sequence identifiers from header
  id.head <- sub("^\\s*Name:\\s+(\\S+).*$","\\1", lines[check2])
  nb.seq <- length(id.head)

  #check duplicated identifiers in header
  if(any(duplicated(id.head)))
    stop("duplicated identifiers in header")

  #get sequence identifiers from alignment
  align <- grep("^\\s*\\S+\\s+[^1-9]+$", lines[limit:length(lines)], value = TRUE)
  id.align <- sub("^\\s*(\\S+)\\s+[^1-9]+$", "\\1", align)

  #localize sequence pieces for each identifier
  loc <- lapply(seq_len(nb.seq), function(i) {which(id.align == id.head[i])})

  #paste and clean sequences
  seq <- sapply(loc, function(i) {paste(sub("^\\s*\\S+\\s+([^1-9]+)$", "\\1", align[i]), collapse = "")})
  seq <- gsub("\\s", "", seq)
  
  #turn aa into upper case
  if (aa.to.upper)
    seq <- toupper(seq)

  #give a list of split sequences
  seq <- strsplit(seq, split = "")
  names(seq) <- id.head

  #turn gap into dash character
  if (gap.to.dash)
    seq <- lapply(seq, function (i) {i[is.gap(i)] <- "-"; return(i)})
  class (seq) <- c("align")
  return(seq)
}