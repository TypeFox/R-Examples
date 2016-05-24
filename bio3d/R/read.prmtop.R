"read.prmtop" <- function(file) {
  cl <- match.call()

  if (missing(file)) {
    stop("read.prmtop: please specify a prmtop 'file' for reading")
  }

  toread <- file.exists(file)
  if (!toread) {
    stop("No input prmtop file found: check filename")
  }

  readformat <- function(s) {
    s=trim(s)
    s=substring(s, 9, nchar(s)-1)
    tmp <- unlist(strsplit(s, ""))
    hmm <- c("a", "E", "I")
    type <- hmm[which(hmm %in% tmp)]

    dims <- unlist(strsplit(s, type))
    if(type=="E")
      dims=c(dims[1], unlist(strsplit(dims[2], "\\.")))
    else
      dims=c(dims, NA)

    ## records per line
    ## record length
    ## type
    return(c(dims, type))
  }


  trim <- function(s) {
    s <- sub("^ +", "", s)
    s <- sub(" +$", "", s)
    s[(s == "")] <- NA
    s
  }

  split.line <- function(x) {
    tmp <- unlist(strsplit(x, split=" "))
    inds <- which(tmp!="")
    return(tmp[inds])
  }

  ## Read and parse file
  raw.lines <- readLines(file)
  flags.ind <- grep("%FLAG", raw.lines)
  
  parse.line <- function(line, fmt) {
    tmp <- seq(1, as.numeric(fmt[1])*as.numeric(fmt[2]), by=as.numeric(fmt[2]))
    substring(line, tmp, c(tmp[2:length(tmp)]-1, nchar(line)))
  }


  all.data <- list()
  for(i in 1:length(flags.ind)) {

    ind.start <- flags.ind[i]
    if(i==length(flags.ind))
      ind.end <- length(raw.lines)
    else
      ind.end <- flags.ind[i+1] - 1
    
    flagname <- split.line(trim( raw.lines[ind.start] ))[2]
    tmp.lines <- raw.lines[ind.start:ind.end][-c(1,2)]

    if(flagname=="TITLE")
      fmt <- c(1, 20, NA, 'a')
    else
      fmt <- readformat(raw.lines[ind.start+1])

    tmp.lines <- trim(unlist(lapply(tmp.lines, parse.line, fmt)))
    
    if(fmt[4]=="I" || fmt[4]=="E")
      tmp.lines=as.numeric(tmp.lines)
  
    tmp.lines=tmp.lines[!is.na(tmp.lines)]
    all.data[[flagname]]=tmp.lines
  }

  all.data$call=cl
  class(all.data)=c("amber", "prmtop")
  return(all.data)
}

