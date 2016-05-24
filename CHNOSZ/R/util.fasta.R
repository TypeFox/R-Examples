# CHNOSZ/util.fasta.R
# read and manipulate FASTA sequence files

grep.file <- function(file,pattern="",y=NULL,ignore.case=TRUE,startswith=">",lines=NULL,grep="grep") {
  # return the line numbers of the file that contain
  # the search term x and optionally don't contain y
  sysgrep <- function(i) {
    # 20091021 changed grep to egrep
    sysexp <- paste(mycat,' "',file,'" | ',grep,' -n ',ic,' "',startswith,pattern[i],'"  | cut -f 1 -d ":"',sep="")
    ix <- as.integer(system(sysexp,intern=TRUE))
    return(ix)
  }
  Rgrep <- function(i) {
    ix <- grep(pattern[i],lines,ignore.case=ignore.case)
    if(!is.null(y)) {
      iy <- grep(y,lines,ignore.case=ignore.case)
      ix <- ix[!ix %in% iy] 
    }
    if(!is.null(startswith)) {
      ihead <- which(substr(lines,1,1)==startswith)
      ix <- ix[ix %in% ihead]
    }
    return(ix)
  }
  # use the system's grep if available and y is NULL
  # TODO: include other *nix
  if(is.null(y) & Sys.info()[[1]]=="Linux" & is.null(lines)) {
    # figure out whether to use 'cat', 'zcat' or 'xzcat'
    suffix <- substr(file,nchar(file)-2,nchar(file))
    if(suffix==".gz") mycat <- "zcat"
    else if(suffix==".xz") mycat <- "xzcat"
    else mycat <- "cat"
    # use the system grep
    if(is.null(startswith)) startswith <- "" else startswith <- paste("^",startswith,".*",sep="")
    if(ignore.case) ic <- "-i" else ic <- ""
    out <- lapply(1:length(pattern), sysgrep)
  } else {
    # use R grep
    if(is.null(lines)) lines <- readLines(file)
    out <- lapply(1:length(pattern), Rgrep)
  }
  # make numeric (NA for ones that aren't matched)
  out <- as.numeric(sapply(out,as.numeric))
  return(as.numeric(out))
}

read.fasta <- function(file, i=NULL, ret="count", lines=NULL, ihead=NULL,
  start=NULL, stop=NULL, type="protein", id=NULL) {
  # read sequences from a fasta file
  # some of the following code was adapted from 
  # read.fasta in package seqinR
  # value of 'i' is what sequences to read 
  # value of 'ret' determines format of return value:
  #   count: amino acid composition (same columns as thermo$protein, can be used by add.protein)
  #        or nucleic acid base composition (A-C-G-T)
  #   seq: amino acid sequence
  #   fas: fasta entry
  # value of 'id' is used for 'protein' in output table,
  #   otherwise ID is parsed from FASTA header (can take a while)
  is.nix <- Sys.info()[[1]]=="Linux"
  if(is.nix & is.null(lines)) {
    msgout("read.fasta: reading ",basename(file),"\n")
    # figure out whether to use 'cat', 'zcat' or 'xzcat'
    suffix <- substr(file,nchar(file)-2,nchar(file))
    if(suffix==".gz") mycat <- "zcat"
    else if(suffix==".xz") mycat <- "xzcat"
    else mycat <- "cat"
    nlines <- as.integer(system(paste(mycat,' "',file,'" | wc -l',sep=""),intern=TRUE))
    ihead <- as.integer(system(paste(mycat,' "',file,'" | grep -n "^>" | cut -f 1 -d ":"',sep=""),intern=TRUE))
    #linefun <- function(i1,i2) as.character(system(paste('sed -n ',i1,',',i2,'p ',file,sep=""),intern=TRUE))
    lines <- system(paste(mycat,' "',file,'"',sep=""),intern=TRUE)
    linefun <- function(i1,i2) lines[i1:i2]
  } else {
    if(is.null(lines)) {
      lines <- readLines(file)
      msgout("read.fasta: reading ",basename(file),"\n")
    }
    nlines <- length(lines)
    if(is.null(ihead)) ihead <- which(substr(lines,1,1)==">")
    linefun <- function(i1,i2) lines[i1:i2]
  }
  # identify the lines that begin and end each sequence
  if(is.null(i)) {
    i <- ihead
    begin <- i + 1
    end <- i - 1
    end <- c(end[-1], nlines)
  } else {
    begin <- i + 1
    iend <- match(i,ihead)
    # we have to be careful about the last record
    iend[iend==ihead[length(ihead)]] <- NA
    end <- ihead[iend+1] - 1
    end[is.na(end)] <- nlines
  } 
  # just return the lines from the file
  if(ret=="fas") {
    iline <- numeric()
    for(i in 1:length(begin)) iline <- c(iline,(begin[i]-1):end[i])
    return(lines[iline])
  }
  # get each sequence from the begin to end lines
  seqfun <- function(i) paste(linefun(begin[i],end[i]),collapse="")
  sequences <- lapply(1:length(i), seqfun)
  # organism name is from file name
  # (basename minus extension)
  bnf <- strsplit(basename(file),split=".",fixed=TRUE)[[1]][1]
  organism <- bnf
  # protein/gene name is from header line for entry
  # (strip the ">" and go to the first space)
  if(is.null(id)) id <- as.character(palply("", 1:length(i), function(j) {
    # get the text of the line
    f1 <- linefun(i[j],i[j])
    # stop if the first character is not ">"
    # or the first two charaters are "> "
    if(substr(f1,1,1)!=">" | length(grep("^> ",f1)>0))
      stop(paste("file",basename(file),"line",j,"doesn't begin with FASTA header '>'."))
    # discard the leading '>'
    f2 <- substr(f1, 2, nchar(f1))
    # keep everything before the first space
    return(strsplit(f2," ")[[1]][1])
  } ))
  if(ret=="count") {
    counts <- count.aa(sequences, start, stop, type)
    ref <- abbrv <- NA
    chains <- 1
    if(type=="protein") {
      colnames(counts) <- aminoacids(3)
      # 20090507 made stringsAsFactors FALSE
      out <- cbind(data.frame(protein=id, organism=organism,
        ref=ref, abbrv=abbrv, chains=chains, stringsAsFactors=FALSE), counts)
    } else if(type %in% c("DNA", "RNA")) {
      out <- cbind(data.frame(gene=id, organism=organism,
        ref=ref, abbrv=abbrv, chains=chains, stringsAsFactors=FALSE), counts)
    }
  } else return(sequences)
}

uniprot.aa <- function(protein, start=NULL, stop=NULL) {
  # download protein sequence information from UniProt
  iprotein <- numeric()
  # construct the initial URL
  proteinURL <- paste("http://www.uniprot.org/uniprot/", protein, sep="")
  msgout("uniprot.aa: trying ", proteinURL, " ...")
  # try loading the URL, hiding any warnings
  oldopt <- options(warn=-1)
  URLstuff <- try(readLines(proteinURL),TRUE)
  options(oldopt)
  if(class(URLstuff)=="try-error") {
    msgout(" failed\n")
    return(NA)
  }
  # 20091102: look for a link to a fasta file
  linkline <- URLstuff[[grep("/uniprot/.*fasta", URLstuff)[1]]]
  # extract accession number from the link
  linkhead <- strsplit(linkline, ".fasta", fixed=TRUE)[[1]][1]
  accession.number <- tail(strsplit(linkhead, "/uniprot/", fixed=TRUE)[[1]], 1)
  msgout(" accession ", accession.number, " ...\n")
  # now download the fasta file
  fastaURL <- paste("http://www.uniprot.org/uniprot/", accession.number, ".fasta", sep="")
  URLstuff <- readLines(fastaURL)
  # show the name of the protein to the user
  header <- URLstuff[[1]]
  header2 <- strsplit(header, paste(protein, ""))[[1]][2]
  header3 <- strsplit(header2, " OS=")[[1]]
  protein.name <- header3[1]
  header4 <- strsplit(header3[2], " GN=")[[1]][1]
  header5 <- strsplit(header4[1], " PE=")[[1]]
  organism.name <- header5[1]
  msgout("uniprot.aa: ", protein.name, " from ", organism.name)
  # 20130206 use read.fasta with lines, start, stop arguments
  aa <- read.fasta(file="", lines=URLstuff, start=start, stop=stop)
  msgout(" (length ", sum(aa[1, 6:25]), ")\n", sep="")
  po <- strsplit(protein, "_")[[1]]
  aa$protein <- po[1]
  aa$organism <- po[2]
  return(aa)
}

count.aa <- function(seq, start=NULL, stop=NULL, type="protein") {
  # count amino acids or DNA bases in one or more sequences given as elements of the list seq
  # put them in alphabetical order (amino acids: same order as in thermo$protein)
  if(type=="protein") letts <- aminoacids(1)
  else if(type=="DNA") letts <- c("A", "C", "G", "T")
  else stop(paste("unknown sequence type", type))
  # to count the letters in each sequence
  countfun <- function(seq, start, stop) {
    # get a substring if one or both of start or stop are given
    # if only one of start or stop is given, get a default value for the other
    if(!is.null(start)) {
      if(is.null(stop)) stop <- nchar(seq)
      seq <- substr(seq, start, stop)
    } else if(!is.null(stop)) {
      seq <- substr(seq, 1, stop)
    }
    # the actual counting
    nnn <- table(strsplit(toupper(seq), "")[[1]])
    # get them in alphabetical order
    ilett <- match(names(nnn), letts)
    # in case any letters aren't in our alphabet
    ina <- is.na(ilett)
    if(any(ina)) msgout(paste("count.aa: unrecognized letter(s) in", type, "sequence:",
      paste(names(nnn)[ina], collapse=" "), "\n"))
    count <- numeric(length(letts))
    count[ilett[!ina]] <- nnn[!ina]
    return(count)
  }
  # counts for each sequence
  a <- palply("", seq, countfun, start, stop)
  a <- t(as.data.frame(a, optional=TRUE))
  # clean up row/column names
  colnames(a) <- letts
  rownames(a) <- 1:nrow(a)
  return(a)
}

