print.fasta <- function(x, alignment=TRUE, ...) {
  
  if ( (inherits(x, "fasta") | inherits(x, "pdbs")) ){ 
    id <- x$id
    ali <- as.matrix(x$ali)
    call <- paste(deparse(x$call), sep = "\n", collapse = "\n")
  } else {
    ali <- as.matrix(x)
    id <- rownames(ali)
    call <- NA
  }

  if(inherits(x, "pdbs"))
    row.desc <- "structures"
  else
    row.desc <- "sequences"
  
  cn      <- class(x)
  nstruct <- length(id)
  dims    <- dim(ali)
  
  gaps <- gap.inspect(ali)
  dims.nongap <- dim(ali[, gaps$f.inds, drop=FALSE])
  dims.gap <- dim(ali[, gaps$t.inds, drop=FALSE])

  if(alignment)
    .print.fasta.ali(ali, id=id, ...)
                     #width=width, 
                     #col.inds=col.inds, numbers=numbers)
  else
    cat("\n")

  if(!is.na(call)) {
    cat("Call:\n  ", call,
        "\n\n", sep = "")
    
    cat("Class:\n  ", paste(cn, collapse=", "),
        "\n\n", sep = "")
  }
  
  ##cat("Number of ", row.desc, ":\n  ", nstruct,
  ##    "\n\n", sep="")
  
  cat("Alignment dimensions:\n  ",
      dims[1L], " sequence rows; ",
      dims[2L], " position columns",
      " (", dims.nongap[2L], " non-gap, ", dims.gap[2L], " gap) ", 
      "\n", sep="")
  
  cat("\n")

 
  ## Attribute summary
  j <- paste(attributes(x)$names, collapse = ", ")
  cat(strwrap(paste(" + attr:",j,"\n"),width=60, exdent=8), sep="\n")

}

.print.fasta.ali <- function(x, ##limit.out=NULL,
                             width=NULL, col.inds=NULL,
                             numbers=TRUE, conservation=TRUE, ...) {
  ##-- Print sequence alignment in a nice formated way
  ##    source("print.aln.R")
  ##    x<-read.fasta("poo.fa")
  ##    print.aln(x)
  ##
  ##    file <- system.file("examples/kif1a.fa",package="bio3d")
  ##    aln  <- read.fasta(file)
  ##    print.aln(aln, width=40)
  ##    print.aln(aln, width=60)
  ##    print.aln(aln, col=c(10,12,14,22,90:100), numbers=F)
  ## ToDo:
  ##   Does not work if alignment contains only one position (one seq?)
  ##     y=x; y$ali=x$ali[,1]

  if ( (inherits(x, "fasta") | inherits(x, "pdbs")) ){ 
    id <- x$id
    ali <- as.matrix(x$ali)
  } else {
    ali <- as.matrix(x)
    id <- rownames(x)
  }

  ## remove any NA values
  ali[is.na(ali)] <- "-"
  
  ##- Trim to 'col.inds' if provided 
  if(!is.null(col.inds)) {
    ali <- ali[,col.inds, drop=FALSE]
  }

  if(nrow(ali)<2)
    conservation <- FALSE
  
  ##- conservation
  cons <- NULL
  if(conservation) {
    tmp1 <- conserv(ali, method="entropy10")
    tmp2 <- conserv(ali, method="identity")
    cons <- rep(" ", ncol(ali))
    cons[ tmp1==1 ] <- "^"
    cons[ tmp2==1 ] <- "*"
  }
  
  ## Check and truncate possilbe long ids
  if(any(nchar(id) > 20)) {
    id <- basename(id)
    if(any(nchar(id) > 17)) {
        id <- substr(id,1,10)
    }
    id <-paste0("[Truncated_Name:", 1:length(id),"]",id)
  }

  ##- Format sequence identifiers
  ids.nchar <- max(nchar(id))+3 ## with a gap of 3 spaces btwn id and sequence
  ids.format <- paste0("%-",ids.nchar,"s")
  ids <- sprintf(ids.format, id)
  
  ## Format for annotation printing (see below)
  pad.format <- paste0("%+",(ids.nchar+1),"s")

  ## Format for conservation annotation printing (see below)
  pad.format2 <- paste0("%+",(ids.nchar),"s")
  
  ##- Scale 'width' of output if not specified in input call
  tput.col <- 85  ## typical terminal width from system("tput cols")
  if(is.null(width)) {
    width <- tput.col - ids.nchar - 4
  }
  
  ## Make sure we end on a 10 block
  width <- floor(width/10)*10
  
  ##- Work out sequence block widths
  nseq <- length(ids)
  nres <- ncol(ali)
  
  block.start <- seq(1, nres, by=width)
  if(nres < width) {
    block.end <- nres
  } else {
    block.end <- unique(c( seq(width, nres, by=width), nres))
  }
  nblocks <- length(block.start)
  
  block.annot  <- rep(" ", width)
  block.annot[ c(1,seq(10, width, by=10)) ] = "."
    
  blocks <- matrix(NA, ncol=nblocks, nrow=nseq) 
  for(i in 1:nblocks) {
    ##- Sequence block
    positions <- block.start[i]:block.end[i]
    blocks[,i] <- paste0(ids, apply(ali[, positions, drop=FALSE], 1, paste, collapse=""))
    
    ##-- Formated Printing of annotations (numbers & ticks) and sequence blocks
    if(numbers) {
      ##- Annotations for each sequence block
      annot = block.annot[1:length(positions)]
      annot[length(annot)] = block.end[i]
      annot[1] = sprintf(pad.format, block.start[i])
      cat(paste(annot, collapse=""),"\n")
    }
    
    ##- Sequence block
    cat(blocks[,i], sep="\n")

    ##- Formated Printing of conservation (stars for conserved columns)
    if(conservation) {
      annot2 <- c("", cons[positions])
      annot2[1] = sprintf(pad.format2, "")
      cat(paste(annot2, collapse=""),"\n")
    }
    
    ##- Ticks + numbers again
    if(numbers) {
      cat(paste(annot, collapse=""),"\n\n")
    } else{ cat("\n") }
  }
  
  ##invisible(blocks) ## Can be useful for plot.fasta() later!!
}

