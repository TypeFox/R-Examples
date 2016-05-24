#$Id: panprep.R 170 2014-07-23 12:57:06Z larssn $

panPrep <- function( in.file, GID.tag, out.file, protein=TRUE, discard=NA ){
  if( protein ){
    alien <- "[^ARNDCQEGHILKMFPSTWYV]"
    minl <- 10
  } else {
    alien <- "[^ACGT]"
    minl <- 30
  }
  fdta <- readFasta( in.file )
  if( !is.na( discard ) ) fdta <- fdta[grep( discard, fdta$Header, invert=T ),]
  fdta$Sequence <- toupper( gsub( "\\*", "", fdta$Sequence ) )
  fdta <- fdta[which( nchar( fdta$Sequence ) >= minl ),]
  fdta$Sequence <- gsub( alien, "X", fdta$Sequence )
  nseq <- dim( fdta )[1]
  tok1 <- paste( GID.tag, paste( "seq", 1:nseq, sep="" ), sep="_" )
  fdta$Header <- paste( tok1, fdta$Header )
  fext <- unlist( gregexpr( "\\.[a-zA-Z]+$", out.file, extract=T ) )
  fname <- paste( gsub( "\\.[a-zA-Z]+$", "", out.file ), "_", GID.tag, fext, sep="" )
  writeFasta( fdta, out.file=fname )
}