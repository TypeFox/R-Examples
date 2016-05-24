#$Id: fasta.R 189 2014-09-06 08:22:07Z larssn $


readFasta <- function( in.file ){
  # Lars Snipen and Kristian Hovde Liland
  # Biostatistics, Norwegian University of Life Sciences
  
  all.lines <- readLines( in.file )
  header.lines <- grep( "^>", all.lines )
  ns <- length( header.lines )
  if( ns == 0 ){
    stop( "no FASTA sequences found" )
  }
  dp <- header.lines + 1
  end <- c(header.lines[-1] - 1, length(all.lines) )
  fdta <- data.frame(
    Header   = gsub( "^>", "", all.lines[header.lines]),
    Sequence = vapply( (1:ns), function(i) {paste( all.lines[dp[i]:end[i]], collapse = "")}, character(1)), stringsAsFactors=F )
  class( fdta ) <- c( "Fasta", "data.frame" )
  return( fdta )
}


writeFasta <- function( fdta, out.file, width=80 ){
  # Lars Snipen and Kristian Hovde Liland
  # Biostatistics, Norwegian University of Life Sciences
  cn <- colnames( fdta )
  if( !("Header" %in% cn) | !("Sequence" %in% cn) ){
    stop( "This is not a fasta object, Header or Sequence is lacking\n" )
  }
  ns <- dim( fdta )[1]
  nc <- nchar( fdta$Sequence )
  h <- paste( ">", fdta$Header, sep="" )
  all.lines <- sapply( 1:ns, function( i ){
    s1 <- seq( from=1, to=nc[i], by=width )
    s2 <- union( seq( from=min( width, nc[i] ), to=nc[i], by=width ), nc[i] )
    z <- substring( fdta$Sequence[i], s1, s2 )
    return( c( h[i], z ) )
  })
  writeLines( unlist( all.lines ), con=out.file, sep="\n" )
}


plot.Fasta <- function( x, col="tan4", border="tan4", ... ){
  Fasta <- x
  nc <- nchar( Fasta$Sequence )
  ns <- length( nc )
  barplot( nc[ns:1], names.arg=ns:1, col=col, border=col, horiz=TRUE, xlab="Sequence length", ylab="Sequence number", main="Fasta object", ... )
}


str.Fasta <- function( object, ... ){
  Fasta <- object
  n.seq <- dim( Fasta )[1]
  cat( "Fasta:", n.seq, "sequences\n" )
}


summary.Fasta <- function( object, ... ){
  Fasta <- object
  n.seq <- dim( Fasta )[1]
  alphabet <- unique( unlist( strsplit( Fasta$Sequence, split="" ) ) )
  cat( "Fasta formatted sequence data containing", n.seq, "sequences\n" )
  cat( "Alphabet:", sort( alphabet ), "\n" )
}
