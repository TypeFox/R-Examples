#$Id: entrez.R 125 2013-06-28 09:03:45Z larssn $

entrezDownload <- function( accession, out.file, verbose=TRUE ){
  if( verbose ) cat( "Downloading genome..." )
  connect <- file( out.file, open="w" )
  for( j in 1:length( accession ) ){
    adr <- paste( "http://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=nucleotide&id=", accession[j], "&retmode=text&rettype=fasta", sep="" )
    entrez <- url( adr, open="rt" )
    if( isOpen( entrez ) ){
      lines <- readLines( entrez )
      writeLines( lines, con=connect )
      close( entrez )
    } else {
      cat( "Download failed: Could not open connection\n" )
    }
  }
  close( connect )
  if( verbose ) cat( "...sequences saved in", out.file, "\n" )
  return( out.file )
}


getAccessions <- function( master.record.accession ){
  adr <- paste( "http://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=nucleotide&id=", master.record.accession, "&retmode=text&rettype=gb", sep="" )
  entrez <- url( adr, open="rt" )
  accessions <- ""
  if( isOpen( entrez ) ){
    lines <- readLines( entrez )
    close( entrez )
    wgs.line <- gsub( "WGS[ ]+", "", lines[grep( "WGS   ", lines )] )
    ss <- unlist( strsplit( wgs.line, split="-" ) )
    head <- gregexpr( "[A-Z]+[0]+", ss[1], extract=T )
    ss.num <- as.numeric( gsub( "^[A-Z]+[0]+", "", ss ) )
    if( length( ss.num ) > 1 ){
      range <- ss.num[1]:ss.num[2]
    } else {
      range <- ss.num
    }
    ns <- ceiling( length( range ) / 500 )
    accessions <- character( ns )
    for( j in 1:ns ){
      s1 <- (j-1)*500 + 1
      s2 <- min( j*500, length( range ) )
      accessions[j] <- paste( paste( head, range[s1]:range[s2], sep="" ), collapse="," )
    }
  } else {
    cat( "Download failed: Could not open connection\n" )
  }
  return( accessions )
}
