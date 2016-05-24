#$Id: prodigal.R 78 2013-05-12 10:05:43Z larssn $

prodigalPredict <- function( genome.file, prot.file, nuc.file=NULL, closed.ends=TRUE, motif.scan=FALSE ){
  command <- paste( "prodigal -i ", genome.file, " -a ", prot.file, " -o prodigal.txt -q", sep="" )
  if( !is.null( nuc.file ) ){
    command <- paste( command, " -d ", nuc.file, sep="" )
  }
  if( closed.ends ){
    command <- paste( command, " -c ", sep="" )
  }
  if( motif.scan ){
    command <- paste( command, " -n ", sep="" )
  }
  system( command )
  file.remove( "prodigal.txt" )
  return( paste( "Prodigal predictions in", prot.file ) )
}
