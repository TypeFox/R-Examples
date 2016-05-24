#$Id: blasting.R 170 2014-07-23 12:57:06Z larssn $

blastAllAll <- function( in.files, out.folder, e.value=1, job=1, verbose=T ){
  for( i in 1:length( in.files ) ){
    command <- paste( "makeblastdb -logfile log.txt -dbtype prot -out blastDB", job, " -in ", in.files[i], sep="" )
    system( command )
    file.remove( "log.txt" )
    gi <- gregexpr( "GID[0-9]+", in.files[i], extract=T )
    for( j in 1:length( in.files ) ){
      gj <- gregexpr( "GID[0-9]+", in.files[j], extract=T )    
      rname <- paste( gj, "_vs_", gi, ".txt", sep="" )
      res.files <- dir( out.folder )
      if( !(rname %in% res.files) ){
        if( verbose ) cat( "blastAllAll: ", rname, "\n" )
        input <- paste( "-query ", in.files[j], sep="" )
        dbase <- paste( "-db blastDB", job, sep="" )
        output <- paste( "-out ", file.path( out.folder, rname ), sep="" )
        command <- paste( "blastp -matrix BLOSUM45 -evalue", e.value, "-outfmt 6", input, dbase, output )
        system( command )
      }
    }
  }
  file.remove( paste( "blastDB", job, ".pin", sep="" ) )
  file.remove( paste( "blastDB", job, ".phr", sep="" ) )
  file.remove( paste( "blastDB", job, ".psq", sep="" ) )
  return( "done" )
}


readBlastTable <- function( blast.file ){
  columnames <- c( "Query", "Hit", "Percent.identity", "Alignment.length", "Mismatches", "Gap.openings", "Query.start", "Query.end", "Hit.start", "Hit.end", "E.value", "Bit.score" )
  data <- read.table( blast.file, sep="\t", quote="", header=F, col.names=columnames, stringsAsFactors=F )
  return( data )
}






