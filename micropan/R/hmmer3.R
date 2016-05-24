#$Id: hmmer3.R 189 2014-09-06 08:22:07Z larssn $


hmmerScan <- function( in.files, db, out.folder, verbose=TRUE ){
  basic <- "hmmscan -o hmmer3.txt --cut_ga --noali --cpu 0 --domtblout"
  db.name <- rev( unlist( strsplit( db, split="/" ) ) )[1]
  for( i in 1:length( in.files ) ){
    gi <- gregexpr( "GID[0-9]+", in.files[i], extract=T )
    rname <- paste( gi, "_vs_", db.name, ".txt", sep="" )
    res.files <- dir( out.folder )
    if( !(rname %in% res.files) ){
      if( verbose ) cat( "hmmerScan: Scanning", in.files[i], "...\n" )
      command <- paste( basic, file.path( out.folder, rname ), db, in.files[i]  )
      system( command )
    }
  }
}



readHmmer <- function( hmmer.file, e.value=1, use.acc=TRUE ){
  al <- readLines( hmmer.file )
  al <- al[which( !grepl( "\\#", al ) )]
  al <- gsub( "[ ]+", " ", al )
  lst <- strsplit( al, split=" " )
  if( use.acc ){
    hit <- sapply( lst, function(x){ x[2] } )
  } else {
    hit <- sapply( lst, function(x){ x[1] } )
  }
  query <- sapply( lst, function(x){ x[4] } )
  ievalue <- as.numeric( sapply( lst, function(x){ x[13] } ) )
  score <- as.numeric( sapply( lst, function(x){ x[14] } ) )
  start <- as.numeric( sapply( lst, function(x){ x[18] } ) )
  stopp <- as.numeric( sapply( lst, function(x){ x[19] } ) )
  desc <- sapply( lst, function(x){ paste( x[23:length( x )], collapse=" " ) } )
  hmmer.table <- data.frame( Query=query, Hit=hit, Evalue=ievalue, Score=score, Start=start, Stop=stopp, Description=desc, stringsAsFactors=F )
  hmmer.table <- hmmer.table[which( hmmer.table$Evalue <= e.value ),]
  return( hmmer.table )
}

