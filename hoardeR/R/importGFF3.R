importGFF3 <- function(gff){
  # Determine the header lines:
  con  <- file(gff, open = "r")
  gff3Header <- c()
  header <- TRUE
  headerLines <- 0
  while (header) {
    oneLine <- readLines(con, n = 1, warn = FALSE)
    firstChr <- substr(oneLine,1,1)
    if(firstChr=="#"){
      headerLines <- headerLines + 1
      gff3Header[headerLines] <- oneLine
    } else {
      header <- FALSE
    }
  }
  close(con)
  
  # Load the file and split the Variable V9
  gff3Loaded <- read.csv(file=gff, sep="\t", header=FALSE, stringsAsFactors=FALSE, skip=headerLines)
  #cuffLoaded <- cuffLoaded[1:100,]
  V9 <- gff3Loaded$V9
  V9 <- strsplit(V9,";")
  
  # Now get the part before the space (=first) and after (=second)
  first <- sapply(V9, function(x)unlist(regmatches( x , gregexpr( "[_a-zA-Z]+=" , x ) )))
  second <- sapply(V9, function(x)unlist(regmatches( x , gregexpr( "=[_a-zA-Z0-9=.]+" , x ) )))
  first <- sapply(first,str_replace_all, "=", "")
  second <- sapply(second,str_replace_all, "=", "")
  
  #  Column names
  tags <- unique( unlist(first) )
  
  #  Intermeidate matrices
  temp <- mapply( cbind , second , first )
  
  #  Match to appropriate columns and coerce to data.frame
  out <- data.frame( do.call( rbind , lapply( temp , function(x) x[ match( tags , x[,2] ) ]  ) ) )
  names(out) <- tags
  out <- cbind(gff3Loaded[,-9],out)
  out
}