readGenBank <-
function( gbkfile ){
  
  all.lines <- readLines( con=gbkfile )
  gen.idx <- grep(" CDS        | tRNA        | rRNA        ", all.lines )
  gen.lines <- gsub( " +", ";", gsub( "^ +", "", all.lines[gen.idx] ) )
  join.idx <- grep( "join\\(", gen.lines )
  if( length( join.idx ) > 0 ){
    cat("Your genbank file contains complex pseudogene structures (keyword: join). CNOGpro does not currently guarantee correct output for these structures.")
    xtra.lines <- NULL
    xtra.idx <- NULL
    for( i in 1:length( join.idx ) ){
      gl <- gen.lines[join.idx[i]]
      tox <- unlist( strsplit( gl, split=";" ) )
      typ <- tox[1]
      jtxt <- tox[2]
      cc <- 1
      while( length( grep( "\\)$", jtxt ) ) == 0 ){
        jtxt <- paste( jtxt, gsub( " +", "", all.lines[gen.idx[join.idx[i]]+cc] ), sep="" )
        cc <- cc + 1
      }
      jtxt <- gsub( "join\\(|\\)$", "", jtxt )
      exon <- unlist( strsplit( jtxt, split="," ) )
      gen.lines[join.idx[i]] <- paste( typ, exon[1], sep=";" )
      for( j in 2:length( exon ) ){
        xtra.lines <- c( xtra.lines, paste( typ, exon[j], sep=";" ) )
        xtra.idx <- c( xtra.idx, gen.idx[join.idx[i]] )
      }
    }
    gen.lines <- c( gen.lines, xtra.lines )
    gen.idx <- c( gen.idx, xtra.idx )
    ixx <- order( gen.idx )
    gen.idx <- gen.idx[ixx]
    gen.lines <- gen.lines[ixx]
  }
  
  loc.idx <- grep( "/locus_tag=", all.lines )
  loc.lines <- gsub( "^ +|/locus_tag=|\"", "", all.lines[loc.idx] )
  idd <- sapply( gen.idx, function(x){min( which( (loc.idx-x)>0 ) )} )
  loc.tag <- loc.lines[idd]
  
  strand <- rep( 1, length( gen.lines ) )
  strand[grep( "complement", gen.lines )] <- -1
  gen.lines <- gsub( "\\.\\.", ";", gsub( "complement\\(|\\)", "", gen.lines ) )
  tmat <- sapply( strsplit( gen.lines, split=";" ), function(x){c(x[1],x[2],x[3])} )
  left <- as.numeric( gsub( "[\\<\\>]", "", tmat[2,] ) )
  right <- as.numeric( gsub( "[\\<\\>]", "", tmat[3,] ) )
  dfr <- data.frame( Type=tmat[1,], Locus=loc.tag, Strand=strand, Left=left, Right=right, Length=right-left+1,
                     stringsAsFactors=F )
  gbkfilename <- basename(gbkfile)
  gb2fasta(source.file=gbkfile, destination.file=paste(tempdir(),"/",gsub(".gbk|.gb|.gff|.genbank|.txt|.text","",gbkfilename),".fasta",sep=""))
  
  # Now getting chromosome name and length
  header <- (gsub(x=all.lines[grep("LOCUS",x=all.lines)],pattern="\\s+",replacement=","))
  header <- sapply( strsplit(header,split=","), function(x){c(x[2],x[3])})

  #dfr <- rbind(c("CHR",header[1,],1,1,header[2,],header[2,]),dfr)
  dfr <- insertRow(dataframe=dfr,newrow=c("CHR",header[1,],1,1,header[2,],header[2,]),1)
  dfr[,3:6] <- sapply(dfr[,3:6],function(x) (as.numeric(x)))
  return( dfr )
}
