########################################################################
#                         oriloc
#
# Prediction of replication boundaries in unannotated genomes
#
########################################################################

oriloc <- function(
  seq.fasta = "ftp://pbil.univ-lyon1.fr/pub/seqinr/data/ct.fasta",
  g2.coord = "ftp://pbil.univ-lyon1.fr/pub/seqinr/data/ct.predict",
  glimmer.version = 3,
  oldoriloc = FALSE,
  gbk = NULL,
  clean.tmp.files = TRUE,
  rot = 0)
{
  aGBKfileWasGiven <- !missing(gbk) && !is.null(gbk)
  if(aGBKfileWasGiven) # Work directly with genbank file
  {
    tmpgbk <- tempfile(pattern = "orilocgbk")
    if(substr(gbk,1,7)=="http://" || substr(gbk,1,6)=="ftp://" || substr(gbk,1,7)=="file://"){
      download.file( gbk, destfile = tmpgbk )
    }
    else{
      file.copy(from = gbk, to = tmpgbk)
    }
    seq.fasta <- tempfile(pattern = "orilocfasta")
    g2.coord <- tempfile(pattern = "orilocg2")
   
    gb2fasta( tmpgbk, seq.fasta )
    gbk2g2( tmpgbk, g2.coord )
    #
    # gbk2g2 yields glimmer version 2.0 files, so force to version 2.0 in this case:
    #
    glimmer.version <- 2
  } 
#  
# Get first sequence from fasta file:
#
  seq <- read.fasta(file = seq.fasta, set.attributes = FALSE)[[1]]
  lseq <- length(seq)
#
# Read CDS coordinate file:
#
  g2 <- readLines( g2.coord )
#
# Patch for glimmer3 version:
#
  if( glimmer.version > 2 ){
    # remove first line:
    g2 <- g2[-1]
    # remove first three characters (i.e. orf)
    g2 <- sapply(g2, function(x) substr(x,4,nchar(x)), USE.NAMES = FALSE)
  }
#
# Extract info from g2.coord file
#
  tokens <- function( string )
  {
    tmp <- unlist(strsplit( string, split = " "))
    tmp[nchar(tmp) > 0 ][seq_len(3)]
  }
  tmp <- sapply( g2, tokens )
  gnum  <- as.numeric(tmp[1, ]) # gene number in g2.coord
  start <- as.numeric(tmp[2, ]) # start positions in bp
  end   <- as.numeric(tmp[3, ]) # end position in bp

  if( length(start) != length(end) )
    stop("start and end vector must be the same length")
#
# Rotate the genome if required
#
  if( rot != 0 )
  {
    #
    # Circular permutation of a vector
    #
    rotate <- function(x, rot = 0) 
    {
      n <- length(x)
      rot <- rot %% n
      if(rot == 0)
      {
        x
      } 
      else 
      {
        c(x[(rot+1):n],x[seq_len(rot)])
      }
    } 
    seq <- rotate(x  = seq, rot = rot )
    start <- start - rot
    end <- end - rot
    
    start[ start < 1 ] <- start[ start < 1 ] + lseq
    end[ end < 1] <- end[ end < 1 ] + lseq
    
    gnum <- gnum[order(start)]
    end <- end[order(start)]

    start <- sort(start)
  }
#
#
#
  pos <- (start + end)/2000    # Mid gene position in Kb
  ncds <- length(pos)
#
# CDS that wrap around the genome
#
  wrap <- abs(end-start) > lseq/2

#
# Compute the DNA walk gene by gene in third codon positions
#
  x <- integer(ncds)
  y <- integer(ncds)
  skew <- numeric(ncds)
  CDS.excess <- integer(ncds)

  for( i in order(pos) )
  {
    # Look for third codon position

    if( !wrap[i] ) # regular cds that do not wrap around the genome
    {
      if( start[i] < end[i] ) # CDS 5'->3' direct strand
      {
        CDS.excess[i] <- 1
        tcp <- seq( from = start[i] + 2, to = end[i], by = 3)
      }
      else # complementary strand
      {
        CDS.excess[i] <- -1
        tcp <- seq( from = start[i] - 2, to = end[i], by = -3)
      }
    }
    else # a cds that wraps around the genome
    {
      if( start[i] > end[i] ) # CDS 5'->3' direct strand
      {
        CDS.excess[i] <- 1
        tcp <- seq( from = start[i] + 2, to = lseq + end[i], by = 3)
        tcp[ tcp > lseq ] <- tcp[ tcp > lseq ] - lseq
      }
      else # CDS 3'->5' complementary strand
      {
        CDS.excess[i] <- -1
        tcp <- seq( from = start[i] - 2, to = -(lseq-end[i]), by = -3)
        tcp[ tcp < 1 ] <- tcp[ tcp < 1 ] + lseq
      }
    }

    tcnucl <- seq[tcp]
    x[i] <- length(tcnucl[tcnucl=="t"]) - length(tcnucl[tcnucl=="a"])
    y[i] <- length(tcnucl[tcnucl=="c"]) - length(tcnucl[tcnucl=="g"])
  }
  x <- cumsum(x)
  y <- cumsum(y)
  CDS.excess <- cumsum(CDS.excess)
#
# Old oriloc program, direct from C without trying to vectorize.
# To reproduce old results.
#
  if( oldoriloc )
  {
    Regression <- function(x, y, Li)
    {
      a <- 0 ; b <- 0 ; c <- 0;
      for( m in seq_len(Li-1) ) # I think this should go to Li included
      {
        b <- b + y[m]^2
        a <- a + x[m]^2
        c <- c + x[m]*y[m]
      }
      alfa1 <- (atan(2*c/(a-b)))/2 
      alfa2 <- alfa1 - pi/2;
      derivate1 <- 2*(a-b)*cos(2*alfa1)+4*c*sin(2*alfa1)

      if(derivate1 < 0)
        return( tan(alfa2) )
      else
        return( tan(alfa1) )
    }

    slope <- Regression( x, y, ncds )
  
    for ( i in seq_len(ncds))
    {
      X.line <- ( y[i] + slope*x[i] )/(2*slope)
      Y.line <- slope*X.line        
      distance <- sqrt( Y.line^2 + X.line^2 )

      if( Y.line < 0 )
        distance <- -distance

      skew[i] <- distance
    }
  }
  else # New oriloc program
  {
#
# Project DNAwalk points (x,y) onto orthogonal regression line
#
    	pca <- ade4::dudi.pca( cbind(x,y), scann = FALSE, nf = 1, scale = FALSE, center = FALSE )
    	rec <- ade4::reconst(pca)
    	skew <- sign(rec$x)*sqrt(rec$x^2+rec$y^2)
  }
#
# Try to get get a correct orientation (same as GC skew)
#
  if( cor(skew, y ) < 0 )
    skew <- -skew
#
# Build result
#
  result <- data.frame( cbind( gnum, start/1000, end/1000,
    CDS.excess, skew, x, y) )
  names(result) <- c("g2num", "start.kb", "end.kb", "CDS.excess",
    "skew","x","y")
#
# Delete temporary files if requested:
#
  if(aGBKfileWasGiven && clean.tmp.files)
  {
    file.remove(tmpgbk)
    file.remove(seq.fasta)
    file.remove(g2.coord)
  } 
#
# the end
#
  return(result)
}
