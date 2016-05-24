"seqidentity" <-
function( alignment , normalize=TRUE, similarity=FALSE, ncore=1, nseg.scale=1) {

  # Parallelized by parallel package (Sun Jul  7 17:35:38 EDT 2013)
  ncore <- setup.ncore(ncore)
  if(ncore > 1) {
     # Issue of serialization problem
     # Maximal number of cells of a double-precision matrix
     # that each core can serialize: (2^31-1-61)/8
     R_NCELL_LIMIT_CORE = 2.68435448e8
     R_NCELL_LIMIT = ncore * R_NCELL_LIMIT_CORE

     if(nseg.scale < 1) {
        warning("nseg.scale should be 1 or a larger integer\n")
        nseg.scale=1
     }
  }

  ids <- NULL
  if(is.list(alignment)) {
    if(inherits(alignment, c("fasta", "pdbs")))
      ids <- alignment$id
    alignment <- alignment$ali
  }
  else {
    ids <- rownames(alignment)
  }

  ## calculate similarity instead of identity?
  if(similarity) {
    alnTo10 <- function(x) {
      new <- rep(NA, length(x))
      aa <- c("V","I","L","M",  "F","W","Y",  "S","T",
              "N","Q",  "H","K","R",  "D","E",
              "A","G",  "P",  "C",  "-","X")
      
      new[ x %in% aa[1:4] ]   = "V"   # Hydrophobic, Aliphatic
      new[ x %in% aa[5:7] ]   = "F"   # Aromatic
      new[ x %in% aa[8:9] ]   = "S"   # Ser/Thr
      new[ x %in% aa[10:11] ] = "N"   # Polar
      new[ x %in% aa[12:14] ] = "R"   # Positive
      new[ x %in% aa[15:16] ] = "D"   # Negative
      new[ x %in% aa[17:18] ] = "A"   # Tiny
      new[ x %in% aa[19] ]    = "P"   # Proline
      new[ x %in% aa[20] ]    = "C"   # Cysteine
      new[ x %in% aa[21:22] ] = "-"   # Gaps
      
      return(matrix(new, nrow=1))
    }
    alignment <- t(apply(alignment, 1, alnTo10) )
  }
  alignment[is.gap(alignment)] = NA

  ide <- function(x, y) {
#### Edit by Heiko Strathmann
#### Wed Aug  4 10:48:16 PDT 2010
#### Fix for bug with all gap sequences
    r <- sum(x==y, na.rm=TRUE)
    t <- sum(complete.cases(cbind(x,y)))
    if (normalize && t != 0) {
      r <- r/t
    }
##################################
    return( round(r, 3) )
  }

  nseq <- nrow(alignment)
  inds <- pairwise( nseq )
  ni <- nrow(inds)

  if(ncore > 1) {
     RLIMIT = R_NCELL_LIMIT
     nDataSeg = floor((ni-1)/RLIMIT) + 1
     nDataSeg = floor(nDataSeg * nseg.scale)
     lenSeg = floor(ni/nDataSeg)
     s = NULL
     for(i in 1:nDataSeg) {
        istart = (i-1)*lenSeg + 1
        iend = if(i<nDataSeg) i*lenSeg else ni
        s <- c(s, mclapply(istart:iend, function(j) {
           ide(alignment[inds[j,1],], alignment[inds[j,2],])
          }) )
     }
     s <- unlist(s)
  } else {
     s <- rep(NA, ni)
   
     for(i in 1:ni) {
       s[i]<-ide(alignment[inds[i,1],], alignment[inds[i,2],])
     }
  }
  ## make 's' into matrix 'sm'
  sm <- matrix(1, ncol=nseq,nrow=nseq)
  sm[inds]<-s
  if(nseq==2) {
    sm[inds[,2], inds[,1]]<-s
  } else {
    sm[inds[,c(2,1)]]<-s 
  }

  if(!is.null(ids)) {
    rownames(sm) <- ids
    colnames(sm) <- ids
  }
  return(sm) # ide matrix
}

