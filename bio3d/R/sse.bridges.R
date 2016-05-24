"sse.bridges" <- function(sse, type="helix", hbond=TRUE, 
                          energy.cut=-1.0 ) {
  if(missing(sse))
    stop("sse missing")

  if(is.null(sse$hbonds))
    stop("sse$hbonds does not exists. run dssp with 'resno=FALSE' and 'full=TRUE'")
  
  natoms <- nrow(sse$hbonds)
  
  if(type=="helix") {
    sse2 <- sse$helix
    stype <- "H"
    lim <- 4
  }
  
  if(type=="sheet") {
    sse2 <- sse$sheet
    stype <- "S"
    lim <- 2
  }

  if(length(sse2$start)==0)
    return(NULL)

  ## character array of SSE membership
  simple.sse <- sse$sse
  simple.sse[ !(sse$sse %in% c("H", "E", "G", "I")) ] <- "L"
  simple.sse[ sse$sse %in% c("E") ]                   <- "S"
  simple.sse[ sse$sse %in% c("H", "I", "G") ]         <- "H"
  
  inds <- NULL
  for ( i in 1:(natoms-lim) ) {
    if(simple.sse[i]!=stype)
      next;

    paired <- NULL
    if(type=="helix") {
      paired <- i+4
      if(simple.sse[paired]!=stype)
        next;
    }

    if (type=="sheet") {
      ## bridge pair info
      paired <- sse$hbonds[i,c("BP1", "BP2")]
      paired <- paired[ !is.na(paired) ]

      if(length(paired)==0)
        next;
    }

    if(hbond) {
      ## H-bond info
      resid <- sse$hbonds[i, c(3,5,7,9)]
      energ <- sse$hbonds[i, c(4,6,8,10)]
      
      energ <- energ[ !is.na(resid) ]
      resid <- resid[ !is.na(resid) ]
      
      inds.tmp <- which(resid %in% paired)
      resid <- resid[inds.tmp]
      energ <- energ[inds.tmp]
      
      dups <- duplicated(resid)
      resid <- resid[ !dups ]
      energ <- energ[ !dups ]

      if(length(resid)==0)
        next;
      
      for( j in 1:length(resid) ) {
        if(energ[j] < energy.cut)
          inds <- c(inds, i, resid[j])
      }
    }

    else {
      for ( j in 1:length(paired) ) {
        inds <- c(inds, i, paired[j])
      }
    }
  }

  if(length(inds)==0)
    return(NULL)
  
  mat      <- matrix(inds, ncol=2, byrow=T)
  mat      <- t(apply(mat, 1, sort))
  pair.ids <- apply(mat, 1, function(x) paste(x, collapse="-"))
  mat      <- matrix(mat[!duplicated(pair.ids), ], ncol=2)
  
  return(mat)
}
