stripclosestpairs <- function(rownamedarray, finallength, scales, fixedNbr=0) { # return rownames!
  krigmax <- blackbox.getOption("krigmax")
  ## the fixed points should be the last fixedNbr rows in rownamedarray !!
  #finallength includes fixed points
  if (nrow(rownamedarray)<=finallength) {
    if (nrow(rownamedarray)<finallength)
      cat("(!) In stripclosestpairs(...): size of input array lower than target final size", "\n")
    return(rownames(rownamedarray))
  } 
  localarray <- sweep(rownamedarray,2L,sqrt(scales),FUN=`/`) # t(t(rownamedarray)/sqrt(scales))
  nnr <- nrow(localarray) ## includes fixed points
  nnc <- ncol(localarray)
  nnf <- nnr-fixedNbr ## only non-fixed points
  if (fixedNbr>0) {
    fixed <- (nnf+1):nnr
    fixedpts <- matrix(localarray[fixed, ], ncol=nnc)
  } else {
    fixed <- c()
  }
  bandwidth <- as.integer((krigmax^2)/nnr)
  if (bandwidth<krigmax) { ## then consider partial distance matrix in banded form
    distvec <- rep(NA, nnf*(bandwidth+fixedNbr))
    pos <- 0 ###
    for (ll in 1:nnf) { ## constructs distances
      tmp <- localarray[ll, ]
      ## distances between nonfixed points
      if (ll<nnf) {
        bandpts <- matrix(localarray[(ll+1):min(ll+bandwidth, nnf), ], ncol=nnc)
        pif <- apply(bandpts, 1, function(v) {sum((v-tmp)^2)})
        ###distvec <- append(distvec, pif) ## distances between nonfixed points
        dl <- length(pif)
        distvec[(pos+1):(pos+dl)] <- pif ###
        pos <- pos+dl
      }
      ## next, case where indexing up to band width would imply distances with fixed points or even with nothing
      pos <- ll*(bandwidth+fixedNbr)-fixedNbr # skips positions which thus stay filled with NA's
      ## next true distances between nonfixed and fixed points
      pif <- apply(fixedpts, 1, function(v) {sum((v-tmp)^2)})
      dl <- length(pif)
      distvec[(pos+1):(pos+dl)] <- pif
      pos <- pos+dl
    } ## final distvec is nnf*(bandwidth+fixedNbr) array stored by rows in a vector, no within-fixed comparisons
    maxv <- max(distvec, na.rm=T)
    distvec[is.na(distvec)] <- maxv ## those involving dummy columns and those 'there'
    distmat <- matrix(distvec, nrow=nnf)
    distvec <- NA
    rownames(distmat) <- rownames(rownamedarray)[1:nnf]
    nr <- nrow(distmat); ## maxi nnf hence not the fixed coordinates
    ## 03/2014: all 'zut' related code is a way to avoid recursively calling order()
    ## it has been tested only in the 'full matrix' case... and it appears slower than recusively cally order()
    # zutindices <- distmat
    # zutindices[TRUE, TRUE] <- seq_len(length(distmat))
    # zutordre <- order(as.vector(distmat)) ## will be updated to give the order of remainig elements with indices referring to positions in the reduced distmat
    while (nr>(finallength-fixedNbr)) {
      ordre <- order(as.vector(distmat))[1:(nr-finallength+fixedNbr)] ## what' being done, but without calling order() every time;
      # ordre <- zutordre[1:(nr-finallength+fixedNbr)]
      removand <- c()
      while (length(ordre)>0) {
        rows <- ordre %% nr # numbers, not rownames
        rows[rows==0] <- nr ## each element of ordre is a position in as.vector(distmat), we recover the matrix positions (i, j)
        cols <- (ordre-1) %/% nr
        cols <- cols+1
        trows <- table(rows) # its names are row indices, not rownames
        choix <- as.numeric(names(trows)[trows==max(trows)])
        if (length(choix)>1) {
          rowmins <- apply(distmat[choix, ], 1, min) #names of this are row indices of distmat
          minrowmins <- min(rowmins)
          namedrowindices <- which(rowmins==minrowmins) ## may be a pair of nonfixed points
          if (length(namedrowindices)>1) { ## further  chose within this pair based on other values in their rows
            tmp <- distmat[namedrowindices, ]
            atuer <- c()
            for (cc in 1:nrow(tmp)) atuer <- append(atuer, min(tmp[cc, tmp[cc, ]!=minrowmins]))
            rownam <- names(namedrowindices)[which.min(atuer)] ## has kept original row name of distmat
            choix <- which(rownames(distmat)==rownam)
          } else choix <- which(rownames(distmat)==names(namedrowindices))
        }
        ordre <- ordre[rows!=choix]
        cols <- cols[rows!=choix]
        ordre <- ordre[cols!=choix]
        removand <- append(removand, choix)
      }
      for (ll in removand) { ## effectively discards 'columns' from further consideration
        if (ll>1) for (lll in 1:min(bandwidth, ll-1)) distmat[ll-lll, lll] <- maxv
      }
      distmat <- distmat[-removand, ]
      nr <- nrow(distmat)
      # dif <- setdiff(zutordre, unique(c(as.vector(zutindices[removand, ]), as.vector(zutindices[, removand]))))
      # zutindices[, removand] <- NA
      # zutindices[removand, ] <- NA
      # zutindices[!is.na(zutindices)] <- seq_len(length(distmat)) ##
      # zutordre <- zutindices[dif] ## order with new indices from order with old indices
      # zutindices <- zutindices[-removand, -removand]
    }
    return(c(rownames(distmat), rownames(rownamedarray)[fixed]))
  } else { ## full matrix
    distvec <- proxy::dist(localarray) ##  dist object lower triangle of the distance matrix stored by columns in a vector but printed as a half matrix...
    maxv <- max(distvec)
    if (fixedNbr>1) distvec[(length(distvec)-fixedNbr*(fixedNbr-1)/2+1):length(distvec)] <- maxv ## compar among fixed points ignored below
    distmat <- as.matrix(distvec) ## constructs a symm matrix from the dist object; not optimal for memory, but otherwise would make it hard to remove columns
    distvec <- NA
    diag(distmat) <- maxv ## following code will also ignore these elements
    rownames(distmat) <- rownames(rownamedarray)
    nr <- nrow(distmat); ## will decrease
    # zutindices <- distmat
    # zutindices[TRUE, TRUE] <- seq_len(length(distmat))
    # zutordre <- order(as.vector(distmat)) ## will be updated to give the order of remainig elements with indices referring to positions in the reduced distmat
    while (nr>finallength) {
      ordre <- order(as.vector(distmat))[1:(2*(nr-finallength))] ## each distance appears twice
      # ordre <- zutordre[1:(2*(nr-finallength))]
      removand <- c()
      while (length(ordre)>0) { ## removes points by small groups... hence the total loop is slow
        rows <- ordre %% nr # numbers, not rownames
        rows[rows==0] <- nr
        cols <- (ordre-1) %/% nr
        cols <- cols+1
        ## takes out rows of fixed coords:
        if (fixedNbr>0) {tmp <- rows %w/o% (nr-fixedNbr+1):nr;} else tmp <- rows
        trows <- table(tmp) # its names are row indices, not rownames
        choix <- as.numeric(names(trows)[trows==max(trows)])  ## chosen for 'deletion'
        if (length(choix)>1) { ## if more than one is chosen, we have to choose on of them
          rowmins <- apply(distmat[choix, ], 1, min) #names of this are row indices of distmat
          minrowmins <- min(rowmins)
          namedrowindices <- which(rowmins==minrowmins) ## may be a pair of nonfixed points
          if (length(namedrowindices)>1) { ## further  chose within this pair based on other values in their rows
            tmp <- distmat[namedrowindices, ]
            atuer <- sapply(seq(nrow(tmp)), function(cc) {min(tmp[cc, tmp[cc, ]!=minrowmins])}) ## (for each row, this took the minimal element among those that are not the minimal one of the whole matrix
            rownam <- names(namedrowindices)[which.min(atuer)] ## namedrowindices has kept original row name of distmat
            choix <- which(rownames(distmat)==rownam)
          } else choix <- which(rownames(distmat)==names(namedrowindices))
        }
        ordre <- ordre[rows!=choix]
        cols <- cols[rows!=choix]
        ordre <- ordre[cols!=choix]
        removand <- append(removand, choix)
      }
      distmat <- distmat[-removand, -removand]
      nr <- nrow(distmat)
      # dif <- setdiff(zutordre, unique(c(as.vector(zutindices[removand, ]), as.vector(zutindices[, removand]))))
      # zutindices[, removand] <- NA
      # zutindices[removand, ] <- NA
      # zutindices[!is.na(zutindices)] <- seq(nr^2) ##
      # zutordre <- zutindices[dif] ## order with new indices from order with old indices
      # zutindices <- zutindices[-removand, -removand]
    }
    return(rownames(distmat))	## includes fixed points
  }
} ##end stripclosestpairs
