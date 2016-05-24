# Automatically generated from all.nw using noweb


#$Log: pedigree.unrelated.q,v $
#Revision 1.2  2010/02/11 22:36:48  sinnwell
#require kinship to be loaded before use
#
#Revision 1.1  2009/11/10 19:21:52  sinnwell
#Initial revision
#
#Revision 1.1  2009/11/03 16:42:27  sinnwell
#Initial revision
#
## Authors: Dan Schaid, Shannon McDonnell
## Updated by Jason Sinnwell

pedigree.unrelated <- function(ped, avail) {
  
  # Requires: kinship function

  # Given vectors id, father, and mother for a pedigree structure,
  # and avail = vector of T/F or 1/0 for whether each subject
  # (corresponding to id vector) is available (e.g.,
  # has DNA available), determine set of maximum number
  # of unrelated available subjects from a pedigree.

  # This is a greedy algorithm that uses the kinship
  # matrix, sequentially removing rows/cols that
  # are non-zero for subjects that have the most
  # number of zero kinship coefficients (greedy
  # by choosing a row of kinship matrix that has
  # the most number of zeros, and then remove any
  # cols and their corresponding rows that are non-zero.
  # To account for ties of the count of zeros for rows,
  # a random choice is made. Hence, running this function
  # multiple times can return different sets of unrelated
  # subjects.

  id <- ped$id
  avail <- as.integer(avail)

  
  kin <- kinship(ped)
  
  ord <- order(id)
  id <- id[ord]
  avail <- as.logical(avail[ord])
  kin <- kin[ord,][,ord]

  rord <- order(runif(nrow(kin)))

  id <- id[rord]
  avail <- avail[rord]
  kin <- kin[rord,][,rord]

  id.avail <- id[avail]
  kin.avail <- kin[avail,,drop=FALSE][,avail,drop=FALSE]

  diag(kin.avail) <- 0

  while(any(kin.avail > 0))
    {
      nr <- nrow(kin.avail)
      indx <- 1:nrow(kin.avail)
      zero.count <- apply(kin.avail==0, 1, sum)
      
      mx <- max(zero.count[zero.count < nr])
      mx.zero <- indx[zero.count == mx][1]

      exclude <- indx[kin.avail[, mx.zero] > 0]

      kin.avail <- kin.avail[- exclude, , drop=FALSE][, -exclude, drop=FALSE]

    }

  choice <- sort(dimnames(kin.avail)[[1]])
  
  return(choice)
}


