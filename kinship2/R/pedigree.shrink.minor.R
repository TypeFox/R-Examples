# Automatically generated from all.nw using noweb

#$Log: pedigree.shrink.minor.q,v $
#Revision 1.5  2009/11/19 18:10:26  sinnwell
#F to FALSE
#
#Revision 1.4  2009/11/19 14:57:13  sinnwell
#*** empty log message ***
#
#Revision 1.3  2009/11/17 23:11:41  sinnwell
#*** empty log message ***
#
#Revision 1.1  2008/07/16 20:22:55  sinnwell
#Initial revision
#


is.parent <- function(id, findex, mindex){
  # determine subjects who are parents
  # assume input of father/mother indices, not ids
  
  father <- mother <- rep(0, length(id))
  father[findex>0] <- id[findex]
  mother[mindex>0] <- id[mindex]
  
  isFather <- !is.na(match(id, unique(father[father!=0])))
  isMother <- !is.na(match(id, unique(mother[mother!=0])))
  isParent <- isFather |isMother
  return(isParent)
}

is.founder <- function(mother, father){
  check <- (father==0) & (mother==0)
  return(check)
}


is.disconnected <- function(id, findex, mindex)
{

  # check to see if any subjects are disconnected in pedigree by checking for
  # kinship = 0 for all subjects excluding self
  father <- id[findex]
  mother <- id[mindex]  
  kinMat <- kinship(id, father, mother)
  diag(kinMat) <- 0
  disconnected <- apply(kinMat==0.0, 1, all)

  return(disconnected)
}

