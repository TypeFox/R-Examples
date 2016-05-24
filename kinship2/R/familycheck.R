# Automatically generated from all.nw using noweb
# This routine checks out a family id, by trying to construct its own
#  and comparing the results
#
# The input argument "newfam" is optional: if you've already created this
#   vector for other reasons, then putting the arg in saves time.
#
# Output is a dataframe with columns:
#   famid: the family id, as entered into the data set
#   n    : number of subjects in the family
#   unrelated: number of them that appear to be unrelated to anyone else 
#          in the entire pedigree set.  This is usually marry-ins with no 
#          children (in the pedigree), and if so are not a problem.
#   split : number of unique "new" family ids.
#            if this is 0, it means that no one in this "family" is related to
#                   anyone else (not good)
#            1 = everythings is fine
#            2+= the family appears to be a set of disjoint trees.  Are you
#                 missing some of the people?
#   join : number of other families that had a unique famid, but are actually
#            joined to this one.  0 is the hope.
#
#  If there are any joins, then an attribute "join" is attached.  It will be
#   a matrix with famid as row labels, new-family-id as the columns, and
#   the number of subjects as entries.  
#
familycheck <- function(famid, id, father.id, mother.id, newfam) {
    if (is.numeric(famid) && any(is.na(famid)))
        stop ("Family id of missing not allowed")
    nfam <- length(unique(famid))

    if (missing(newfam)) newfam <- makefamid(id, father.id, mother.id)
    else if (length(newfam) != length(famid))
        stop("Invalid length for newfam")

    xtab <- table(famid, newfam)
    if (any(newfam==0)) {
        unrelated <- xtab[,1]
        xtab <- xtab[,-1, drop=FALSE] 
        ## bug fix suggested by Amanda Blackford 6/2011
      }
    else unrelated <-  rep(0, nfam)

    splits <- apply(xtab>0, 1, sum)
    joins  <- apply(xtab>0, 2, sum)

    temp <- apply((xtab>0) * outer(rep(1,nfam), joins-1), 1, sum)

    out <- data.frame(famid = dimnames(xtab)[[1]],
                      n = as.vector(table(famid)),
                      unrelated = as.vector(unrelated),
                      split = as.vector(splits),
                      join = temp,
                      row.names=1:nfam)
    if (any(joins >1)) {
      tab1 <- xtab[temp>0,]  #families with multiple outcomes
      tab1 <- tab1[,apply(tab1>0,2,sum) >0] #only keep non-zero columns
      dimnames(tab1) <- list(dimnames(tab1)[[1]], NULL)
      attr(out, 'join') <- tab1
    }
    
    out
  }


