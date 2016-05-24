# Read in an ibd file, and convert it to a bdsmatrix object
#
# The real work is essentially the same as makefamid -- we need to
#  figure out who makes up a family block
#
#  Each row of the data is a pair of id's, followed by the ibd value.
#  The optional idlist is an alternate dimnames

bdsmatrix.ibd <- function(id1, id2, x, idmap, diagonal) {
    nc <- ncol(id1)
    if (length(nc) && nc==3) {
        id2 <- id1[,2]
        x <- id1[,3]
        id1 <- id1[,1]
        }

    # The line below was later dropped -- someone might be
    #  making an entirely differnt type of matrix
    #    if (any(x <0 | x >1)) stop ("Invalid X values")
    keep <- (x != 0)
    if (!all(keep)) {
        id1 <- id1[keep]
        id2 <- id2[keep]
        x   <- x[keep]
        }

    # assign small integers to each
    idlist <- sort(unique(c(id1, id2)))
    if (missing(idmap)) {
        idmap <- idlist
        }
    else {
        temp <- ncol(idmap)
        if (is.null(temp) || temp !=2) stop("idmap must have 2 columns")
        temp <- match(idlist, idmap[,1])
        if (any(is.na(temp)))
            stop("Values appear in id1 or id2 that are not in idmap")
        idmap <- idmap[temp,2]
        }

    id1 <- match(id1,idlist)
    id2 <- match(id2, idlist)
    if (any(is.na(id1) |is.na(id2))) stop("idlist is not complete")
    
    # make sure the diagonal element is correct
    if (any(id1==id2)) {
        temp <- range(x[id1==id2])
        if (temp[1] != temp[2])
            warning("X values for the diagonal are not constant")
        temp <- median(x[id1==id2])
        if (!missing(diagonal) && diagonal!= temp)
            warning("Specified diagonal value disagrees with the data")
        if (missing(diagonal)) diagonal <- temp
        }
    else {
	if (missing(diagonal)) 
	    stop(paste("No diagonal elements in the data,", 
		       "and no diagonal argument was given"))
	}

    # 
    # If "diagonal" was specified, ensure that everyone is in the final output
    #  by adding a dummy line for them into the data set
    #
    if (!missing(diagonal)) {
	n <- length(idlist)
	id1 <- c(id1, 1:n)
	id2 <- c(id2, 1:n)
	x   <- c(x, rep(diagonal, n))
	}

    # put the smaller first in the list 
    #  remember, the output is a symmetric matrix
    temp <- pmin(id1, id2)
    id2  <- pmax(id1, id2)
    id1  <- temp

    #
    # Remove duplicate data. Note that if the input had
    # two entries for an element, say (1,2,10) and (2,1,12), then
    #   this will remove the latter, and never notice the inconsistent
    #   data values of 10 vs 12.
    dup <- duplicated(cbind(id1, id2))
    if (any(dup)) {
	id1 <- id1[!dup]
	id2 <- id2[!dup]
	x   <- x[!dup]
	}

    # Now, finally we get to go to work
    # Basic algorithm: iteratively set famid=min(id[members of family])
    #   Initially, everyone is a unique family id
    #   At each step, compare them to the family id's of "sibs"
    # I really don't think it will take many iterations -- test cases
    #   are all 3-4.  Worst case is a tri-diagonal submatrix of dimension
    #   k, where it takes k-1 iterations.  For the large breast data
    #   pedigree: on iteration 1 all of the blood relatives of each family
    #   collapse to a single id, and all the marry-ins with children map
    #   to the old id of that blood relative.  On the second iteration
    #   all the ids are final, and then one more for it to recognize that
    #   it is done.  
    #
    famid <- 1:length(idlist) 
    for (i in 1:length(idlist)) {
        newfam <- tapply(famid[c(id1,id2)], famid[c(id2,id1)], min)
	indx <- as.numeric(names(newfam))
	# at this point indx= old family id, newfam = new family id
	#  they will differ for families that are about to be merged
	if (all(indx == newfam)) break
        famid <-  newfam[match(famid, indx)]  #give everyone their new id
       }
    if (i>= length(idlist)) stop("Routine failed with an infinite loop")

    #
    # Now build a bdsmatrix
    #   newid will be the dimname
    # The remaining routine shares a lot with makekinship.
    newid <- idmap  #gives it the right length and class, to start
    
    counts <- table(famid)
    famlist <- sort(unique(famid))  #labels of the counts vector
    unrelated <- (counts==1)
    if (any(unrelated)) {
        nzero <- sum(unrelated)
	who <- !is.na(match(famid, famlist[unrelated]))
        newid[1:nzero] <- idmap[who]
        famlist <- famlist[!unrelated]
        counts <- counts[!unrelated]
        cumcount <- cumsum(counts) + nzero
        }
    else {
        cumcount <- cumsum(counts)
        nzero <- 0
        }

    blockn <- counts*(counts+1)/2   #size of storage for each block
    n2 <- sum(blockn)       # total amount needed
    bdata <- double(n2)
    j <- cumsum(blockn)     
    for (i in 1:length(counts)) {
        who <- (famid == famlist[i])
        n <- counts[i]           #number of people in this "family"
        #rows of data which apply
        whichrows <- !is.na(match(famid[id1], famlist[i])) 
        whichid <- sort(unique(c(id1[whichrows], id2[whichrows]))) #member ids
        fid1 <- match(id1[whichrows], whichid)
        fid2 <- match(id2[whichrows], whichid)
        temp <- matrix(0.0, n, n)
        temp[cbind(fid1, fid2)] <- x[whichrows]
        temp[cbind(fid2, fid1)] <- x[whichrows]
        diag(temp) <- diagonal
        
        bdata[seq(to=j[i], length=blockn[i])] <- temp[row(temp)>=col(temp)]
        newid[seq(to=cumcount[i], length=counts[i])] <- idmap[who]
        }

    bdsmatrix(blocksize=c(rep(1,nzero), counts),
              blocks =  c(rep(diagonal,nzero), bdata),
              dimnames=list(newid, newid))
    }


