#
# An attempt at a founder score
#
fscore.pedigree <- function(ped) {
    n <- length(ped$depth)
    nid <- 1:n
    dad<- match(ped$dadid, ped$id, nomatch=0)
    mom<- match(ped$momid, ped$id, nomatch=0) 
    sex<- 1*(ped$sex=='female')
    level <- ped$depth +1
	


    # Currently, the alignment routine requires that you have either
    #  0 parents or 2 parents, not 1.
    if (any(dad==0 & mom>0) || any(dad>0 & mom==0))
            stop("Everyone must have 0 parents or 2 parents, not just one")
    

    pfun1 <- function(id, dad, mom, sex) { 
        # This function returns a list of all founders above a given id
        if (dad[id]==0 && mom[id]==0) {
            if (sex[id]==1) return(id)  #I am a founder
            else return(NULL)
            }
        else return(c(Recall(dad[id], dad, mom, sex), 
                      Recall(mom[id], dad, mom, sex)))
        }

    #
    # Now, find all marriages 
    #  (This code is essentially the same as that in align.pedigree)
    #   Mothers are row 1, fathers row 2, each column a spouse pair
    # Either having children, having a hint, or a spousal
    #   relationship gets you on the list (that is, any of the things that
    #   would cause this marriage to part of the plotted pedigree).
    # Start with a hashlist: fatherid * max(id) + mother
    hash <- (dad*n + mom)[mom>0 & dad>0]
    sorder <- ped$hints[,2]
    temp1 <- abs(sorder)
    temp1 <- ifelse(ped$sex=='male', (1:n)*n + temp1, 1:n + n*temp1)
    hash <- c(hash, temp1[sorder!=0]) #those with a spouse hint
    if (!is.null(ped$relation) && any(ped$relation[,3]==4)) {
	who <- (ped$relation[,3]==4)  # add spouses from relationship list
	indx <- ped$relation[who,1]   # id of the first spouse
	temp1 <- ifelse(ped$sex[indx]=='male', n*indx + ped$relation[who,2],
			                       indx + n*ped$relation[who,2])
	hash <- c(temp1, hash) #being first is important -- it controls plot
	                       # order per the documentation
	}

    hash <- unique(hash)              #eliminate duplicates
    hash <- hash -1                   #change to range of 0 to n-1 (for %%)
    spouselist <- rbind(1+ hash%%n, floor(hash/n))  # uncompress it
    temp <- matrix(mom[spouselist] + dad[spouselist], nrow=2)
    spouselist <- spouselist[, colSums(temp) > 0]  # toss founders off list

    # Create the founder sets for each spouse pair
    #  It will be a list of npair elements, each of which is a list
    #  containing 2 vectors
    npair <- ncol(spouselist)
    flist <- vector('list', npair)
    for (i in 1:npair) {
        flist[[i]] <- list(f1=pfun1(spouselist[1,i], dad, mom, sex),
                           f2=pfun1(spouselist[2,i], dad, mom, sex))
        }

    pfun2 <- function(flist, ord) {
        # This function computes the badness score for the founders
        #  of each parental pair, along with a report of which of the
        #  4 terms is the "worst".
        tfun <- function(x, ord) {
            t1 <- ord[x[[1]]]
            t2 <- ord[x[[2]]]
            c(sum(t2 < max(t1)), sum(t1 > min(t2)),
              sum(t1 < max(t2)), sum(t2 > min(t1)))
            }
        temp <- matrix(unlist(lapply(flist, tfun, ord=ord)), nrow=4)
        temp.min <- apply(temp,2,min)
        temp.which <- apply(temp== matrix(rep(temp.min,each=4), nrow=4), 2, 
                            function(x) min(c(which(x))))
        cbind(temp.min, temp.which)
        }

    temp <- pfun2(flist, 1:n)
    browser()
    npair
    }

    
