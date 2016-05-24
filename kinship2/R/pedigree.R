# Automatically generated from all.nw using noweb
pedigree <- function(id, dadid, momid, sex, affected, status, relation,
                     famid, missid) {
    n <- length(id)
    if (length(momid) != n) stop("Mismatched lengths, id and momid")
    if (length(dadid) != n) stop("Mismatched lengths, id and momid")
    if (length(sex  ) != n) stop("Mismatched lengths, id and sex")


    # Allow for character/numeric/factor in the sex variable
    if(is.factor(sex))
            sex <- as.character(sex)
    codes <- c("male","female", "unknown", "terminated")
    if(is.character(sex)) sex<- charmatch(casefold(sex, upper = FALSE), codes, nomatch = 3)        

    # assume either 0/1/2/4 =  female/male/unknown/term, or 1/2/3/4
    #  if only 1/2 assume no unknowns
    if(min(sex) == 0) {
      warning("Sex values contain 0, but expected codes 1-4.\n Setting 0=male, 1=female, 2=unknown, 3=terminated. \n")
      sex <- sex + 1
    }
    sex <- ifelse(sex < 1 | sex > 4, 3, sex)
    if(all(sex > 2))
            stop("Invalid values for 'sex'")
        else if(mean(sex == 3) > 0.25)
                warning("More than 25% of the gender values are 'unknown'")
    sex <- factor(sex, 1:4, labels = codes)
    if (missing(missid)) {
      if (is.numeric(id)) missid <- 0
      else missid <- ""
    }
    ## Don't allow missing id values
    if (any(is.na(id))) stop("Missing value for the id variable")
    if (!is.numeric(id)) {
      id <- as.character(id)
      dadid <- as.character(dadid)
      momid <- as.character(momid)
      missid <- as.character(missid)
      if (length(grep('^ *$', id)) > 0)
        stop("A blank or empty string is not allowed as the id variable")
    }

    nofather <- (is.na(dadid) | dadid==missid)
    nomother <- (is.na(momid) | momid==missid)

    if (!missing(famid)) {
        if(is.factor(famid)) famid <- as.character(famid)
        
        if (any(is.na(famid))) stop("The family id cannot contain missing values")
        if (is.factor(famid) || is.character(famid)) {
            if (length(grep('^ *$', famid)) > 0)
                stop("The family id cannot be a blank or empty string")
            }
        #Make a temporary new id from the family and subject pair
        oldid <-id
        id <- paste(as.character(famid), as.character(id), sep='/')
        dadid <- paste(as.character(famid), as.character(dadid), sep='/')
        momid <- paste(as.character(famid), as.character(momid), sep='/')
        }

    if (any(duplicated(id))) {
      duplist <- id[duplicated(id)]
      msg.n <- min(length(duplist), 6)
      stop(paste("Duplicate subject id:", duplist[1:msg.n]))
    }
    findex <- match(dadid, id, nomatch = 0)
    mindex <- match(momid, id, nomatch = 0)
        
    ## if all non-matches of fathers/mothers are on same people, those are founders
    maybe.founder <- (mindex==0 & !nomother) & (findex==0 & !nofather)    
    founder.momid <- unique(momid[maybe.founder])
    founder.dadid <- unique(dadid[maybe.founder])
    if(length(founder.dadid)==1 & length(founder.momid)==1 && founder.momid == founder.dadid) {
      warning("Using ", founder.momid, " as momid and dadid of founders.\n")
      nofather[maybe.founder] <- nomother[maybe.founder] <- TRUE
    } 

    if(sum(nofather & nomother)<2) 
      stop("pedigree must have at least 2 founders\n")
        
    if(any(sex[findex] != "male")) {
        who <- unique((id[findex])[sex[findex] != "male"])
        msg.n <- 1:min(5, length(who))  #Don't list a zillion
        stop(paste("Id not male, but is a father:", 
                   paste(who[msg.n], collapse= " ")))
        }

    if (any(findex==0 & !nofather)) {
        who <- dadid[which(findex==0 & !nofather)]
        msg.n <- 1:min(5, length(who))  #Don't list a zillion
        stop(paste("Value of 'dadid' not found in the id list", 
                   paste(who[msg.n], collapse= " ")))
        }
        
    if(any(sex[mindex] != "female")) {
        who <- unique((id[mindex])[sex[mindex] != "female"])
        msg.n <- 1:min(5, length(who))
        stop(paste("Id not female, but is a mother:", 
                   paste(who[msg.n], collapse = " ")))
        }

    if (any(mindex==0 & !nomother)) {
        who <- momid[which(mindex==0 & !nomother)]
        msg.n <- 1:min(5, length(who))  #Don't list a zillion
        stop(paste("Value of 'momid' not found in the id list", 
                   paste(who[msg.n], collapse= " ")))
        }

    if (any(mindex==0 & findex!=0) || any(mindex!=0 & findex==0)) {
        who <- id[which((mindex==0 & findex!=0) |(mindex!=0 & findex==0))] 
        msg.n <- 1:min(5, length(who))  #Don't list a zillion
        stop(paste("Subjects must have both a father and mother, or have neither",
                   paste(who[msg.n], collapse= " ")))
    }

    if (!missing(famid)) {
        if (any(famid[mindex] != famid[mindex>0])) {
            who <- (id[mindex>0])[famid[mindex] != famid[mindex>0]]
            msg.n <- 1:min(5, length(who))
            stop(paste("Mother's family != subject's family", 
                       paste(who[msg.n], collapse=" ")))
            }
        if (any(famid[findex] != famid[findex>0])) {
            who <- (id[findex>0])[famid[findex] != famid[findex>0]]
            msg.n <- 1:min(5, length(who))
            stop(paste("Father's family != subject's family", 
                       paste(who[msg.n], collapse=" ")))
            }
        }
    if (missing(famid))
        temp <- list(id = id, findex=findex, mindex=mindex, sex=sex)
    else temp<- list(famid=famid, id=oldid, findex=findex, mindex=mindex, 
                     sex=sex)
    if (!missing(affected) && !is.null(affected)) {
        if (is.matrix(affected)){
            if (nrow(affected) != n) stop("Wrong number of rows in affected")
            if (is.logical(affected)) affected <- 1* affected
            } 
        else {
            if (length(affected) != n)
                stop("Wrong length for affected")

            if (is.logical(affected)) affected <- as.numeric(affected)
            if (is.factor(affected))  affected <- as.numeric(affected) -1
            }
        if(max(affected, na.rm=TRUE) > min(affected, na.rm=TRUE)) 
          affected <- affected - min(affected, na.rm=TRUE)
        if (!all(affected==0 | affected==1 | is.na(affected)))
                    stop("Invalid code for affected status")
        temp$affected <- affected
        }
    if(!missing(status) && !is.null(status)) {
        if(length(status) != n)
            stop("Wrong length for affected")
        if (is.logical(status)) status <- as.integer(status)
        if(any(status != 0 & status != 1))
            stop("Invalid status code")
        temp$status <- status
        }
    if (!missing(relation) && !is.null(relation)) {
        if (!missing(famid)) {
            if (is.matrix(relation)) {
                if (ncol(relation) != 4) 
                    stop("Relation matrix must have 3 columns + famid")
                id1 <- relation[,1]
                id2 <- relation[,2]
                code <- relation[,3]
                famid <- relation[,4]
                }
            else if (is.data.frame(relation)) {
                id1 <- relation$id1
                id2 <- relation$id2
                code <- relation$code
                famid <- relation$famid
                if (is.null(id1) || is.null(id2) || is.null(code) ||is.null(famid)) 
                stop("Relation data must have id1, id2, family id and code")
                }
            else stop("Relation argument must be a matrix or a dataframe")
            }
        else {
            if (is.matrix(relation)) {
                if (ncol(relation) != 3) 
                    stop("Relation matrix must have 3 columns")
                id1 <- relation[,1]
                id2 <- relation[,2]
                code <- relation[,3]
                }
            else if (is.data.frame(relation)) {
                id1 <- relation$id1
                id2 <- relation$id2
                code <- relation$code
                if (is.null(id1) || is.null(id2) || is.null(code)) 
                    stop("Relation data frame must have id1, id2, and code")
                }
            else stop("Relation argument must be a matrix or a list")
            }
        
        if (!is.numeric(code))
            code <- match(code, c("MZ twin", "DZ twin", "UZ twin", "spouse"))
        else code <- factor(code, levels=1:4,
                            labels=c("MZ twin", "DZ twin", "UZ twin", "spouse"))
        if (any(is.na(code)))
            stop("Invalid relationship code")
         
        # Is everyone in this relationship in the pedigree?
        if (!missing(famid)) {
            temp1 <- match(paste(as.character(famid), as.character(id1), sep='/'), 
                           id, nomatch=0)
            temp2 <- match(paste(as.character(famid), as.character(id2), sep='/'),
                           id, nomatch=0)
          }
        else {
            temp1 <- match(id1, id, nomatch=0)
            temp2 <- match(id2, id, nomatch=0)
          }
        
        if (any(temp1==0 | temp2==0))
            stop("Subjects in relationships that are not in the pedigree")
        if (any(temp1==temp2)) {
            who <- temp1[temp1==temp2]
            stop(paste("Subject", id[who], "is their own spouse or twin"))
            }

        # Check, are the twins really twins?
        ncode <- as.numeric(code)
        if (any(ncode<3)) {
            twins <- (ncode<3)
            if (any(momid[temp1[twins]] != momid[temp2[twins]]))
                stop("Twins found with different mothers")
            if (any(dadid[temp1[twins]] != dadid[temp2[twins]]))
                stop("Twins found with different fathers")
            }
        # Check, are the monozygote twins the same gender?
        if (any(code=="MZ twin")) {
            mztwins <- (code=="MZ twin")
            if (any(sex[temp1[mztwins]] != sex[temp2[mztwins]]))
                stop("MZ Twins with different genders")
            }

        ##Use id index as indx1 and indx2
        if (!missing(famid)) {
            temp$relation <- data.frame(famid=famid, indx1=temp1, indx2=temp2, code=code)
            
          }
        else temp$relation <- data.frame(indx1=temp1, indx2=temp2, code=code)
        }
    if (missing(famid)) class(temp) <- 'pedigree'
    else class(temp) <- 'pedigreeList'
    temp
    }
"[.pedigreeList" <- function(x, ..., drop=F) {
    if (length(list(...)) != 1) stop ("Only 1 subscript allowed")
    ufam <- unique(x$famid)
    if (is.factor(..1) || is.character(..1)) indx <- ufam[match(..1, ufam)]
    else indx <- ufam[..1]
        
    if (any(is.na(indx))) 
            stop(paste("Familiy", (..1[is.na(indx)])[1], "not found"))

    keep <- which(x$famid %in% indx)  #which rows to keep
    for (i in c('id', 'famid', 'sex'))
        x[[i]] <- (x[[i]])[keep]
    
    kept.rows <- (1:length(x$findex))[keep]
    x$findex <- match(x$findex[keep], kept.rows, nomatch=0)
    x$mindex <- match(x$mindex[keep], kept.rows, nomatch=0)
    
    #optional components
    if (!is.null(x$status)) x$status <- x$status[keep]
    if (!is.null(x$affected)) {
        if (is.matrix(x$affected)) x$affected <- x$affected[keep,,drop=FALSE]
        else x$affected <- x$affected[keep]
        }
    if (!is.null(x$relation)) {
        keep <- !is.na(match(x$relation$famid, indx))
       if (any(keep)) {
            x$relation <- x$relation[keep,,drop=FALSE]
            ##Update twin id indexes
            x$relation$indx1 <- match(x$relation$indx1, kept.rows, nomatch=0)
            x$relation$indx2 <- match(x$relation$indx2, kept.rows, nomatch=0)
            ##If only one family chosen, remove famid
            if (length(indx)==1) {x$relation$famid <- NULL}
            }
        else x$relation <- NULL  # No relations matrix elements for this family
        }
    
    if (length(indx)==1)  class(x) <- 'pedigree'  #only one family chosen
    else class(x) <- 'pedigreeList'
    x
    }
"[.pedigree" <- function(x, ..., drop=F) {
    if (length(list(...)) != 1) stop ("Only 1 subscript allowed")
    if (is.character(..1) || is.factor(..1)) i <- match(..1, x$id)
    else i <- (1:length(x$id))[..1]
    
    if (any(is.na(i))) paste("Subject", ..1[which(is.na(i))][1], "not found")

    z <- list(id=x$id[i],findex=match(x$findex[i], i, nomatch=0),
              mindex=match(x$mindex[i], i, nomatch=0),
             sex=x$sex[i])
    if (!is.null(x$affected)) {
        if (is.matrix(x$affected)) z$affected <- x$affected[i,, drop=F]
        else z$affected <- x$affected[i]
    }
    if (!is.null(x$famid)) z$famid <- x$famid[i]

   
    if (!is.null(x$relation)) {
      indx1 <- match(x$relation$indx1, i, nomatch=0)
      indx2 <- match(x$relation$indx2, i, nomatch=0)
      keep <- (indx1 >0 & indx2 >0)  #keep only if both id's are kept
      if (any(keep)) {
        z$relation <- x$relation[keep,,drop=FALSE]
        z$relation$indx1 <- indx1[keep]
        z$relation$indx2 <- indx2[keep]
      }
    }
    
    if (!is.null(x$hints)) {
        temp <- list(order= x$hints$order[i])
        if (!is.null(x$hints$spouse)) {
            indx1 <- match(x$hints$spouse[,1], i, nomatch=0)
            indx2 <- match(x$hints$spouse[,2], i, nomatch=0)
            keep <- (indx1 >0 & indx2 >0)  #keep only if both id's are kept
            if (any(keep))
                temp$spouse <- cbind(indx1[keep], indx2[keep],
                                     x$hints$spouse[keep,3])
            }
        z$hints <- temp
        }

    if (any(z$findex==0 & z$mindex>0) | any(z$findex>0 & z$mindex==0))
        stop("A subpedigree cannot choose only one parent of a subject")
    class(z) <- 'pedigree'
    z
    }
dim.pedigree <- function(x) length(x$id)
