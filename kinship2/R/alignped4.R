# Automatically generated from all.nw using noweb
alignped4 <- function(rval, spouse, level, width, align) {
    if (is.logical(align)) align <- c(1.5, 2)  #defaults
    maxlev <- nrow(rval$nid)
    width <- max(width, rval$n+.01)   # width must be > the longest row

    n <- sum(rval$n)  # total number of subjects
    myid <- matrix(0, maxlev, ncol(rval$nid))  #number the plotting points
    for (i in 1:maxlev) {
        myid[i, rval$nid[i,]>0] <-  cumsum(c(0, rval$n))[i] + 1:rval$n[i]
        }

    # There will be one penalty for each spouse and one for each child
    npenal <- sum(spouse[rval$nid>0]) + sum(rval$fam >0) 
    pmat <- matrix(0., nrow=npenal+1, ncol=n)

    indx <- 0
    # Penalties to keep spouses close
    for (lev in 1:maxlev) {
        if (any(spouse[lev,])) {
            who <- which(spouse[lev,])
            indx <- max(indx) + 1:length(who)
            pmat[cbind(indx, myid[lev,who])] <-  sqrt(align[2])
            pmat[cbind(indx, myid[lev,who+1])] <- -sqrt(align[2])
            }
        }

    # Penalties to keep kids close to parents
    for (lev in (1:maxlev)[-1])  { # no parents at the top level
        families <- unique(rval$fam[lev,])
        families <- families[families !=0]  #0 is the 'no parent' marker
        for (i in families) {  #might be none
            who <- which(rval$fam[lev,] == i)
            k <- length(who)
            indx <- max(indx) +1:k   #one penalty per child
            penalty <- sqrt(k^(-align[1]))
            pmat[cbind(indx, myid[lev,who])] <- -penalty
            pmat[cbind(indx, myid[lev-1, rval$fam[lev,who]])] <- penalty/2
            pmat[cbind(indx, myid[lev-1, rval$fam[lev,who]+1])] <- penalty/2
            }
        }
    maxrow <- min(which(rval$n==max(rval$n)))
    pmat[nrow(pmat), myid[maxrow,1]] <- 1e-5
    ncon <- n + maxlev    # number of constraints
    cmat <- matrix(0., nrow=ncon, ncol=n)
    coff <- 0  # cumulative constraint lines so var
    dvec <- rep(1., ncon)
    for (lev in 1:maxlev) {
        nn <- rval$n[lev]
        if (nn>1) {
            for (i in 1:(nn-1)) 
                cmat[coff +i, myid[lev,i + 0:1]] <- c(-1,1)
            }

        cmat[coff+nn,   myid[lev,1]]  <- 1     #first element >=0
        dvec[coff+nn] <- 0
        cmat[coff+nn+1, myid[lev,nn]] <- -1    #last element <= width-1
        dvec[coff+nn+1] <- 1-width
        coff <- coff + nn+ 1
        }

    if (exists('solve.QP')) {
         pp <- t(pmat) %*% pmat + 1e-8 * diag(ncol(pmat))
         fit <- solve.QP(pp, rep(0., n), t(cmat), dvec)
         }
    else stop("Need the quadprog package")

    newpos <- rval$pos
    #fit <- lsei(pmat, rep(0, nrow(pmat)), G=cmat, H=dvec)
    #newpos[myid>0] <- fit$X[myid]           
    newpos[myid>0] <- fit$solution[myid]
    newpos
    }
