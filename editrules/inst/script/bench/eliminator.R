
# Eliminate all variables one-by-one and gather some details.
#
elimInfo <- function(e, elimfun){
#    cat('R version:',R.Version()$version.string,'\n')
#    cat('editrules version:',as.character(packageVersion('editrules')),'\n')
    v <- getVars(e)

    L <- lapply(v, function(s) elimfun(e,s))
    names(L) <- v
    results <- data.frame(variable=v,t(sapply(L,`[[`,'exitstatus')),stringsAsFactors=FALSE)
    results$numcat <- sapply(editrules:::getInd(e),length)
    results$numiter <- sapply(L, function(d) sum(d$dat$nA>0))
    results$contained <- sapply(L, function(d) d$dat$nA[1])
    results <- cbind(results,t(sapply(L,function(d) colSums(d$dat[,4:7]))))
    results$maxedits <- sapply(L,function(d) max(d$dat$nA))
    results$editsreturned <- sapply(L,function(d) nrow(d$E))
    results
}




###########################################################################

# Here, eliminated rules are not removed on-the-fly
# 
algorithmA.follow <- function(E, var, ...){
    # do not bother with edits not containing var
    I <- contains(E,var)

    ind <- editrules:::getInd(E)
    J <- ind[[var]]
    # START ADMINISTRATION
    nAt <- nA <- Created <- Empty <- Subset <- Eliminated <- numeric(length(J))
    iter <- 1:length(J)
    exitstatus <-  c(nothingtoeliminate=FALSE,noiteration=FALSE)
    # nothing to eliminate...
   
    if ( sum(I) == 0 ){ 
        ### ADMIN
          exitstatus['nothingtoeliminate'] <- TRUE
        ### END ADMIN
        return(
            list(E=E,
                exitstatus=exitstatus,
                dat =  data.frame(
                    iter=iter, 
                    nAt=nAt, 
                    nA=nA, 
                    Created=Created,
                    Empty=Empty,
                    Subset=Subset,
                    Eliminated=Eliminated)
            )
        ) 
    }

    # if elimination is not possible... (at least one category cannot be resolved)
    if ( any(colSums(E[I,J,drop=FALSE])==0) ){ 
        ### ADMIN
          exitstatus['noiteration'] <- TRUE
        ### END ADMIN
        return(
            list(E=E[!I,],
                exitstatus = exitstatus,
                dat =  data.frame(
                    iter=iter, 
                    nAt=nAt, 
                    nA=nA, 
                    Created=Created,
                    Empty=Empty,
                    Subset=Subset,
                    Eliminated=Eliminated)
            )
        )
    }

    A <- editrules:::getArr(E)
    At <- A[!I,,drop=FALSE]
    A <- A[I,,drop=FALSE]

    iii <- 0   
    k <- integer(0)
    nel <- 0
    for ( j in J ){
        ### ADMIN
          iii <- iii + 1
          nA[iii] <- nrow(A)
          nAt[iii] <- nrow(At)
        ### END ADMIN
        if (nrow(A) <= 1) break
        aPlus <- A[ A[,j],,drop=FALSE]
        aMin  <- A[!A[,j],,drop=FALSE]  
        n <- nrow(aPlus)
        m <- nrow(aMin)
        if ( m == 0 ) next
        if ( n == 0 ) break
        B <- array(FALSE,dim=c(n*m,ncol(A)))
        I1 <- rep(1:n,times=m)
        I2 <- rep(1:m,each=n)
        B[I1,J] <-  aPlus[I1,J,drop=FALSE] | aMin[I2,J,drop=FALSE]
        B[I1,-J] <- aPlus[I1,-J,drop=FALSE] & aMin[I2,-J,drop=FALSE]
        ### ADMIN 
          Created[iii] <- nrow(B)
        ### END ADMIN
        B <- rbind(B,aPlus)
        iEmpty <- editrules:::isRedundant.boolmat(B,ind) 
        A <- B[!iEmpty,,drop=FALSE]
        iSubset <- editrules:::isSubset.boolmat(A)  
        if ( any(iSubset) ) A <- A[!iSubset,,drop=FALSE]
        iSubset2 <- editrules:::isSubsetWrt.boolmat(A,At)
        if ( any(iSubset2) ) A <- A[!iSubset2,,drop=FALSE]
        ### ADMIN
            el <- apply(A[,J,drop=FALSE],1,all) 
            Subset <- sum(iSubset) + sum(iSubset2)
            Empty[iii] <- sum(iEmpty)
            Eliminated <- sum(el) #- nel
            nel <- sum(el)
        ### END ADMIN
    }
    E <- editrules:::neweditarray(rbind(At,A),ind=ind,
        sep=editrules:::getSep(E), 
        levels=editrules:::getlevels(E))
    # remove edits containing the variable
    E <- E[!contains(E,var),,drop=FALSE]
    list(E=E,
        exitstatus = exitstatus,
        dat =  data.frame(
            iter=iter, 
            nAt=nAt, 
            nA=nA, 
            Created=Created,
            Empty=Empty,
            Subset=Subset,
            Eliminated=Eliminated)
    )
}

# like eliminate.editarray in the package, but with some extra administration.
# 
algorithmB.follow <- function(E, var, ...){
    # do not bother with edits not containing var
    I <- contains(E,var)

    ind <- editrules:::getInd(E)
    J <- ind[[var]]
    # START ADMINISTRATION
    nAt <- nA <- Created <- Empty <- Subset <- Eliminated <- numeric(length(J))
    iter <- 1:length(J)
    exitstatus <-  c(nothingtoeliminate=FALSE,noiteration=FALSE)
    # nothing to eliminate...
   
    if ( sum(I) == 0 ){ 
        ### ADMIN
          exitstatus['nothingtoeliminate'] <- TRUE
        ### END ADMIN
        return(
            list(E=E,
                exitstatus=exitstatus,
                dat =  data.frame(
                    iter=iter, 
                    nAt=nAt, 
                    nA=nA, 
                    Created=Created,
                    Empty=Empty,
                    Subset=Subset,
                    Eliminated=Eliminated)
            )
        ) 
    }

    # if elimination is not possible... (at least one category cannot be resolved)
    if ( any(colSums(E[I,J,drop=FALSE])==0) ){ 
        ### ADMIN
          exitstatus['noiteration'] <- TRUE
        ### END ADMIN
        return(
            list(E=E[!I,],
                exitstatus = exitstatus,
                dat =  data.frame(
                    iter=iter, 
                    nAt=nAt, 
                    nA=nA, 
                    Created=Created,
                    Empty=Empty,
                    Subset=Subset,
                    Eliminated=Eliminated)
            )
        )
    }

    A <- editrules:::getArr(E)
    At <- A[!I,,drop=FALSE]
    A <- A[I,,drop=FALSE]

    iii <- 0   
    k <- integer(0)
    for ( j in J ){
        ### ADMIN
          iii <- iii + 1
          nA[iii] <- nrow(A)
          nAt[iii] <- nrow(At)
        ### END ADMIN
        if (nrow(A) <= 1) break
        aPlus <- A[ A[,j],,drop=FALSE]
        aMin  <- A[!A[,j],,drop=FALSE]  
        n <- nrow(aPlus)
        m <- nrow(aMin)
        if ( m == 0 ) next
        if ( n == 0 ) break
        B <- array(FALSE,dim=c(n*m,ncol(A)))
        I1 <- rep(1:n,times=m)
        I2 <- rep(1:m,each=n)
        B[I1,J] <-  aPlus[I1,J,drop=FALSE] | aMin[I2,J,drop=FALSE]
        B[I1,-J] <- aPlus[I1,-J,drop=FALSE] & aMin[I2,-J,drop=FALSE]
        ### ADMIN 
          Created[iii] <- nrow(B)
        ### END ADMIN
        B <- rbind(B,aPlus)
        iEmpty <- editrules:::isRedundant.boolmat(B,ind) 
        A <- B[!iEmpty,,drop=FALSE]
        iSubset <- editrules:::isSubset.boolmat(A)  
        if ( any(iSubset) ) A <- A[!iSubset,,drop=FALSE]
        iSubset2 <- editrules:::isSubsetWrt.boolmat(A,At)
        if ( any(iSubset2) ) A <- A[!iSubset2,,drop=FALSE]
        el <- apply(A[,J,drop=FALSE],1,all)
        ### ADMIN
            Subset <- sum(iSubset) + sum(iSubset2)
            Empty[iii] <- sum(iEmpty)
            Eliminated <- sum(el)
        ### END ADMIN
        At <- rbind(At,A[el,,drop=FALSE])     
        A <- A[!el,,drop=FALSE]
    }

    E <- editrules:::neweditarray(At,ind=ind,
        sep=editrules:::getSep(E), 
        levels=editrules:::getlevels(E))
    E <- E[!contains(E,'var'),,drop=FALSE]
    list(E=E,
        exitstatus = exitstatus,
        dat =  data.frame(
            iter=iter, 
            nAt=nAt, 
            nA=nA, 
            Created=Created,
            Empty=Empty,
            Subset=Subset,
            Eliminated=Eliminated)
    )
}





###########################################################################
# ELIMINATION ALGORITHMS, no administration overhead
#
algorithmA <- function(E, var, ...){
    # do not bother with edits not containing var
    I <- contains(E,var)

    ind <- editrules:::getInd(E)
    J <- ind[[var]]

    # nothing to eliminate
    if ( sum(I) == 0 ) return(E)
    # elimination not necessary (at least 1 category cannot be resolved)
    if ( any(colSums(E[I,J,drop=FALSE])==0)  ) return(E[!I,,drop=FALSE])

    A <- editrules:::getArr(E)
    At <- A[!I,,drop=FALSE]
    A <- A[I,,drop=FALSE]

    iii <- 0   
    k <- integer(0)
    nel <- 0
    for ( j in J ){
        if (nrow(A) <= 1) break
        aPlus <- A[ A[,j],,drop=FALSE]
        aMin  <- A[!A[,j],,drop=FALSE]  
        n <- nrow(aPlus)
        m <- nrow(aMin)
        if ( m == 0 ) next
        if ( n == 0 ) break
        B <- array(FALSE,dim=c(n*m,ncol(A)))
        I1 <- rep(1:n,times=m)
        I2 <- rep(1:m,each=n)
        B[I1,J] <-  aPlus[I1,J,drop=FALSE] | aMin[I2,J,drop=FALSE]
        B[I1,-J] <- aPlus[I1,-J,drop=FALSE] & aMin[I2,-J,drop=FALSE]
        B <- rbind(B,aPlus)
        iEmpty <- editrules:::isRedundant.boolmat(B,ind) 
        A <- B[!iEmpty,,drop=FALSE]
        iSubset <- editrules:::isSubset.boolmat(A)  
        if ( any(iSubset) ) A <- A[!iSubset,,drop=FALSE]
        iSubset2 <- editrules:::isSubsetWrt.boolmat(A,At)
        if ( any(iSubset2) ) A <- A[!iSubset2,,drop=FALSE]
    }
    E <- editrules:::neweditarray(rbind(At,A),ind=ind,
        sep=editrules:::getSep(E), 
        levels=editrules:::getlevels(E))
    # remove edits containing the variable
    E[!contains(E,var),]
}

# like eliminate.editarray in the package, but with some extra administration.
# 
algorithmB <- function(E, var, ...){
    # do not bother with edits not containing var
    I <- contains(E,var)

    ind <- editrules:::getInd(E)
    J <- ind[[var]]
    # nothing to eliminate...
    if ( sum(I) == 0 ) return(E)
    # if elimination is not possible... (at least one category cannot be resolved)
    if ( any(colSums(E[I,J,drop=FALSE])==0) ) return(E[!I,,drop=FALSE])

    A <- editrules:::getArr(E)
    At <- A[!I,,drop=FALSE]
    A <- A[I,,drop=FALSE]

    iii <- 0   
    k <- integer(0)
    for ( j in J ){
        if (nrow(A) <= 1) break
        aPlus <- A[ A[,j],,drop=FALSE]
        aMin  <- A[!A[,j],,drop=FALSE]  
        n <- nrow(aPlus)
        m <- nrow(aMin)
        if ( m == 0 ) next
        if ( n == 0 ) break
        B <- array(FALSE,dim=c(n*m,ncol(A)))
        I1 <- rep(1:n,times=m)
        I2 <- rep(1:m,each=n)
        B[I1,J] <-  aPlus[I1,J,drop=FALSE] | aMin[I2,J,drop=FALSE]
        B[I1,-J] <- aPlus[I1,-J,drop=FALSE] & aMin[I2,-J,drop=FALSE]
        B <- rbind(B,aPlus)
        iEmpty <- editrules:::isRedundant.boolmat(B,ind) 
        A <- B[!iEmpty,,drop=FALSE]
        iSubset <- editrules:::isSubset.boolmat(A)  
        if ( any(iSubset) ) A <- A[!iSubset,,drop=FALSE]
        iSubset2 <- editrules:::isSubsetWrt.boolmat(A,At)
        if ( any(iSubset2) ) A <- A[!iSubset2,,drop=FALSE]
        el <- apply(A[,J,drop=FALSE],1,all)
        At <- rbind(At,A[el,,drop=FALSE])     
        A <- A[!el,,drop=FALSE]
    }

    E <- editrules:::neweditarray(At,ind=ind,
        sep=editrules:::getSep(E), 
        levels=editrules:::getlevels(E))
    # remove edits containing the variable
    E[!contains(E,var),,drop=FALSE]
}






