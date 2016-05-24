classAgreement <-
    function (tab, match.names=FALSE)
{
    n <- sum(tab)
    ni <- apply(tab, 1, sum)
    nj <- apply(tab, 2, sum)

    ## patch for matching factors
    if (match.names && !is.null(dimnames(tab))) {
        lev <- intersect (colnames (tab), rownames(tab))
        p0 <- sum(diag(tab[lev,lev]))/n
        pc <- sum(ni[lev] * nj[lev])/n^2
    } else { # cutoff larger dimension
        m <- min(length(ni), length(nj))
        p0 <- sum(diag(tab[1:m, 1:m]))/n
        pc <- sum((ni[1:m] / n) * (nj[1:m] / n))
    }
    n2 <- choose(n, 2)
    rand <- 1 + (sum(tab^2) - (sum(ni^2) + sum(nj^2))/2)/n2
    nis2 <- sum(choose(ni[ni > 1], 2))
    njs2 <- sum(choose(nj[nj > 1], 2))
    crand <- (sum(choose(tab[tab > 1], 2)) -
                  (nis2 * njs2)/n2)/((nis2 + njs2)/2 -
                                         (nis2 * njs2)/n2)
    list(diag = p0, kappa = (p0 - pc)/(1 - pc), rand = rand,
         crand = crand)
}

matchClasses <- function(tab, method = "rowmax", iter=1, maxexact=9,
                         verbose=TRUE){

    methods <- c("rowmax", "greedy", "exact")
    method <- pmatch(method, methods)

    rmax <- apply(tab,1,which.max)
    myseq <- 1:ncol(tab)
    cn <- colnames(tab)
    rn <- rownames(tab)
    if(is.null(cn)){
        cn <- myseq
    }
    if(is.null(rn)){
        rn <- myseq
    }

    if(method==1){
        retval <- rmax
    }
    if(method==2 | method==3){
        if(ncol(tab)!=nrow(tab)){
            stop("Unique matching only for square tables.")
        }


        dimnames(tab) <- list(myseq, myseq)
        cmax <- apply(tab,2,which.max)
        retval <- rep(NA, ncol(tab))
        names(retval) <- colnames(tab)

        baseok <- cmax[rmax]==myseq
        for(k in myseq[baseok]){
            therow <- (tab[k,])[-rmax[k]]
            thecol <- (tab[, rmax[k]])[-k]
            if(max(outer(therow, thecol, "+")) < tab[k, rmax[k]]){
                retval[k] <- rmax[k]
            }
            else{
                baseok[k] <- FALSE
            }
        }



        if(verbose){
            cat("Direct agreement:", sum(baseok),
                "of", ncol(tab), "pairs\n")
        }

        if(!all(baseok)){
            if(method==3){
                if(sum(!baseok)>maxexact){
                    method <- 2
                    warning(paste("Would need permutation of", sum(!baseok),
                                  "numbers, resetting to greedy search\n"))
                }
                else{
                    iter <- gamma(ncol(tab)-sum(baseok)+1)
                    if(verbose){
                        cat("Iterations for permutation matching:", iter, "\n")
                    }
                    perm <- permutations(ncol(tab)-sum(baseok))
                }
            }

            ## rest for permute matching
            if(any(baseok)){
                rest <- myseq[-retval[baseok]]
            }
            else{
                rest <- myseq
            }

            for(l in 1:iter){
                newretval <- retval
                if(method == 2){
                    ok <- baseok
                    while(sum(!ok)>1){
                        rest <- myseq[!ok]
                        k <- sample(rest, 1)
                        if(any(ok)){
                            rmax <- tab[k, -newretval[ok]]
                        }
                        else{
                            rmax <- tab[k,]
                        }
                        newretval[k] <- as.numeric(names(rmax)[which.max(rmax)])
                        ok[k] <- TRUE
                    }
                    newretval[!ok] <- myseq[-newretval[ok]]
                }
                else{
                    newretval[!baseok] <- rest[perm[l,]]
                }

                if(l>1){
                    agree <- sum(diag(tab[,newretval]))/sum(tab)
                    if(agree>oldagree){
                        retval <- newretval
                        oldagree <- agree
                    }
                }
                else{
                    retval <- newretval
                    agree <- oldagree <- sum(diag(tab[,newretval]))/sum(tab)
                }
            }
        }
    }

    if(verbose){
        cat("Cases in matched pairs:",
            round(100*sum(diag(tab[,retval]))/sum(tab), 2), "%\n")
    }

    if(any(as.character(myseq)!=cn)){
        retval <- cn[retval]
    }
    names(retval) <- rn

    retval
}

compareMatchedClasses <- function(x, y,
                                  method="rowmax", iter=1, maxexact=9,
                                  verbose=FALSE)
{
    if(missing(y)){
        retval <- list(diag=matrix(NA, nrow=ncol(x), ncol=ncol(x)),
                       kappa=matrix(NA, nrow=ncol(x), ncol=ncol(x)),
                       rand=matrix(NA, nrow=ncol(x), ncol=ncol(x)),
                       crand=matrix(NA, nrow=ncol(x), ncol=ncol(x)))
        for(k in 1:(ncol(x)-1)){
            for(l in (k+1):ncol(x)){
                tab <- table(x[,k], x[,l])
                m <- matchClasses(tab, method=method, iter=iter,
                                  verbose=verbose, maxexact=maxexact)
                a <- classAgreement(tab[,m])
                retval$diag[k,l] <- a$diag
                retval$kappa[k,l] <- a$kappa
                retval$rand[k,l] <- a$rand
                retval$crand[k,l] <- a$crand
            }
        }
    }
    else{
        x <- as.matrix(x)
        y <- as.matrix(y)
        retval <- list(diag=matrix(NA, nrow=ncol(x), ncol=ncol(y)),
                       kappa=matrix(NA, nrow=ncol(x), ncol=ncol(y)),
                       rand=matrix(NA, nrow=ncol(x), ncol=ncol(y)),
                       crand=matrix(NA, nrow=ncol(x), ncol=ncol(y)))
        for(k in 1:ncol(x)){
            for(l in 1:ncol(y)){
                tab <- table(x[,k], y[,l])
                m <- matchClasses(tab, method=method, iter=iter,
                                  verbose=verbose, maxexact=maxexact)
                a <- classAgreement(tab[,m])
                retval$diag[k,l] <- a$diag
                retval$kappa[k,l] <- a$kappa
                retval$rand[k,l] <- a$rand
                retval$crand[k,l] <- a$crand
            }
        }
    }

    retval
}


permutations <- function(n) {

    if(n ==1)
        return(matrix(1))
    else if(n<2)
        stop("n must be a positive integer")

    z <- matrix(1)
    for (i in 2:n) {
        x <- cbind(z, i)
        a <- c(1:i, 1:(i - 1))
        z <- matrix(0, ncol=ncol(x), nrow=i*nrow(x))
        z[1:nrow(x),] <- x
        for (j in 2:i-1) {
            z[j*nrow(x)+1:nrow(x),] <- x[, a[1:i+j]]
        }
    }
    dimnames(z) <- NULL
    z
}
