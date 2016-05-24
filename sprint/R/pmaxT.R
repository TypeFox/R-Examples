.naNUM <- -93074815
.BLIM <- 2^30

pmaxT <- function(X, classlabel, test="t", side="abs", fixed.seed.sampling="y",
                  B=10000, na=.naNUM, nonpara="n")
{

    if(is.factor(classlabel))
        classlabel<-unclass(classlabel)-1

    if ( is.numeric(classlabel) ) {
        extra <- max(classlabel)+1
    } else {
        warning(paste("\nClasslabel needs to be a numeric vector.\n"))
        return(FALSE)
    }

    if( !checkothers(side=side, fixed.seed.sampling=fixed.seed.sampling, B=B, na=na, nonpara=nonpara) )
        return(FALSE)

    if( !is.list(tmp <- transformX(X, classlabel, test, na, nonpara)) )
        return(FALSE)

    # Returns a vector of two elements. The value of B and if we are doing complete
    # or random permutations
    # newB[1]  -->  number of permutations
    # newB[2]  -->  0 for random, 1 for complete
    newB <- getmaxB(classlabel, test, B)

    # Check that newB is computed correctly
    if( is.na(newB[1]) )
        return(FALSE)

    # If we are doing complete permutations the fixed.seed.sampling is
    # going to be always 'n', the permutations cannot be stored in memory
    # since they are too many. Better to compute them on the fly
    if( newB[2] == 1 )
      fixed.seed.sampling<-"n"

    options <- c(test, side, fixed.seed.sampling)

    res <- .C("pmaxT", as.double(tmp$X), as.integer(tmp$m), n=as.integer(tmp$n), as.integer(tmp$classlabel),
              as.double(na), t=double(tmp$m), p=double(tmp$m), adjP=double(tmp$m), as.integer(newB[1]),
              index=integer(tmp$m), as.character(options), as.integer(extra), as.integer(newB[2]), PACKAGE="sprint")

    # If the value of n is negative then it means that MPI is not initialized and
    # the function should abort and return FALSE
    if ( res$n == -1 ) {
        warning(paste("MPI is not initialized. Function is aborted.\n"))
        return(FALSE)
    } else {

        res <- cbind(index=res$index, teststat=res$t, rawp=res$p, adjp=res$adjP)
        niceres(res, X, res[,1])
    }
}


checkothers <- function(side="abs", fixed.seed.sampling="y", B=10000, na=.naNUM, nonpara="n")
{
    if((length(B)>1) || !(is.integer(as.integer(B))) || (!is.vector(B))) {
        warning(paste("B needs to be just a integer\n","your B=",B,"\n"))
        return(FALSE)
    }

    if(B<0) {
        warning(paste("The number of Permutations (B) needs to be positive\n, If you want to complete permutation, just specify B as any number greater than the maximum number of permutation\n","your B=",B))
        return(FALSE)
    }

    if((length(na)>1) || !(is.numeric(na)) || (!is.vector(na))) {
        warning(paste("The 'na' needs to be just a number\n","your na=",na,"\n"))
        return(FALSE)
    }

    if((!is.character(side)) || (!is.vector(side)) || (length(side)>1) || (!any(side==c("upper","abs","lower")))) {
        warning(paste("The 'side' needs to be a single character from c('upper','abs','lower')\n","your side=",side,"\n"))
        return(FALSE)
    }

    if((!is.character(fixed.seed.sampling)) || (!is.vector(fixed.seed.sampling)) || (length(fixed.seed.sampling)>1) || (!any(fixed.seed.sampling==c("y","n")))) {
        warning(paste("The 'fixed.seed.sampling' needs to be a single character from c('y','n')\n","your fixed.sampling=",fixed.seed.sampling,"\n"))
        return(FALSE)
    }

    if((!is.character(nonpara)) || (!is.vector(nonpara)) || (length(nonpara)>1) || (!any(nonpara==c("y","n")))) {
        warning(paste("The 'nonpara' needs to be a single character from c('y','n')\n","your nonpara=",nonpara,"\n"))
        return(FALSE)
    }

    return(TRUE)
}


transformX <- function(X, classlabel, test, na, nonpara)
{
    X <- number2na(data.matrix(X), na)
    if( !checkX(X, classlabel, test) )
        return(FALSE)
    n <- ncol(X)
    if(test == "pairt") {
        if(n%%2 == 1) {
            warning(paste("The number of columns for X must be an even number in the pair t test\n","your X=",X,"\n your classlabel=",classlabel,
                       "\n your test=",test,"\n"))
            return(FALSE)
        }
        halfn <- n%/%2;
        evendata <- X[,c(1:halfn)*2]
        odddata <- X[,c(1:halfn)*2-1]
        vecX <- (evendata-odddata)
        vecX <- data.matrix(vecX)
    }else{
        vecX <- data.matrix(X)
    }
    if(test == "wilcoxon" || nonpara == "y"){
        for(i in c(1:nrow(vecX))){
          vecX[i,] <- rank(vecX[i,])
        }
    }
    vecX <- na2number(c(vecX), na)
    if( !is.list(newL <- transformL(classlabel, test)) )
        return(FALSE)

    list(X=vecX, m=nrow(X), n=newL$n, classlabel=newL$classlabel)
}

checkX <- function(X, classlabel, test) {

    if((!is.matrix(X)) || !(is.numeric(X))) {
        warning(paste("X needs to be a matrix\n","your X = ",X,"\n"))
        return(FALSE)
    }

    if(ncol(X) != length(classlabel)) {
        warning(paste("The number of column of X needs to be the same as the length of 'classlabel'\n","your X=",X,"\n your classlabel is",classlabel,"\n"))
        return(FALSE)
    }

    if( !checkclasslabel(classlabel, test) )
        return(FALSE)

    return(TRUE)
}

transformL <- function(classlabel, test)
{
    classlabel <- as.integer(classlabel)
    if( !checkclasslabel(classlabel, test) )
        return(FALSE)
    n <- length(classlabel)
    newL <- classlabel
    if(test == "pairt"){
        if(n%%2 == 1) {
            warning(paste("The length of 'classlabel' must be an even number in the pair t\n","your classlabel=",classlabel,"\n your test=",test="\n"))
            return(FALSE)
        }
        halfn <- n%/%2;
        n <- halfn
        newL <- rep(0,n);
        for(i in c(1:n)){
            newL[i] <- classlabel[2*i]
        }
    }
    list(classlabel=newL, n=n)
}


na2number <- function(x, na) {
    y <- x
    y[is.na(y)] <- na
    y
}

number2na <- function(x, na) {
    y <- x
    y[y == na] <- NA
    y
}

niceres <- function(res, X, index){
    newres <- res
    name <- rownames(X, do.NULL=FALSE, prefix="")
    if(missing(index)) {
        rownames(newres)<-name
    }else {
        rownames(newres)<-name[index]
    }
    newres[abs(newres)>=0.9*1e20]<-NA
    data.frame(newres)
}

getmaxB<-function(classlabel, test, B, verbose=FALSE)
{
    if(B > .BLIM) {
        warning(paste("The setting of B=",B,"is too large, Please set B<",.BLIM,"\n"))
        return(c(NA, NA))
    }

    n <- length(classlabel)

    if(test == "pairt") {
        maxB <- 2^(n %/% 2)
    }

    if(any(test == c("t", "f", "wilcoxon", "t.equalvar"))) {
        k<-max(classlabel)
        maxB <- 1
        curn <- n
        for(i in c(0:k)) {
            nk <- sum(classlabel == i)
            for(j in c(1:nk)) {
                maxB <- maxB*curn/j
                curn <- curn-1
            }
        }
    }

    if(test == "blockf") {
        k <- max(classlabel)
        maxB <- 1
        for(i in c(1:(k+1))) {
            maxB <- maxB*i
        }
        maxB <- maxB^(n%/%(k+1))
    }

    # Finished the computing of maxB
    # Decide if we are doing complete *or* random permutations

    # Check if user wants complete permutations *AND* that the complete enumeration
    # is below the limit
    if((B==0) & (maxB>.BLIM)) {
        warning(paste("The complete enumeration is too big: ", maxB, ". Please set random permutation\n"))
        return(c(NA, NA))
    }

    # Check if the user requested complete permutations
    # If the user requested more than the available complete permutations then
    # just set the permutations to all available and proceed
    if((B>maxB) || (B==0)) {
        if(verbose)
            cat("We'll do complete enumerations\n")
        return(c(maxB, 1))
    }

    # Perform random permutations
    return(c(B, 0))
}

checkclasslabel <- function(classlabel, test)
{
    classlabel <- as.integer(classlabel)
    if((!is.character(test)) || (!is.vector(test)) || (length(test)>1) || (!any(test==c("t", "f", "blockf", "pairt", "wilcoxon", "t.equalvar")))) {
        warning(paste("\nYour setting of 'test' is", test, "\nThe 'test' needs to be a single character from c('t',f','blockf','pairt','wilcoxon','t.equalvar')\n"))
        return(FALSE)
    }

    if((!is.integer(as.integer(classlabel))) ||(!is.vector(classlabel))) {
        warning("\nThe 'classlabel' needs to be just a vector of integers\n")
        return(FALSE)
    }

    if(any(test==c("t", "wilcoxon", "t.equalvar"))) {
        x <- sum(classlabel == 0)
        y <- sum(classlabel == 1)
        if((x==0) || (y==0) || (x+y<length(classlabel))) {
            warning(paste("\nIn 't' test, every number in class label needs to be 0 or 1 and neither of the 0 set or 1 set can be empty set\n",
                          "The folllowing is your setting of 'classlabel'", classlabel, "\n"))
            return(FALSE)
        }
    }

    if(test == "f") {
        tab <- table(classlabel)
        tab <- tab[tab>0]
        if(length(tab) < 2) {
            warning(paste("\nIn F test, we need at least two groups\n",
                          "Your setting of 'classlabel' is", classlabel, "\n"))
            return(FALSE)
        }
        if(sum(tab)-length(tab) < 2) {
            warning(paste("\nInsufficient df for denominator of F the settings are", classlabel, "\n"))
            return(FALSE)
        }
    }

    if(test == "pairt") {
        K <- max(classlabel)
        if(K != 1) {
            warning(paste("\nIn paired t test, we only handle two groups\n","your 'classlabel' = ", classlabel,"\n"))
            return(FALSE)
        }
        if(length(classlabel)%%2 == 1) {
            warning(paste("\nThe 'classlabel' length must be an even number in the paired t\n","your 'classlabel'=",classlabel,"\n"))
            return(FALSE)
        }
        halfn <- length(classlabel) %/% 2
        for(i in c(1:halfn)) {
            cur <- classlabel[(2*i-1):(2*i)]
            if((sum(cur==0)==0) || (sum(cur==1)==0)) {
                warning(paste("\nSome errors in specifying 'classlabel' for the paired t test for the block",i,"located at","(",2*i-1,2*i,")\n",
                              "your 'classlabel'=",classlabel,"\n"))
                return(FALSE)
            }
        }
    }

    if(test == "blockf") {
        K<-max(classlabel)
        if(K < 1) {
            warning(paste("\nIn blockF test, we need at least two groups\n",
                          "your 'classlabel'=",classlabel,"\n"))
            return(FALSE)
        }
        if(length(classlabel)%%(K+1) > 0) {
            warning(paste("The classlabel length must be the multiple of the number of treatments in the block test\n","your 'classlabel'=",classlabel,"\n"))
            return(FALSE)
        }
        B <- length(classlabel) %/% (K+1)
        for(i in c(1:B)) {
            cur<-classlabel[c((K+1)*(i-1)+1):((K+1)*i)]
            # To check if cur is a permutation of c(0,1,..,K)
            for(j in c(0:K)) {
                if(sum(cur==j) == 0) {
                    warning(paste("The 'classlabel' has some errors for the blockf test at block",i,"located at",
                                  "(",(K+1)*(i-1)+1,(K+1)*i,")","There is no elements =",j,"within this block\n","your 'classlabel'=",classlabel,"\n"))
                    return(FALSE)
                }
            }
        }
    }

    return(TRUE)
}


