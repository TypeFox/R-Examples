## functions for supporting generation of large resolution V designs
## some functions are also useful for dealing with generators in general
## it would make sense to tidy up utilcat together with some of these


revdigits <- function(obj, ndigits = 1 + floor(log2(max(obj)))){
    aus <- digitsBase(obj, ndigits=ndigits)[rev(1:ndigits),]
    if (!is.matrix(aus)) aus <- matrix(aus, ncol=1)
    aus
}

indexcalc <- function(obj, k = 1 + floor(log2(max(obj)))){
    ## obj is numeric vector with positive column numbers of Yates matrix
    digmat <- revdigits(obj, ndigits=k)
    ## Brief Yates column names, because large designs otherwise not possible
    namen <- apply(digmat, 2, function(obj) paste(Letters[which(obj==1)],collapse=""))
    zahlen <- lapply(as.data.frame(digmat), function(obj) which(obj==1))
    names(zahlen) <- namen
    zahlen
    }

gencalc <- function(gen){
    ## function to return Yates column numbers (= Walsh indices)
    ## from all types of admissible generators

    ## obj can be
    ##       a numeric vector of column numbers (remains unchanged then),
    ##  OR   a character vector with factor letter combinations (from Letters;
    ##       works for up to 50 base factors only, which does not currently seem a restriction)
    ##  OR   a list with numerical vectors of the base factor indices
    ##       in the generators must be used

    if (is.numeric(gen)) aus <- gen
    else{
    ## k is the minimum number needed for current generator
    if (is.list(gen)) k <- max(unlist(gen)) ## list entry for combinations
         else k <- max(which(Letters %in% unlist(strsplit(gen,""))))
    if (is.character(gen)) {hilf <- strsplit(gen,"")
    hilf <- lapply(hilf, function(obj1) which(Letters %in% unlist(obj1)))}
    else hilf <- gen
    aus <- sapply(hilf, function(obj2) sum(2^((1:k)-1)*(1:k %in% obj2)))
    }
    aus
}

YatesFly <- function(walshindex, k=NULL){
    ## the walshindex numbers include the numbers for the base columns
    ## in line with what Sanchez and Sanchez gave in their article

    ## it is possible but usually not wise to increase the dimension
    ##    which will then loose the first k-n.base.cols non-base columns as generators

    ## function for implementing Sanchez Sanchez and other non-catalogued designs
    ## dimension if not given
    if (is.null(k)) k <- 1 + floor(log2(max(walshindex)))
    ## base columns
    base.facs <- which(log2(walshindex)%%1==0)
    if (length(base.facs) > k) stop("k must be at least ", length(base.facs))
    if (length(base.facs) < k) {
       warning("A smaller design would have been possible.\nThe first ", k-length(base.facs), " non-base generators have been omitted for completing the base set to k=", k)
       walshindex <- c(2^((1:k)-1), walshindex[-base.facs][-(1:(k-length(base.facs)))])
       ## recalculate base.facs for later use
       base.facs <- which(log2(walshindex)%%1==0)
    }

    nfactors <- length(walshindex)
    ## check dimension vs. base
    ## removed this check, allowing larger matrix for fewer factors
    ## not usually useful, but ...
        ##if (!length(base.facs)==k) stop("invalid walshindex vector")
    ## index positions
    bf <- 1:k; names(bf) <- Letters[1:k]
    zahlen <- c(as.list(bf),indexcalc(walshindex[-base.facs]))
    mat <- matrix(0, nrow=2^k, ncol=nfactors)
    for (i in 1:k) mat[,i] <- rep(c(-1,1), each=2^(i-1), times=2^(k-i))
    for (i in (k+1):nfactors) mat[,i] <- apply(mat[,zahlen[[i]]],1,prod)
    colnames(mat) <- names(zahlen)
    return(list(mat=mat, Yates.brief=zahlen, k=k, nfactors=nfactors, nruns=2^k,
           gen=walshindex[-base.facs]))
}

