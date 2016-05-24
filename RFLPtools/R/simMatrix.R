###############################################################################
## Compute similarity matrix for BLAST data
###############################################################################

## x: local blast output
## sequence.range: logical: w/o sequence range
## Min: minimum used in case of sequence.range
## Max: maximum used in case of sequence.range
simMatrix <- function(x, sequence.range = FALSE, Min, Max){
    stopifnot(is.data.frame(x))
    vars <- c("query.id", "subject.id", "alignment.length", "identity")
    if(!all(vars %in% names(x))){
        stop("Variables", vars[which(!(vars %in% names(x)))], "are not present in 'x'!")
    }
    if(!all(sid <- x[,"subject.id"] %in% x[,"query.id"])){
        warning("The following 'subject.id's do not occur as 'query.id's: ",
                paste(x[!sid,"subject.id"], collapse = ", "),
                "\nThese IDs have been removed before computing the similarity matrix!")
        x <- x[sid,]
    }
    ind <- x[,"query.id"] == x[,"subject.id"]

    x1 <- x[ind,]
    x2 <- split(x1, factor(x1[,"query.id"]))
    rm(x1)
    gc(verbose = FALSE)
    
    fun1 <- function(x){
        ind <- which.max(x[,"alignment.length"])
        x[ind,]
    }
    x3 <- do.call(rbind, lapply(x2, fun1))
    rm(x2)
    gc(verbose = FALSE)
    
    if(sequence.range){
        if(missing(Min) & missing(Max))
            stop("'Min' and 'Max' are missing. At least one has to be specified!")
        len <- x3[,"alignment.length"]
        if(!missing(Min)){
            if(Min > max(len)) 
                stop("'Min' is larger than maximum sequence length!")
        }else{
            Min <- min(len)
        }
        if(!missing(Max)){
            if(Max < min(len)) 
                stop("'Max' is smaller than maximum sequence length!")
        }else{
            Max <- max(len)
        }
        if(Min >= Max)
            stop("'Min' >= 'Max'!")
        ind.range <- (len >= Min) & (len <= Max)
        x31 <- x3[ind.range,]
        sequenzen <- x31[,"query.id"]
        x32 <- x[x[,"query.id"] %in% sequenzen,]
        x33 <- x32[x32[,"subject.id"] %in% sequenzen,]
        x4 <- rbind(x33, x31)
        rm(x31, x32)
        gc(verbose = FALSE)
    }else{
        x4 <- rbind(x[!ind,], x3)
    }

    mB <- round(x4[,"identity"]*x4[,"alignment.length"]/100, 0)
    x5 <- data.frame(x4, mB)
    rm(x4)
    gc(verbose = FALSE)    

    x6 <- split(x5, factor(x5[,"query.id"]))
    rm(x5)
    gc(verbose = FALSE)
    fun2 <- function(x, y){
        ind <- y[,"query.id"] == x[1,"query.id"]
        S <- x$mB/y[ind,"alignment.length"]
        cbind(x, S)
    }
    x7 <- do.call(rbind, lapply(x6, fun2, y = x3))
    rm(x3, x6)
    gc(verbose = FALSE)

    x8 <- split(x7, factor(x7[,"query.id"]))
    fun3 <- function(x, ids){
        y <- split(x, factor(x[,"subject.id"]))
        f <- function(x){
            ind <- which.max(x[,"S"])
            x[ind,]
        }
        res <- do.call(rbind, lapply(y, f))
        mis <- ids[!(ids %in% res[,"subject.id"])]
        if(length(mis) != 0){
            M <- matrix(0, ncol = ncol(res), nrow = length(mis))
            colnames(M) <- names(res)
            M <- as.data.frame(M)
            M[,"query.id"] <- res[1,"query.id"]
            M[,"subject.id"] <- mis
            return(rbind(res, M))
        }else{
            return(res)
        }
    }
    if(sequence.range){
        ids <- unique(x33[,"query.id"])
	rm(x33)
    }else{
        ids <- unique(x[,"query.id"])
    }
    x9 <- do.call(rbind, lapply(x8, fun3, ids = ids))
    rm(x8)
    gc(verbose = FALSE)

    x10 <- x9[order(x9[,"query.id"], x9[,"subject.id"]),]
    rm(x9)
    gc(verbose = FALSE)
    
    Sim <- matrix(x10[,"S"], ncol = length(ids), byrow = TRUE)
    Sim <- pmax(Sim, t(Sim))
    colnames(Sim) <- sort(ids)
    rownames(Sim) <- sort(ids)
    return(Sim)
}
