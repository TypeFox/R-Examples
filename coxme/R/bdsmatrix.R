# Automatically generated from all.nw using noweb
#Functions for moving back and forth between Matrix and bdsmatrix objects
rowTocol <- function(bs) {  #bs = size of block
    n <- (bs*(bs+1))/2
    indx <- integer(n)
    offset <- c(0L, cumsum(seq.int(bs-1, 1)))
    k <- 1L
    for (i in seq.int(1,bs)) {
        for (j in 1:i ) {
            indx[k] <- i + offset[j]
            k <- k+1
            }
        }
    indx                              
}

colTorow <- function(bs) {  #bs = size of block
    n <- (bs*(bs+1))/2
    indx <- integer(n)
    offset <- c(0L, cumsum(seq.int(1, bs-1)))
    k <- 1L
    for (i in seq.int(1,bs)) {
        for (j in seq.int(i, bs) ) {
            indx[k] <- i + offset[j]
            k <- k+1
            }
        }
    indx                              
}
setAs("bdsmatrix", "dsCMatrix", function(from) {
    if (length(from@blocks)==0) symmpart(Matrix(from@rmat))
    else {
        temp <- .Call("bds_dsc",
                      from@blocksize,
                      from@blocks,
                      from@rmat,
                      from@Dim)
         new("dsCMatrix", 
                 i = temp$i,
                 p = temp$p,
                 x=  temp$x,
                 Dim = from@Dim,
                 Dimnames= lapply(from@Dimnames, as.character),
                 uplo='U',
                 factors=list())
    }
})             
setAs("gchol.bdsmatrix", "dtCMatrix", function(from) {
    dd <- sqrt(diag(from))  #the multiplication factor
    rownum <- function(z) unlist(lapply(1:z, function(r) 1:r))
    nb <- from@blocksize* (from@blocksize+1)/2  #elements per block
    if (length(from@blocks)>0){
        temp <- vector('list', length(nb))
        bstart <- c(0, cumsum(from@blocksize)) #offset of each block
        for (i in 1:length(nb)) 
            temp[[i]] <- rownum(from@blocksize[i]) + bstart[i]
        m.i <- unlist(temp)
        m.p <- unlist(lapply(from@blocksize, function(x) seq(1,x)))

        xindx <- unlist(sapply(from@blocksize, rowTocol)) +
                 rep.int(c(0, cumsum(nb))[1:length(nb)], nb)      
        m.x <- from@blocks[xindx]

        if (length(from@rmat >0)) {
            nc <- ncol(from@rmat)  #number of columns in rmat
            nr <- nrow(from)     #total number of rows
            ii <- seq(to=nr, length=nc)
            m.i <- c(m.i, unlist(lapply(ii, function(r) 1:r)))
            m.p <- c(m.p, ii)
            m.x <- c(m.x, from@rmat[row(from@rmat) <= nr + col(from@rmat) -nc])
        }
    }
    else {
        nc <- ncol(from@rmat)  #number of columns in rmat
        nr <- nrow(from)     #total number of rows
        ii <- seq(to=nr, length=nc)
        m.i <- unlist(lapply(ii, function(r) 1:r))
        m.p <- ii
        m.x <- from@rmat[row(from@rmat) <= nr + col(from@rmat) -nc]
    } 
        
    #Modify x
    m.x <- m.x * dd[m.i]  # fixes the off diagonals
    m.x[rep(1:length(m.p), m.p) ==m.i] <- dd  #diagonals
    new("dtCMatrix", 
        i = as.integer(m.i-1),
        p = as.integer(c(0, cumsum(m.p))),
        Dim= dim(from),
        Dimnames=from@Dimnames,
        x= m.x,
        uplo='U',
        diag='N')
})             
# Code to find the subset myself
#  This ONLY works for the special case below
getblock <- function(x, start, end) {
    nrow <-  as.integer(1+end-start)
    xp <- x@p[start:(end+1)]
    if (nrow==1) return(x@x[xp[1]+1])  #singleton element

    keep <- (1+min(xp)):max(xp)   
    new("dsCMatrix", i=x@i[keep]+ 1L - as.integer(start), 
        p= xp- min(xp), 
        Dim=c(nrow, nrow), Dimnames=list(NULL, NULL),
        x = x@x[keep], uplo=x@uplo, factors=list())
}

setAs("dsCMatrix", "bdsmatrix", function(from) {
    dd <- dim(from)
    if (dd[1] != dd[2]) stop("Variance matrices must be square")

    nc <- ncol(from)
    minrow <- from@i[from@p[1:nc] +1] +1
    minrow <- rev(cummin(rev(minrow)))
    block.start <- which(1:nc == minrow)
    block.end <- c(block.start[-1] -1, ncol(from))
    nblock <- length(block.start)
                     
    blocks <- vector('list', nblock)  
    for (i in 1:nblock) {
#        indx <- block.start[i]:block.end[i]
#        blocks[[i]] <- as.matrix(from[indx, indx])
        blocks[[i]] <- as.matrix(getblock(from, block.start[i], block.end[i]))
    }
    bdsmatrix(blocksize=sapply(blocks, nrow), blocks=unlist(blocks), 
              dimnames=dimnames(from))
    })
