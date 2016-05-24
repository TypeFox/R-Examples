"blockrand" <-
function(n, num.levels=2, levels=LETTERS[seq(length=num.levels)],
                      id.prefix, stratum, block.sizes=1:4, block.prefix,
         uneq.beg=FALSE, uneq.mid=FALSE, uneq.min=0,uneq.maxit=10){

    treat <- vector(mode(levels))
    block.id <- numeric(0)
    block.size <- numeric(0)

    i <- 1
    while(length(treat) < n) {

        if ( (uneq.beg && (i == 1)) ||
                 (uneq.mid && ( length(treat) > n/2)) ) {
            n2 <- max(block.sizes)*length(levels) + 1
            ii <- uneq.maxit
            while( ii > 0 ) {
                tmp <- sample(levels, n2, replace=TRUE)
                if ( diff(range(table(tmp))) >= uneq.min ) { ii <- 0 }
                ii <- ii -1
            }
            treat <- c(treat, tmp)
            block.id <- c( block.id, rep(i,n2) )
            block.size <- c(block.size, rep(n2,n2))

            i <- i + 1
            if( length(treat) > n/2 ) { uneq.mid <- FALSE }
        }

        block.n <- sample(block.sizes,1)
        tmp <- rep(levels,block.n)
        n2 <- length(tmp)

        treat <- c(treat, sample(tmp))
        block.id <- c(block.id, rep(i, n2) )
        block.size <- c(block.size, rep(n2,n2))

        i <- i + 1
    }

    n3 <- length(treat)
    if(missing(id.prefix)){
        id <- seq(length.out=n3)
    } else if( is.numeric(id.prefix) ){
        id <- seq(id.prefix, length.out=n3)
    } else {
        w <- floor(logb(n3,10))+1
        id <- paste(id.prefix, formatC(seq(length.out=n3),flag='0',width=w),
                    sep='')
    }

    out <- data.frame(id=I(id))

    if(!missing(stratum)){
        out$stratum <- rep(stratum,n3)
    }

    if(missing(block.prefix)) {
        out$block.id <- factor(block.id)
    } else if (is.numeric(block.prefix)) {
        out$block.id <- factor( block.id+block.prefix )
    } else {
        w <- floor(logb(max(block.id),10))+1
        out$block.id <- factor( paste(block.prefix,
                                      formatC(block.id,flag='0',width=w),
                                      sep='') )
    }

    out$block.size <- block.size

    out$treatment <- factor(treat)

    return(out)
}

