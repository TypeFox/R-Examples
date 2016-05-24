#' @title  Model matrix for point-process responses
#'
#' @description
#' Constructs  a \code{data.frame} to be fitted using \code{\link{pro}}.
#' Reference: X Luo, S Gee, V Sohal, D Small (In Press). A Point-process Response Model for Optogenetics Experiments on Neural Circuits. _Statistics in Medicine_.
#'
#' @export
#' 
#' @param spike A binary vector represents spiking (1) or no spiking (0).
#' @param flash A binary vector of the same length of \code{spike}, 1 for flashing and 0 for non-flashing.
#' @param fixed Whether a fixed time window of spike/flash history should be used. If  it is \code{NULL}, a varying time window of history will be used as described in the reference. If it is a integer \code{j}, a fixed window from index \code{t-j} to \code{t} will be used. 
#' @param kv Whether the history dependence model in Kass and Ventura (2001) (A Spike-Train Probability Model, Neural Computation 13, 1713-1720) should be employed. This differs from the history dependence model in the reference.
#' @return  a \code{data.frame} of the three response functions (PF, CF, SF) and other intermediate functions (for future modeling use).
#' @examples
#' n <- 500
#' set.seed(100)
#' re <- sim.lif(n, rbinom(n, 1, 0.14), 7, 3)
#' d <- model.pro(re$sbin, re$I)
#' d[1:10, ]
model.pro <- function(spike, flash, fixed=NULL, kv=F) {
    ## Internal functions
    findInd <- function(sn, tf, first=T) {
        tmp <- sn[tf]
        if (sum(tf) == 0) {
            return(NA)
        } else {
            if (first) {
                return(tmp[1])
            } else {
                return(tmp[length(tmp)])
            }
        }
    }
    

    fixedflashspacing <- function(ff) {
        ## ff is a vector of legnth fixed, 0 or 1s
        re <- 0
        sff <- sum(ff)
        fseq <- 1:length(ff) 
        if (sff == 0 ) {
            return(length(ff)^2)
        } else if (sff == 1) {
            af <- fseq[ff==1]
            return(af^2 + (length(ff)-af)^2)
        } else {
            fftimes <- fseq[ff==1]
            return(sum(diff(fftimes)^2) + fftimes[1]^2 + (length(ff)-fftimes[length(fftimes)])^2)
        }
    }

    
    nsf <- length(spike)
    t0 <- 1:nsf
    t0 <- t0[spike==1]
    if (length(t0) < 2 ) {
        warning("only one spike in sweep!")
        return(NULL)
    }

    if (is.null(fixed)) {
        t0 <- t0[1]                           #burn-in time
    } else {
        t0 <- fixed
    }
    
    ## find past flash
    tmp <- flash[1:t0]
    if (sum(tmp)> 0) {
        ## if there is at least on eflash before current spike [1, t0] 
        f0 <- findInd(seq_along(tmp), tmp==1, first=F)
    } else {
        ## if no flash, set it to 1 
        f0 <- 1 
    }
    
    nout <- nsf-t0
    re <- matrix(0, nrow=nout, ncol=7)
    curspike <- t0
    curflash <- f0
    curflashspike <- t0

    cumflash <- 0
    cumpftime <- 0
    cumsqpf <- 0
    sqcount <- 0
    curcount <- 0
    zeroout <- 1

      
    if (!is.null(fixed)) {
        fixedcumspike <- sum(spike[1:t0])
        fixedcumflash <- sum(flash[1:t0])
        fixedcumsqcount <- fixedflashspacing(flash[1:t0])
    }
    

    for (j in (t0+1):nsf) {
        if(zeroout) {
            cumflash <- 0
            cumpftime <- 0
            zeroout <- 0
            sqcount <- 0
            curcount <- 0
        }
        pstime <- j-curspike                #before revealing whether it is a spike or not
        if (spike[j] ==1) {
            curspike <- j
            ## should not reset. after is not significant
            if (kv) curflashspike <- j                     #20110810: reset pf to zero after spike
            zeroout <- 1
        } 
        if (flash[j] ==1) {
            ##  current flash before previous spike
            curflash <- j
            curflashspike <- j 
            cumflash <- cumflash+1
        }
        cumpftime <- cumpftime+j-curflash
        if (j-curflash == 0) sqcount <- sqcount+curcount
        if (kv) {
            curcount <- (j-curflashspike)^2
        } else {
            curcount <- (j-curflash)^2
        }
        cumsqpf <- sqcount+curcount

        ## if fixed
        if (is.null(fixed)) {
            re[j-t0,] <- c(spike[j], flash[j], pstime, j-curflash, cumflash, cumpftime, cumsqpf)
        } else {
            re[j-t0,] <- c(spike[j], flash[j], pstime, j-curflash, sum(flash[(j-fixed):j]), cumpftime, fixedflashspacing(flash[(j-fixed):j])
                           )
        }
    }

    re <- as.data.frame(re)
    
    names(re) <- c("spike", "flash", "pstime", "pftime", "cumflash", "cumpftime", "cumsqpf")
    
    re <- data.frame(re, logpf=log(re$pftime+1), logcumf=log(re$cumflash+1), logcumpf=log(re$cumpftime+1), logcumsqpf=log(re$cumsqpf+1))

    re$loglogcumsqpf <- log(1 +re$logcumsqpf)
    
    return(re)
}

