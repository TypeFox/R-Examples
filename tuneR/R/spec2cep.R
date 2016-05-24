# This code is based on the Matlab implementations of PLP and Rasta
# feature calculations by Daniel P. W. Ellis of Columbia University /
# International Computer Science Institute.  For more details, see:
# http://www.ee.columbia.edu/~dpwe/resources/matlab/rastamat/

spec2cep <- function(spec, ncep=12, type=c("t2", "t1", "t3", "t4")){

    if(!(is.numeric(spec) && is.matrix(spec)))
      stop("'spec' has to be a numeric matrix")

    if(!(ncep==as.integer(ncep) && ncep > 0))
        stop("'ncep' has to be a positive integer")

    # DCT Matrix
    dctm23 <- function(spec, ncep){
        srow <- nrow(spec)
        # Orthogonal type
        return(cos( (seq(ncep)-1) * matrix(rep(seq(1, (2*srow-1), 2), each=ncep),
                    nrow=ncep)/(2 * srow) * pi) * sqrt(2/srow))
    }
    dctm2 <- function(spec, ncep){
        dctm <- dctm23(spec, ncep)
        # Make it unitary
        dctm[1,] <- dctm[1,]/sqrt(2)
        return(dctm)
    }
    dctm4 <- function(spec, ncep){
        srow <- nrow(spec)
        # Type 1 with implicit repeating of first, last bins
        dctm <- cos( (seq(ncep)-1) * matrix(rep(seq(srow), each=ncep),
                    nrow=ncep)/(srow+1) * pi) * 2
        dctm[,1] <- dctm[,1] + 1
        dctm[,srow] <- dctm[,srow] + ((-1)^(seq(ncep)-1))
        dctm <- dctm / (2*(srow+1))
        return(dctm)
    }
    dctm1 <- function(spec, ncep){      
        srow <- nrow(spec)
        dctm <- cos( (seq(ncep)-1) * matrix(rep(0:(srow-1), each=ncep),
                    nrow=ncep)/(srow-1) * pi) * 2/(2*(srow-1))
        dctm[,c(1, srow)] <- dctm[,c(1, srow)]/2
        return(dctm)
    }

    type <- match.arg(type)
    dctm <- switch(type,
                t1 = dctm1(spec=spec, ncep=ncep),
                t2 = dctm2(spec=spec, ncep=ncep),
                t3 = dctm23(spec=spec, ncep=ncep),
                t4 = dctm4(spec=spec, ncep=ncep))

    cep <- dctm %*% log(spec)

    return(list(cep=cep, dctm=dctm))
}
