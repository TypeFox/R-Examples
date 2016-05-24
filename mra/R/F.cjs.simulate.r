F.cjs.simulate <- function( super.p, super.s, fit, N1=1000, births.per.indiv="constant.popln", R=100 ){

#   Generate R capture history matricies that follow the model in fit.

if( missing(super.p) | missing(super.s) ){
    if( missing( fit ) ){
        stop( "Either super.p and super.s must be specified, or fit must be specified")
    }

    super.p <- fit$p.hat
    super.s <- fit$s.hat[,-ncol(fit$s.hat)]

    #   Set initial capture probability, which cannot be estimated by model in fit
    #   This makes initial population size equal to number of rows in p, i.e. M_t
    super.p[,1] <- 1
} else {
    #   Make super.p and super.s matrices if they are vectors
    if(!is.matrix(super.p)) super.p <- matrix(super.p, 1)
    if(!is.matrix(super.s)) super.s <- matrix(super.s, 1)

    if( any(is.na(super.p)) | any(is.na(super.s[,-ncol(super.s)])) ){
        stop("NA's not allowed in super.p nor in all but final column of super.s")
    }
    #   super.p must have one more column than super.s, but allow them to be specified with equal numbers of columns
    #   if same size, assume last super.s column is nonsense
    if( ncol(super.p) == ncol(super.s) ){
        super.s <- super.s[,-ncol(super.s)]
    }
    if( ncol(super.p) != (ncol(super.s) + 1) ){
        stop("Number of columns in super.p must be one more than number of columns in super.s")
    }
}

if( length(births.per.indiv) != ncol(super.s) & (births.per.indiv[1] != "constant.popln") ){
    stop(paste("births.per.indiv vector must equal 'constant.popln', or length must equal", ncol(super.s)))
}

# Keeping nrows of super.p and super.s separate allows super p and super s to have different number of rows
N.super.p <- nrow(super.p)
N.super.s <- nrow(super.s)
ns <- ncol(super.p)

#   Allocate space for answer
sim.hists <- vector("list", R)

#   Main loops,  I tried to do this with apply, but...
for( k in 1:R ){

    popln <- matrix(1, N1, 1)

    #   Sample N1 from the super population for the initial realized population
    super.samp <- sample( 1:N.super.p, N1, replace=TRUE )
    p <- super.p[super.samp, ]
    hists <- rbinom( N1, 1, p[,1] )

    super.samp <- sample( 1:N.super.s, N1, replace=TRUE )
    s <- super.s[super.samp, ]

    for( j in 2:ns ){
        alive <- popln[,j-1] == 1
        survived <- captured <- rep(0,nrow(popln))
        survived[alive] <- rbinom( sum(alive), 1, s[alive,j-1] )
        captured[survived == 1] <- rbinom( sum(survived), 1, p[survived == 1,j] )
        popln <- cbind( popln, survived )
        hists <- cbind( hists, captured )

        if( births.per.indiv[1] == "constant.popln" ){
            births <- N1 - sum(survived)
        } else {
            # births.per.indiv is a vector of rates.  Multiply by current number in population
            births <- round(sum(survived) * births.per.indiv[j-1])
        }

        if( births > 0 ){
            new.popln <- new.caps <- matrix(0, births, ncol(popln))    # ncol(popln) should equal j here
            new.popln[,j] <- 1

            super.samp <- sample(1:N.super.p, births, replace=TRUE)
            new.p <- matrix(super.p[super.samp,], births, ns)   # need the matrix(...) because births can = 1, in which case new.p would be vector and only 1 dimension.
            new.caps[,j] <- rbinom( births, 1, new.p[,j] )
            popln <- rbind(popln, new.popln)
            hists <- rbind(hists, new.caps)
            p <- rbind(p, new.p)

            super.samp <- sample(1:N.super.s, births, replace=TRUE)
            s <- rbind(s, super.s[super.samp,])
        }
    }

    labs <- paste( "t", 1:ns, sep="" )
    dimnames( hists ) <- list(NULL, labs)
    dimnames( popln ) <- list(NULL, labs)

    #   Remove animals that were never captured
    hists <- hists[ rowSums(hists) > 0, ]

    #   Wrap it up
    sim.hists[[k]] <- list(hists=hists, popln.n=colSums(popln))
}

sim.hists
}
