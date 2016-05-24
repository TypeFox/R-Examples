null_basis <- function(object){

    .null_basis <- function(object){
        S <- svd( object )
        id <- S$d<1e-15
        if (any(id)){
            null.basis <- S$v[,id, drop=FALSE]
            null.basis
        }
        else {
            matrix(nrow=ncol(object), ncol=0)
        }
    }
    
    if (is.character(object)){
        stop("'object' of type 'character' not valid")
    }
    
    if (is.vector(object)){
        object <- matrix(object, ncol=1)
    }

    if (class(object) %in% c("matrix","Matrix")){
        .null_basis(object)
    } else {
        m <- try(model.matrix(object), silent=TRUE)
        if (class(m) != "try-error")
            .null_basis(m)
        else
            stop("Can not find null basis for 'object'")
    }
}


is_estimable <- function(K, null.basis){
    if (is.null(null.basis) ||
        (is.matrix(null.basis) && ncol(null.basis)==1) ){
        rep(TRUE, nrow(K))
    } else {
        out <- lapply(1:nrow(K),
                      function(i){
                          k <- K[i,]
                          all(abs(apply(null.basis, 2, function(x) sum(k * x))) < 1e-04)
                      })
        unlist( out )
    }
}



#' ## FIXME: do_contrast is never used
#' .do_contrast <- function(KK, bhat, V0, ddf, is.est){
#'     used       <- which(!is.na(bhat))
#'     not.used   <- which(is.na(bhat))
#'     bhat.used  <- bhat[used]
#'     VV <- V0
#'     ddfm <- function(kk, se) ddf
#'     res <- matrix(NA, nrow=nrow(KK), ncol=3)
#'     for (ii in 1:nrow(res)){
#'         if (is.est[ii]){
#'             kk   <- KK[ii,used]
#'             est  <- sum(kk*bhat.used)
#'             se   <- sqrt(sum(kk * (VV %*% kk)))
#'             df2  <- ddfm(kk, se)
#'             res[ii,] <- c(est, se, df2)
#'         }
#'     }
#'     colnames(res) <- c("estimate","SE","df")
#'     res
#' }

#' ## FIXME: do_pairs is never used
#' .do_pairs <- function(KK){ ## pairwise differences of lsmeans
#'     NN <- nrow(KK)
#'     MM <- ncol(KK)
#'     SS <- as.matrix(attr(KK,"grid")) ## todo: check if null
#'     EE <- vector("list", NN*(NN-1)/2)
#'     DD <- matrix(0, nrow=NN*(NN-1)/2, ncol=NN)
#'     kk <- 1
#'     for (ii in 1:(NN-1)){
#'         for (jj in (ii+1):NN){
#'             DD[kk, ii] <- 1
#'             DD[kk, jj] <- -1
#'             EE[[kk]] <- list(ff=SS[ii,], .ff=SS[jj,])
#'             kk <- kk + 1
#'         }
#'     }
#'     ff <- do.call(rbind, lapply(EE, function(x) x$ff))
#'     .ff <- do.call(rbind, lapply(EE, function(x) x$.ff))
#'     colnames(.ff) <- paste(".", colnames(.ff),sep='')
#'     pair <- as.data.frame(cbind(ff, .ff))
#'     list(DD=DD, grid=pair)
#' }
