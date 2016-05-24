linearGPmapanalysis <-
function (gmap, reference = "F2", freqmat = NULL, max.level = NULL, 
    S_full = NULL) 
{
    "kronecker_vector" <- function(matlist, kronind, row = TRUE) {
        vector <- 1
        for (i in 1:length(kronind)) {
            if (row) {
                vector <- kronecker(matlist[[i]][as.integer(kronind[i]), 
                  ], vector)
            }
            else {
                vector <- kronecker(matlist[[i]][, as.integer(kronind[i])], 
                  vector)
            }
        }
        return(vector)
    }
    gmap <- as.matrix(gmap)
    nloci <- round(log(length(gmap), 3))
    if (is.null(S_full)) {
        if (is.null(max.level) && nloci <= 8) {
            S_full = TRUE
        }
        else {
            S_full = FALSE
        }
    }
    if ((!is.null(max.level) || nloci > 8) && S_full) {
        print("Number of loci > 8 or max.level option used , setting S_full to FALSE")
        S_full = FALSE
    }
    ans <- preparelinearGPmap(gmap, reference = reference, freqmat = freqmat, 
        S_full = S_full)
    class(ans) <- "noia.linear.gpmap"
    if (S_full) {
        ans$E <- ans$sinv %*% as.matrix(ans$gmap)
        nn <- NULL
        nn <- colnames(ans$smat)
        names(ans$E) <- nn
        ans$variances <- effectsVariances(ans)
    }
    else {
        nn <- NULL
        nn <- effectsNamesGeneral(ans$nloc, max.level)
        genoindex <- expand.grid(rep(list(c(1, 2, 3)), as.integer(ceiling(ans$nloc))))
        genoindex <- genoindex[rowSums(genoindex == 1) >= (ans$nloc - 
            max.level), ]
        freq <- as.vector(ans$genofreq)
        ans$E <- NULL
        ans$variances <- NULL
        for (i in 1:dim(genoindex)[1]) {
            Sinvrow <- kronecker_vector(ans$sinv, genoindex[i, 
                ])
            Scol <- kronecker_vector(ans$smat, genoindex[i, ], 
                row = FALSE)
            ans$E[i] <- Sinvrow %*% as.matrix(ans$gmap)
            ans$variances[i] <- sum(Scol * Scol * freq) * ans$E[i]^2
        }
        names(ans$E) <- nn
        ans$variances[1] <- 0
        names(ans$variances) <- nn
    }
    ans$V_G <- sum((ans$gmap - ans$E[1])^2 * as.vector(ans$genofreq))
    return(ans)
}
