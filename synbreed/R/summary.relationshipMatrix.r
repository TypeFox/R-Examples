summary.relationshipMatrix <- function(object,...){
     relMat <- object
     offdiag <- relMat[lower.tri(relMat,diag=FALSE)]
     ans <- list(dim=c(nrow=nrow(relMat),ncol=ncol(relMat)),
                 rank=try(qr(relMat)$rank, silent=TRUE),
                 range.off.diagonal=c(min=min(offdiag, na.rm=TRUE), max=max(offdiag, na.rm=TRUE)),
                 mean.diag=mean(diag(relMat), na.rm=TRUE),
                 mean.off.diag=mean(offdiag, na.rm=TRUE),
                 nUnique=length(unique(as.vector(offdiag))),
                 diag.val=summary(as.vector(diag(relMat))),
                 empty = sum(is.na(relMat)))
     class(ans) <- "summary.relationshipMatrix"
     ans
}

print.summary.relationshipMatrix <- function(x,...){
    if(class(x$rank) == "try-error") warning("\n\n  There are ", x$empty, " of ", prod(x$dim), " values missing in your relationshipMatrix!\n  ",
                                             "Computation is done with removed 'NA' values.\n", immediate.=TRUE)
    cat(" dimension                    ",x$dim[1],"x",x$dim[2],"\n")
    if(class(x$rank) == "try-error"){
        cat(" rank                          not computable because of the missing values!\n")
    } else {
        cat(" rank                         ",x$rank,"\n")
    }
    cat(" range of off-diagonal values ",x$range.off.diagonal[1],"--",x$range.off.diagonal[2],"\n")
    cat(" mean off-diagonal values     ",x$mean.off.diag,"\n")
    cat(" range of diagonal values     ",x$diag.val[1],"--",x$diag.val[6],"\n")
    cat(" mean diagonal values         ",x$mean.diag,"\n")
    cat(" number of unique values      ",x$nUnique,"\n")
}
