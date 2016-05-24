recomputeG <-
function (reg, loc) 
{
    newE <- isolateLoc(reg, loc)
    fullG <- reg$smat %*% newE
    newvcovE <- matrix(0, ncol = length(newE), nrow = length(newE))
    rownames(newvcovE) <- colnames(newvcovE) <- names(newE)
    oldvcovE <- vcov(reg$regression)
    rownames(oldvcovE) <- colnames(oldvcovE) <- names(reg$E)
    newvcovE[names(reg$E)[newE != 0], names(reg$E)[newE != 0]] <- oldvcovE[names(reg$E)[newE != 
        0], names(reg$E)[newE != 0]]
    stdG <- as.matrix(sqrt(diag(reg$smat %*% newvcovE %*% t(reg$smat))))
    gg <- as.character(genNames(length(loc)))
    gg <- matrix(unlist(strsplit(gg, "")), byrow = TRUE, nrow = length(gg))
    fullG.names <- matrix("1", nrow = nrow(gg), ncol = reg$nloc)
    fullG.names[, loc] <- gg
    fullG <- cbind(fullG[apply(fullG.names, 1, paste, collapse = ""), 
        ], stdG[apply(fullG.names, 1, paste, collapse = ""), 
        ])
    rownames(fullG) <- genNames(length(loc))
    return(fullG)
}
