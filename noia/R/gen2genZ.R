gen2genZ <-
function (gen) 
{
    ans <- NULL
    n.gen <- apply(!is.na(gen), 2, "sum")
    freqs <- rbind(apply(gen == noia::genotypesNames[1], 2, "sum", 
        na.rm = TRUE)/n.gen, apply(gen == noia::genotypesNames[2], 
        2, "sum", na.rm = TRUE)/n.gen, apply(gen == noia::genotypesNames[3], 
        2, "sum", na.rm = TRUE)/n.gen)
    for (c in 1:(ncol(gen))) {
        tmp <- NULL
        for (g in gen[, c]) {
            t <- NULL
            if (is.na(g)) {
                t <- c(freqs[1, c], freqs[2, c], freqs[3, c])
            }
            else if (g == noia::genotypesNames[1]) {
                t <- c(1, 0, 0)
            }
            else if (g == noia::genotypesNames[2]) {
                t <- c(0, 1, 0)
            }
            else if (g == noia::genotypesNames[3]) {
                t <- c(0, 0, 1)
            }
            else {
                stop("Genotype ", g, " unknown.")
            }
            tmp <- rbind(tmp, t)
        }
        ans <- cbind(ans, tmp)
    }
    rownames(ans) <- rownames(gen)
    n <- NULL
    for (l in colnames(gen)) {
        n <- c(n, paste(l, noia::genotypesNames, sep = "-"))
    }
    colnames(ans) <- n
    return(ans)
}
