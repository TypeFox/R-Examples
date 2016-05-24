plot.FamiliasPedigree <- function (x, y, ...) 
{
    ped <- x
    id                  <- ped$id
    npers               <- length(id)
    dadid               <- rep(NA, npers)
    dadid[ped$findex>0] <- id[ped$findex[ped$findex>0]]
    momid               <- rep(NA, npers)
    momid[ped$mindex>0] <- id[ped$mindex[ped$mindex>0]]
    sex                 <- ped$sex
    n                   <- 1
    if(all(is.na(dadid)) & all(is.na(momid))) 
        stop("Cannot plot a pedigree without any relations.")
    for (i in 1:npers) {
        if (is.na(dadid[i]) && !is.na(momid[i])) {
            id        <- c(id, paste("added", n))
            dadid     <- c(dadid, NA)
            momid     <- c(momid, NA)
            sex       <- c(sex, "male")
            dadid[i]  <- id[npers + n]
            n         <- n+1
        } else if (!is.na(dadid[i]) && is.na(momid[i])) {
            id        <- c(id, paste("added", n))
            dadid     <- c(dadid, NA)
            momid     <- c(momid, NA)
            sex       <- c(sex, "female")
            momid[i]  <- id[npers + n]
            n         <- n+1
        }
    }
    plot(pedigree(id, dadid, momid, sex), ...)
}
