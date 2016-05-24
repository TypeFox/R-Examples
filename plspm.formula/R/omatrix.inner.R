omatrix.inner <-
function (latents, latlist, mat = TRUE, iplot = TRUE) 
{
    N <- length(latents)
    vldroite <- unique(unlist(latlist[[2]]))
    vlgauche <- latlist[[1]]
    calc.nfois.exo <- function(vlat) {
        nbfois.exo <- function(vtot) {
            return(sum(vlat %in% vtot))
        }
        return(sum(sapply(latlist[[2]], nbfois.exo)))
    }
    calc.nfois.exo <- Vectorize(calc.nfois.exo)
    vlexo <- latents[1 - as.numeric(latents %in% vldroite)]
    calc.equ.exo <- function(vlat) {
        indx <- which(latlist[[1]] == vlat)
        if (length(indx) < 1) {
            res <- 0
        }
        else {
            res <- sum(as.numeric(vlexo %in% latlist[[2]][[indx]]))
        }
        return(res)
    }
    calc.equ.exo <- Vectorize(calc.equ.exo)
    ntotF <- sapply(latlist[[2]], length)
    calc.nbequ <- function(vlat) {
        indx <- which(latlist[[1]] == vlat)
        if (length(indx) < 1) {
            res <- 0
        }
        else {
            res <- ntotF[indx]
        }
        return(res)
    }
    calc.nbequ <- Vectorize(calc.nbequ)
    thetaj <- 1 - as.numeric(latents %in% vldroite)
    Ej <- as.vector(calc.nbequ(latents))
    gammaj <- 1 - as.numeric(latents %in% vlgauche)
    Fj <- as.vector(calc.nfois.exo(latents))
    Kj <- as.vector(calc.equ.exo(latents))
    muj <- 10^7*Ej*thetaj+(10^5*Ej-10*Fj+Kj)*(1-thetaj-gammaj)-(10^3*Fj-Kj)*gammaj
    olatents <- latents[order(muj)]
    reslist <- list(mu = muj, ordre = olatents)
    if (mat) {
        matlist <- function(vect) {
            return(as.numeric(olatents %in% vect))
        }
        Mlist <- lapply(latlist[[2]], matlist)
        mat.vect <- function(j) {
            indj <- which(latlist[[1]] == olatents[j])
            if (length(indj) < 1) {
                return(rep(0, N))
            }
            else {
                return(unlist(Mlist[indj]))
            }
        }
        mat.vect <- Vectorize(mat.vect)
        Mat <- t(mat.vect(1:N))
        rownames(Mat) <- olatents
        reslist <- c(reslist, list(matrice = Mat))
    }
    if (iplot) {
        innerplot(Mat)
    }
    return(reslist)
}
