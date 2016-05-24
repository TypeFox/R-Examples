# incluye estandarizacion de las marcas de los patrones simulados

syrjala.test <- 
function (ppp1, ppp2, nsim = 999) 
{
    
    datanames <- c(deparse(substitute(ppp1)), deparse(substitute(ppp2)))
    gammaf <- function(cosa.ppp) {
        gamma <- NULL
        for (k in 1:cosa.ppp$n) {
            gamma <- c(gamma, sum(cosa.ppp$marks[cosa.ppp$x <= 
                cosa.ppp$x[k] & cosa.ppp$y <= cosa.ppp$y[k]]))
        }
        return(gamma)
    }
    psi <- function(ppp1, ppp2) {
        gamma1 <- gammaf(ppp1)
        gamma2 <- gammaf(ppp2)
        psi <- sum((gamma1 - gamma2)^2)
        return(psi)
    }
    psimean <- function(ppp1, ppp2) {
        
        psi1 <- psi(ppp1, ppp2)
        psi2 <- psi(affine(ppp1, mat = diag(c(-1, -1))), affine(ppp2, 
            mat = diag(c(-1, -1))))
        psi3 <- psi(affine(ppp1, mat = diag(c(1, -1))), affine(ppp2, 
            mat = diag(c(1, -1))))
        psi4 <- psi(affine(ppp1, mat = diag(c(-1, 1))), affine(ppp2, 
            mat = diag(c(-1, 1))))
        psimean <- (psi1 + psi2 + psi3 + psi4)/4
        return(psimean)
    }
    rpsi <- function(pepesim) {
        marcas <- cbind(ppp1$marks, ppp2$marks)
        marcas <- t(apply(marcas, 1, sample))
        ppp1$marks <- marcas[, 1]/sum(marcas[, 1])
        ppp2$marks <- marcas[, 2]/sum(marcas[, 2])
        psi.sim <- psimean(ppp1, ppp2)
    }
    ppp1$marks <- ppp1$marks/sum(ppp1$marks)
    ppp2$marks <- ppp2$marks/sum(ppp2$marks)
    psi.obs <- psimean(ppp1, ppp2)
    pepesim <- list(pepesim = list(ppp1 = ppp1, ppp2 = ppp2))
    pepesim <- rep(pepesim, nsim)
    psi.sim <- sapply(pepesim, rpsi)
    result <- list(psi.obs = psi.obs, psi.sim = psi.sim, datanames = datanames, 
        nsim = nsim)
    class(result) <- c("ecespa.syrjala", class(result))
    return(result)
}
