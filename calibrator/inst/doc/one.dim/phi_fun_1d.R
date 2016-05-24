# This file is intented to be called by calex_1d.R.  It specifies a
# phi.fun for use in the 1-d problem.>


 phi.fun.1d <- 
function (rho, lambda, psi1, psi1.apriori, psi2, psi2.apriori, 
    theta.apriori) 
{

      pdm.maker.psi1 <- function(psi1) {
        jj.omega_x <- diag(psi1[1],nrow=1)
        rownames(jj.omega_x) <- "x"
        colnames(jj.omega_x) <- "x"
        jj.omega_t <- diag(psi1[2],ncol=1)
        rownames(jj.omega_t) <- "A"
        colnames(jj.omega_t) <- "A"
        sigma1squared <- psi1[3]
        return(list(omega_x = jj.omega_x, omega_t = jj.omega_t, 
            sigma1squared = sigma1squared))
    }
    pdm.maker.psi2 <- function(psi2) {
        jj.omegastar_x <- diag(psi2[1],ncol=1)
        rownames(jj.omegastar_x) <- "x"
        colnames(jj.omegastar_x) <- "x"
        sigma2squared <- psi2[2]
        return(list(omegastar_x = jj.omegastar_x, sigma2squared = sigma2squared))
    }
    # CHANGES END: remainder of function unaltered
    jj.mean <- theta.apriori$mean
    jj.V_theta <- theta.apriori$sigma
    jj.discard.psi1 <- pdm.maker.psi1(psi1)
    jj.omega_t <- jj.discard.psi1$omega_t
    jj.omega_x <- jj.discard.psi1$omega_x
    jj.sigma1squared <- jj.discard.psi1$sigma1squared
    jj.discard.psi2 <- pdm.maker.psi2(psi2)
    jj.omegastar_x <- jj.discard.psi2$omegastar_x
    jj.sigma2squared <- jj.discard.psi2$sigma2squared
    jj.omega_t.upper <- chol(jj.omega_t)
    jj.omega_t.lower <- t(jj.omega_t.upper)
    jj.omega_x.upper <- chol(jj.omega_x)
    jj.omega_x.lower <- t(jj.omega_x.upper)
    jj.a <- solve(solve(jj.V_theta) + 2 * jj.omega_t, solve(jj.V_theta, 
        jj.mean))
    jj.b <- t(2 * solve(solve(jj.V_theta) + 2 * jj.omega_t) %*% 
        jj.omega_t)
    jj.c <- jj.sigma1squared/sqrt(det(diag(nrow = nrow(jj.V_theta)) + 
        2 * jj.V_theta %*% jj.omega_t))
    names(jj.c) <- "ht.fun.precalc"
    jj.A <- solve(jj.V_theta + solve(jj.omega_t)/4)
    jj.A.upper <- chol(jj.A)
    jj.A.lower <- t(jj.A.upper)
    list(rho = rho, lambda = lambda, psi1 = psi1, psi1.apriori = psi1.apriori, 
        psi2 = psi2, psi2.apriori = psi2.apriori, theta.apriori = theta.apriori, 
        omega_x = jj.omega_x, omega_t = jj.omega_t, 
        omegastar_x = jj.omegastar_x, sigma1squared = jj.sigma1squared, 
        sigma2squared = jj.sigma2squared, omega_x.upper = jj.omega_x.upper, 
        omega_x.lower = jj.omega_x.lower, omega_t.upper = jj.omega_t.upper, 
        omega_t.lower = jj.omega_t.lower, a = jj.a, b = jj.b, 
        c = jj.c, A = jj.A, A.upper = jj.A.upper, A.lower = jj.A.lower)
}
