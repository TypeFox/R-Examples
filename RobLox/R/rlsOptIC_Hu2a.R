###############################################################################
## psi function
###############################################################################
.Hu2arlsGetpsi <- function(x, k1, k2){
    return(sign(x)*pmax(k1,pmin(abs(x),k2)))
}

###############################################################################
## computation of bias
###############################################################################
.Hu2arlsGetbias <- function(x, k1, k2){
    beta.k12 <- (k1^2*(2*pnorm(k1) - 1) - 2*(k2*dnorm(k2)-k1*dnorm(k1)) 
                    + 2*(pnorm(k2)-pnorm(k1)) + 2*k2^2*pnorm(-k2))
    return(sqrt(.Hu2arlsGetpsi(x = x, k1 = k1, k2 = k2)^2/(2*(pnorm(k2)-pnorm(k1)))^2 
                + ((.Hu2arlsGetpsi(x = x, k1 = k1, k2 = k2)^2 - beta.k12)^2/
                  (4*k1*dnorm(k1)-4*k2*dnorm(k2)+4*(pnorm(k2)-pnorm(k1)))^2)))
}


###############################################################################
## computation of asymptotic variance
###############################################################################
.Hu2arlsGetvar <- function(k1, k2){
    beta.k12 <- (k1^2*(2*pnorm(k1) - 1) - 2*(k2*dnorm(k2) - k1*dnorm(k1)) 
                    + 2*(pnorm(k2)-pnorm(k1)) + 2*k2^2*pnorm(-k2))
    Var.loc <- beta.k12/(2*(pnorm(k2)-pnorm(k1)))^2

    h1 <- (k1^4*(2*pnorm(k1)-1) - 2*(k2^3*dnorm(k2)-k1^3*dnorm(k1))
           - 6*(k2*dnorm(k2)-k1*dnorm(k1)) + 6*(pnorm(k2)-pnorm(k1)) 
           + 2*k2^4*pnorm(-k2))
    Var.sc <- (h1 - beta.k12^2)/(4*k1*dnorm(k1)-4*k2*dnorm(k2)+4*(pnorm(k2)-pnorm(k1)))^2

    return(Var.loc + Var.sc)
}

###############################################################################
## computation of maximum asymptotic MSE
###############################################################################
.Hu2arlsGetmse <- function(k12, r, MAX){
    k1 <- k12[1]; k2 <- k12[2]

    # constraints
    if(k1 < 0 || k2 < 0 || k2 <= k1) return(MAX)

    beta.k12 <- (k1^2*(2*pnorm(k1) - 1) - 2*(k2*dnorm(k2) - k1*dnorm(k1)) 
                    + 2*(pnorm(k2)-pnorm(k1)) + 2*k2^2*pnorm(-k2))
    A.loc <- 1/(2*(pnorm(k2)-pnorm(k1)))
    A.sc <- 1/(4*k1*dnorm(k1)-4*k2*dnorm(k2)+4*(pnorm(k2)-pnorm(k1)))

    bias <- max(.Hu2arlsGetbias(x = 0, k1 = k1, k2 = k2), 
                .Hu2arlsGetbias(x = k1, k1 = k1, k2 = k2), 
                .Hu2arlsGetbias(x = k2, k1 = k1, k2 = k2), 
                .Hu2arlsGetbias(x = sqrt(beta.k12), k1 = k1, k2 = k2), 
                .Hu2arlsGetbias(x = sqrt(max(0,beta.k12-0.5*A.loc^2/A.sc^2)), 
                                k1 = k1, k2 = k2))

    Var <- .Hu2arlsGetvar(k1 = k1, k2 = k2)

    return(Var + r^2*bias^2)
}

###############################################################################
## optimal IC
###############################################################################
rlsOptIC.Hu2a <- function(r, k1.start = 0.25, k2.start = 2.5, delta = 1e-6, MAX = 100){
    res <- optim(c(k1.start, k2.start), .Hu2arlsGetmse, method = "Nelder-Mead", 
                control = list(reltol=delta), r = r, MAX = MAX)

    k1 <- res$par[1]; k2 <- res$par[2]

    beta.k12 <- (k1^2*(2*pnorm(k1) - 1) - 2*(k2*dnorm(k2) - k1*dnorm(k1)) 
                    + 2*(pnorm(k2)-pnorm(k1)) + 2*k2^2*pnorm(-k2))
    A.loc <- 1/(2*(pnorm(k2)-pnorm(k1)))
    A.sc <- 1/(4*k1*dnorm(k1)-4*k2*dnorm(k2)+4*(pnorm(k2)-pnorm(k1)))

    bias <- max(.Hu2arlsGetbias(x = 0, k1 = k1, k2 = k2), 
                .Hu2arlsGetbias(x = k1, k1 = k1, k2 = k2), 
                .Hu2arlsGetbias(x = k2, k1 = k1, k2 = k2), 
                .Hu2arlsGetbias(x = sqrt(beta.k12), k1 = k1, k2 = k2), 
                .Hu2arlsGetbias(x = sqrt(max(0,beta.k12-0.5*A.loc^2/A.sc^2)), 
                                k1 = k1, k2 = k2))

    fct1 <- function(x){ A.loc*sign(x)*pmax(k1,pmin(abs(x),k2)) }
    body(fct1) <- substitute({ A.loc*sign(x)*pmax(k1,pmin(abs(x),k2)) },
                        list(k1 = k1, k2 = k2, A.loc = A.loc))
    fct2 <- function(x){ A.sc*(pmax(k1,pmin(abs(x),k2))^2 - beta.k12) }
    body(fct2) <- substitute({ A.sc*(pmax(k1,pmin(abs(x),k2))^2 - beta.k12) },
                        list(k1 = k1, k2 = k2, beta.k12 = beta.k12, A.sc = A.sc))

    return(IC(name = "IC of Hu2a type", 
              Curve = EuclRandVarList(RealRandVariable(Map = list(fct1, fct2), Domain = Reals())),
              Risks = list(asMSE = res$value, asBias = bias, asCov = res$value - r^2*bias^2), 
              Infos = matrix(c("rlsOptIC.Hu2a", "optimally robust IC for Hu2a estimators and 'asMSE'",
                               "rlsOptIC.Hu2a", paste("where k1 =", round(k1, 3), "and k2 =", round(k2, 3))), 
                           ncol=2, byrow = TRUE, dimnames=list(character(0), c("method", "message"))), 
              CallL2Fam = call("NormLocationScaleFamily")))
}
