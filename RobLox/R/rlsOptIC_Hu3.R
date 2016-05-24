###############################################################################
## chi function
###############################################################################
.Hu3rlsGetchi <- function(x, c1, c2){
    beta.c12 <- (c1^2*(2*pnorm(c1) - 1) - 2*(c2*dnorm(c2)-c1*dnorm(c1)) 
                 + 2*(pnorm(c2)-pnorm(c1)) + 2*c2^2*pnorm(-c2))
    Ind1 <- (abs(x)<c1); Ind2 <- (abs(x)>c2)
    return(c1^2*Ind1 + x^2*(1-Ind1)*(1-Ind2) + c2^2*Ind2 - beta.c12)
}

###############################################################################
## Computation of bias
###############################################################################
.Hu3rlsGetbias <- function(x, k, c1, c2){
    return(sqrt(pmin(k^2, x^2)/(2*pnorm(k)-1)^2 + .Hu3rlsGetchi(x = x,c1 = c1, c2 = c2)^2/
            (4*c1*dnorm(c1)-4*c2*dnorm(c2)+4*(pnorm(c2)-pnorm(c1)))^2))
}


###############################################################################
## computation of asymptotic variance
###############################################################################
.Hu3rlsGetvar <- function(k, c1, c2){
    beta.c12 <- (c1^2*(2*pnorm(c1) - 1) - 2*(c2*dnorm(c2) - c1*dnorm(c1)) 
                 + 2*(pnorm(c2)-pnorm(c1)) + 2*c2^2*pnorm(-c2))
    Var.loc <- (2*pnorm(k) - 1 - 2*k*dnorm(k) + 2*k^2*(1-pnorm(k)))/(2*pnorm(k)-1)^2

    hilf <- (c1^4*(2*pnorm(c1)-1) - 2*(c2^3*dnorm(c2)-c1^3*dnorm(c1)) 
             - 6*(c2*dnorm(c2)-c1*dnorm(c1)) + 6*(pnorm(c2)-pnorm(c1)) 
             + 2*c2^4*pnorm(-c2))
    Var.sc <- (hilf - beta.c12^2)/(4*c1*dnorm(c1)-4*c2*dnorm(c2)+4*(pnorm(c2)-pnorm(c1)))^2

    return(Var.loc + Var.sc)
}

###############################################################################
## computation of maximum asymptotic MSE
###############################################################################
.Hu3rlsGetmse <- function(kc12, r, MAX){
    k <- kc12[1]; c1 <- kc12[2]; c2 <- kc12[3]

    #constraints
    if(k < 0 || c1 < 0 || c2 < 0 || c2 <= c1) return(MAX)

    beta.c12 <- (c1^2*(2*pnorm(c1) - 1) - 2*(c2*dnorm(c2) - c1*dnorm(c1)) 
                 + 2*(pnorm(c2)-pnorm(c1)) + 2*c2^2*pnorm(-c2))
    A.loc <- 1/(2*pnorm(k)-1)
    A.sc <- 1/(4*c1*dnorm(c1)-4*c2*dnorm(c2)+4*(pnorm(c2)-pnorm(c1)))

    bias <- max(.Hu3rlsGetbias(x = 0, k = k, c1 = c1, c2 = c2), 
                .Hu3rlsGetbias(x = k, k = k, c1 = c1, c2 = c2), 
                .Hu3rlsGetbias(x = c1, k = k, c1 = c1, c2 = c2), 
                .Hu3rlsGetbias(x = c2, k = k, c1 = c1, c2 = c2), 
                .Hu3rlsGetbias(x = sqrt(beta.c12), k = k, c1 = c1, c2 = c2), 
                .Hu3rlsGetbias(x = sqrt(max(0, beta.c12-0.5*A.loc^2/A.sc^2)), 
                                k = k, c1 = c1, c2 = c2))

    Var <- .Hu3rlsGetvar(k = k, c1 = c1, c2 = c2)

    return(Var + r^2*bias^2)
}


###############################################################################
## optimal IC
###############################################################################
rlsOptIC.Hu3 <- function(r, k.start = 1.0, c1.start = 0.1, c2.start = 0.5, 
                        delta = 1e-6, MAX = 100){
    res <- optim(c(k.start, c1.start, c2.start), .Hu3rlsGetmse, method = "Nelder-Mead", 
                control = list(reltol=delta), r = r, MAX = MAX)

    k <- res$par[1]; c1 <- res$par[2]; c2 = res$par[3]

    beta.c12 <- (c1^2*(2*pnorm(c1) - 1) - 2*(c2*dnorm(c2) - c1*dnorm(c1)) 
                 + 2*(pnorm(c2)-pnorm(c1)) + 2*c2^2*pnorm(-c2))
    A.loc <- 1/(2*pnorm(k)-1)
    A.sc <- 1/(4*c1*dnorm(c1)-4*c2*dnorm(c2)+4*(pnorm(c2)-pnorm(c1)))

    bias <- max(.Hu3rlsGetbias(x = 0, k = k, c1 = c1, c2 = c2), 
                .Hu3rlsGetbias(x = k, k = k, c1 = c1, c2 = c2), 
                .Hu3rlsGetbias(x = c1, k = k, c1 = c1, c2 = c2), 
                .Hu3rlsGetbias(x = c2, k = k, c1 = c1, c2 = c2), 
                .Hu3rlsGetbias(x = sqrt(beta.c12), k = k, c1 = c1, c2 = c2), 
                .Hu3rlsGetbias(x = sqrt(max(0, beta.c12-0.5*A.loc^2/A.sc^2)), 
                                k = k, c1 = c1, c2 = c2))

    fct1 <- function(x){ A.loc*sign(x)*pmin(abs(x), k) }
    body(fct1) <- substitute({ A.loc*sign(x)*pmin(abs(x), k) },
                        list(k = k, A.loc = A.loc))
    fct2 <- function(x){ Ind1 <- (abs(x)<c1); Ind2 <- (abs(x)>c2)
                         A.sc*(c1^2*Ind1 + x^2*(1-Ind1)*(1-Ind2) + c2^2*Ind2 - beta.c12) }
    body(fct2) <- substitute({ Ind1 <- (abs(x)<c1); Ind2 <- (abs(x)>c2) 
                               A.sc*(c1^2*Ind1 + x^2*(1-Ind1)*(1-Ind2) + c2^2*Ind2 - beta.c12) },
                        list(c1 = c1, c2 = c2, A.sc = A.sc))

    return(IC(name = "IC of Hu3 type", 
              Curve = EuclRandVarList(RealRandVariable(Map = list(fct1, fct2), Domain = Reals())),
              Risks = list(asMSE = res$value, asBias = bias, asCov = res$value - r^2*bias^2), 
              Infos = matrix(c("rlsOptIC.Hu3", "optimally robust IC for Hu3 estimators and 'asMSE'",
                               "rlsOptIC.Hu3", paste("where k =", round(k, 3), ", c1 =", round(c1, 3),
                               "and c2 =", round(c2, 3))), 
                           ncol=2, byrow = TRUE, dimnames=list(character(0), c("method", "message"))), 
              CallL2Fam = call("NormLocationScaleFamily")))
}
