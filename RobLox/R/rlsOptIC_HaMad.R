###############################################################################
## computation of asymptotic variance
###############################################################################
.HaMadrlsGetvar <- function(a, b, c0){
    h1 <- 2*(-a*dnorm(a) + pnorm(a) - 0.5 + a^2*(pnorm(b)-pnorm(a)) 
             + (a/(c0-b))^2*(c0^2*(pnorm(c0)-pnorm(b)) 
             + 2*c0*(dnorm(c0)-dnorm(b)) - c0*dnorm(c0) + b*dnorm(b) 
             + pnorm(c0) - pnorm(b)))
    A.loc <- 1/(2*(pnorm(a) - 0.5 - a/(c0-b)*(pnorm(c0)-pnorm(b))))
    a.mad <- qnorm(0.75)
    b.mad <- 1/(4*a.mad*dnorm(a.mad))

    return(h1*A.loc^2 + b.mad^2)
}

###############################################################################
## computation of maximum asymptotic MSE
###############################################################################
.HaMadrlsGetmse <- function(abc0, r, MAX){
    a <- abc0[1]; b <- abc0[2]; c0 <- abc0[3]

    #constraints
    if(a < 0 || a > b || b > c0) return(MAX)

    A.loc <- 1/(2*(pnorm(a) - 0.5 - a/(c0-b)*(pnorm(c0)-pnorm(b))))
    a.mad <- qnorm(0.75)
    b.mad <- 1/(4*a.mad*dnorm(a.mad))
    bias.2 <- a^2*A.loc^2 + b.mad^2

    Var <- .HaMadrlsGetvar(a = a, b = b, c0 = c0)

    return(Var + r^2*bias.2)
}

###############################################################################
## optimal IC
###############################################################################
rlsOptIC.HaMad <- function(r, a.start = 0.25, b.start = 2.5, c.start = 5.0, 
                          delta = 1e-6, MAX = 100){
    res <- optim(c(a.start, b.start, c.start), .HaMadrlsGetmse, method = "Nelder-Mead", 
                control = list(reltol=delta), r = r, MAX = MAX)

    a <- res$par[1]; b <- res$par[2]; c0 <- res$par[3]
    A.loc <- 1/(2*(pnorm(a) - 0.5 - a/(c0-b)*(pnorm(c0)-pnorm(b))))
    a.mad <- qnorm(0.75)
    b.mad <- 1/(4*a.mad*dnorm(a.mad))
    bias <- sqrt(a^2*A.loc^2 + b.mad^2)

    fct1 <- function(x){ Ind1 <- (abs(x) < a); Ind2 <- (abs(x) < b); Ind3 <- (abs(x) < c0)
                         A.loc*(x*Ind1 + a*sign(x)*(Ind2-Ind1) + a*(c0-abs(x))/(c0-b)*sign(x)*(Ind3-Ind2))}
    body(fct1) <- substitute({ Ind1 <- (abs(x) < a); Ind2 <- (abs(x) < b); Ind3 <- (abs(x) < c0)
                         A.loc*(x*Ind1 + a*sign(x)*(Ind2-Ind1) + a*(c0-abs(x))/(c0-b)*sign(x)*(Ind3-Ind2)) },
                        list(A.loc = A.loc, a = a, b = b, c0 = c0))
    fct2 <- function(x){ b.mad*sign(abs(x) - a.mad) }
    body(fct2) <- substitute({ b.mad*sign(abs(x) - a.mad) },
                        list(a.mad = a.mad, b.mad = b.mad))

    return(IC(name = "IC of HaMad type", 
              Curve = EuclRandVarList(RealRandVariable(Map = list(fct1, fct2), Domain = Reals())),
              Risks = list(asMSE = res$value, asBias = bias, asCov = res$value - r^2*bias^2), 
              Infos = matrix(c("rlsOptIC.HaMad", "optimally robust IC for HaMad estimators and 'asMSE'",
                               "rlsOptIC.HaMad", paste("where a =", round(a, 3), ", b =", round(b, 3),
                                                     ", c =", round(c0, 3))), 
                           ncol=2, byrow = TRUE, dimnames=list(character(0), c("method", "message"))), 
              CallL2Fam = call("NormLocationScaleFamily")))
}
