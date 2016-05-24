J11 <- function(fx, fy, d){

# Computes first partial derivative of J w.r.t. to fx and fy

z <- rep(0, length(fx))
df <- fx - fy
Q <- exp(-abs(df) / d)

# fx = fy    
noStep <- which(df == 0)
z[noStep] <- exp(fx[noStep]) * (d[noStep] - 1) * (d[noStep] + 1) / (6 * d[noStep])

# fx != fy
isStep <- which(df != 0) 

# Taylor approximation if |1-Q| < 10^-3
TT <- which(abs(1 - Q[isStep]) < 1e-3)
if (length(TT) > 0){
    dT <- d[isStep[TT]]
    R <- Q[isStep[TT]] - 1
    z[isStep[TT]] <- exp(max(fx[isStep[TT]], fy[isStep[TT]])) * 
        ((dT ^ 2 - 1) / (6 * dT) 
        + R * (dT ^ 2 - 1) / 12 
        + R ^ 2 * (dT ^ 2 - 1) * (dT - 2) * (3 * dT + 1) / (120 * dT)
        + R ^ 3 * (dT ^ 2 - 1) * (dT - 2) * (dT - 3) * (2 * dT + 1) / (360 * dT))
    }

# exact if |Q - 1| is not small
EE <- which(abs(1 - Q[isStep]) >= 1e-3)    
if (length(EE) > 0){
    dE <- d[isStep[EE]]
    Q <- Q[isStep[EE]]
    z[isStep[EE]] <- exp(pmax(fx[isStep[EE]], fy[isStep[EE]])) * ((dE - 1) * Q ^ (dE + 2) - (dE + 1) * Q ^ (dE + 1) + (dE + 1) * Q ^ 2 + (1 - dE) * Q) / (dE ^ 2 * (Q - 1) ^ 3) 
    }

return(z)
}
