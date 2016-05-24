J00 <- function(fx, fy, d){

########################################################################
# Compute auxiliary function J as described in Weyermann (2008):
#
# J00(fx, fy, d) := sum_{j = 0} ^ {d} exp[(1 - j / d) * fx + (j / d) * fy]
# 
# where d = y - x for integers x < y.
# 
# Matlab version 20.11.07, Kathrin Weyermann
# Ported to R by Kaspar Rufibach, October 2010
########################################################################

# initialization
z <- rep(1, length(fx))
df <- fx - fy
Q <- exp(-abs(df) / d)

# fx = fy    
noStep <- which(df == 0)
if (length(noStep) > 0){z[noStep] <- exp(fx[noStep]) * (d[noStep] + 1)}

# fx != fy
isStep <- which(df != 0) 

if (length(isStep) > 0){
    
    # use Taylor approximation if |1 - Q| < 10^-10
    TT <- which(abs(1 - Q[isStep]) < 1e-10)
    if (length(TT) > 0){
        dT <- d[isStep[TT]]
        R <- Q[isStep[TT]] - 1
        z[isStep[TT]] <- exp(pmax(fx[isStep[TT]], fy[isStep[TT]])) * (dT + 1 + R * (dT ^ 2 + dT) / 2 + R ^ 2 * (dT ^ 3 - dT) / 6)
        }

    # exact if |Q - 1| is not small
    EE <- which(abs(1 - Q[isStep]) >= 1e-10)    
    if (length(EE) > 0){
        dE <- d[isStep[EE]]
        Q <- Q[isStep[EE]]
        z[isStep[EE]] <- exp(pmax(fx[isStep[EE]], fy[isStep[EE]]) ) * (1 - Q ^ (dE + 1)) / (1 - Q)
        }
    }
    
return(z)
}
