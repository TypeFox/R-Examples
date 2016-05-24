J10 <- function(fx, fy, d){

########################################################################
# Compute first partial derivative of J w.r.t. to fx
#
# J10(fx, fy, d) := sum_{j=0} ^ d (1 - j / d) * exp[(1 - j / d) * fx + (j / d) * fy]
#
# where d = y - x for x < y.
# 
# Matlab version 20.11.07, Kathrin Weyermann
# Ported to R by Kaspar Rufibach, October 2010
########################################################################

z <- rep(0, length(fx))
df <- fx - fy
Q <- exp(-abs(df) / d)

# fx = fy    
noStep <- which(df == 0)
z[noStep] <- exp(fx[noStep]) * (d[noStep] + 1) / 2

# fx < fy
JJ <- which(df < 0)
 
# Taylor approximation if |1-Q| < 10^-5
TT <- which(abs(1 - Q[JJ]) < 1e-5)
if (length(TT) > 0){
    dT <- d[JJ[TT]]
    R <- Q[JJ[TT]] - 1
    z[JJ[TT]] <- exp(fy[JJ[TT]]) * ((dT + 1) / 2 + R * (dT + 1) * (2 * dT + 1) / 6 + R ^ 2 * (dT ^ 2 - 1) * (3 * dT + 2) / 24 + R ^ 3 * (dT * (dT ^ 2 - 1) * (dT - 2)))
    }

# exact if |Q - 1| is not small
EE <- which(abs(1 - Q[JJ]) >= 1e-5) 
if (length(EE) > 0){
    dE <- d[JJ[EE]]
    QE <- Q[JJ[EE]]
    z[JJ[EE]] <- exp(fy[JJ[EE]]) * (dE * QE ^ (dE + 2) - (dE + 1) * QE ^ (dE + 1) + QE) / (dE * (1 - QE) ^ 2)
    }
            
# fx > fy
JJ <- which(df > 0) 

# Taylor approximation if |1 - Q| < 10^-5
TT <- which(abs(1 - Q[JJ]) < 1e-5)
dT <- d[JJ[TT]]
R <- Q[JJ[TT]] - 1
z[JJ[TT]] <- exp(fx[JJ[TT]]) * ((dT + 1) / 2 + R * (dT ^ 2 - 1) / 6 + R ^ 2 * (dT ^ 3 / 2 - dT ^ 2 - dT / 2 + 1) / 12)

# exact if |Q - 1| is not close to 1
EE <- which(abs(1 - Q[JJ]) >= 1e-5)    
dE <- d[JJ[EE]]
QE <- Q[JJ[EE]]
z[JJ[EE]] <- exp(fx[JJ[EE]]) * (QE ^ (dE + 1) - (dE + 1) * QE + dE) / (dE * (1 - QE) ^ 2)
    
return(z)
}
