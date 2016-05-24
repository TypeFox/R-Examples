J20 <- function(fx, fy, d){

# Computes second partial of J w.r.t. to fx

z <- rep(0, length(fx))
df <- fx - fy
Q <- exp(-abs(df) / d)

# fx = fy    
noStep <- which(df == 0)
z[noStep] <- exp(fx[noStep]) * (d[noStep] + 1) * (2 * d[noStep] + 1) / (6 * d[noStep])

# fx < fy
JJ <- which(df < 0) 

# Taylor approximation if |1-Q| < 10^-5
TT <- which(abs(1 - Q[JJ]) < 1e-5)
dT <- d[JJ[TT]]
R <- Q[JJ[TT]] - 1
z[JJ[TT]] <- exp(fy[JJ[TT]]) * (40 * dT ^ 2 + 60 * dT + 20 + R * (30 * dT ^ 3 + 60 * dT ^ 2 + 30 * dT) + 
    R ^ 2 * (12 * dT ^ 4 + 15 * dT ^ 3 - 10 * dT ^ 2 - 15 * dT - 2) + R ^ 3 *(3 * dT ^ 5 + dT ^ 4 - 20 * dT ^ 3 + 5 * dT ^ 2 + 17 * dT - 6) 
    + R ^ 4 * (-6 * dT ^ 6 - 5 * (dT ^ 5 + dT ^ 4 - dT ^ 3) - 6 * dT ^ 2)) / (120 * dT)

# exact if |Q - 1| is not small
EE <- which(abs(1 - Q[JJ]) >= 1e-5)    
dE <- d[JJ[EE]]
QE <- Q[JJ[EE]]
z[JJ[EE]] <- exp(fy[JJ[EE]]) * (dE ^ 2 * QE ^ (dE + 3) + (1 - 2 * dE - 2 * dE ^ 2) * QE ^ (dE + 2) + 
    (dE + 1) ^ 2 * QE ^ (dE + 1) - QE ^ 2 - QE) / (dE ^ 2 * (QE - 1) ^ 3)
            
# fx > fy
JJ <- which(df > 0) 

# Taylor approximation if |1 - Q| < 10^-5
TT <- which(abs(1 - Q[JJ]) < 1e-5)
dT <- d[JJ[TT]]
R <- Q[JJ[TT]] - 1
z[JJ[TT]] <- exp(fx[JJ[TT]]) * ((dT + 1) * (2 * dT + 1) / (6 * dT) + R * (dT ^ 2 - 1) / 12 + R ^ 2 * (dT ^ 2 - 1) * (dT - 2) * (2 * dT - 1) / (120 * dT) + 
    R ^ 3 * (dT ^ 2 - 1) * (dT - 2) * (dT - 2) / (360 * dT))

# exact if |Q - 1| is not small
EE <- which(abs(1 - Q[JJ]) >= 1e-5)    
dE <- d[JJ[EE]]
QE <- Q[JJ[EE]]
z[JJ[EE]] <- exp(fx[JJ[EE]]) * (QE ^ (dE + 2) + QE ^ (dE + 1) - (dE + 1) ^ 2 * QE ^ 2 + (2 * dE ^ 2 + 2 * dE - 1) * QE - dE ^ 2) / (dE ^ 2 * (QE - 1) ^ 3)

return(z)
}
