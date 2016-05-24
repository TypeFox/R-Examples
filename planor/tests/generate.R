library("planor")
#---------------------------------------------------------------------------
# EXAMPLES OF AUTOMATIC MODEL GENERATION FROM RESOLUTION
#---------------------------------------------------------------------------
M <- planor.model(resolution=4, factors=planor.factors(c(LETTERS[1:4]),  nlevels=rep(2,4)))

# ----------------------------------------------------------
# 4 factors, resolution 4
K <- planor.designkey(factors=c(LETTERS[1:4]),  nlevels=rep(2,4),  nunits=2^3, resolution=4, max.sol=2)
P <- planor.design(key=K, select=1)
print(P)
resum <- summary(K[1])
# ----------------------------------------------------------
# 5 factors, resolution 3
Km <- planor.designkey(factors=c(LETTERS[1:4], "block"),nlevels=rep(2,5),resolution=3,  nunits=2^4)

# ----------------------------------------------------------
# 6 factors, resolution 5
K <- planor.designkey(factors=planor.factors( factors=c(LETTERS[1:6]), nlevels=rep(2,6)),  nunits=2^5, resolution=5, max.sol=2)
P <- planor.design(key=K, select=1)
print(P)
print(summary(K))

