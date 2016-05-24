###################################################
### code chunk Chap14init
###################################################
options(width=65, digits=5, show.signif.stars = FALSE)   
date()
packageVersion("nlmeU")
packageVersion("nlme")
sessionInfo()
data(armd, package = "nlmeU")

###################################################
### code chunk: R14.1a
###################################################
library(nlme)
(pdCS0    <- pdCompSymm(~agex))               
isInitialized(pdCS0)                               # Not initialized        


###################################################
### code chunk: R14.1b
###################################################
mtxUN    <- matrix(c(4, 1, 1, 9), nrow = 2)        # pdSymm matrix
dt1      <- data.frame(agex = c(15, 45, 71, 82))   # Numeric age
(pdSm    <- pdSymm(mtxUN, ~agex, data = dt1))               
isInitialized(pdSm)                             # Initialized


###################################################
### code chunk: R14.1c
###################################################
mtxCS    <- matrix(4 * diag(3) + 1, nrow = 3)      # CompSymm matrix
dt2      <- data.frame(agef = c("Y", "M", "O", "O")) # Factor age 
(pdCSf   <- pdCompSymm(mtxCS, ~-1 + agef, data = dt2))               


###################################################
### code chunk: R14.2a
###################################################
summary(pdSm)               # Summary
formula(pdSm)               # Formula
Names(pdSm)                 # Row/col names
(Dmtx <- as.matrix(pdSm))   # D matrix
Dim(pdSm)                   # Dimensions of D  
logDet(pdSm)                # log|D^(1/2)|

# VarCorr(pdSm)             # Variances, correlation coefficients
# corMatrix(pdSm)           # Corr(D)


###################################################
### code chunk: 14.2b
###################################################
Names(pdCSf)                 # Row/col names
as.matrix(pdCSf)             # D matrix


###################################################
### code chunk: 14.3a
###################################################
coef(pdSm, unconstrained = FALSE)   # Constrained coefficients   
coef(pdSm)                          # Unconstrained coefficients  


###################################################
### code chunk: 14.3b
###################################################
coef(pdCSf, unconstrained = FALSE)  # Constrained coefficients  
coef(pdCSf)                         # Unconstrained coefficients  
log(5)/2                            # First coefficient verified 
rho <- 0.2                          # rho
nc  <- 3                            # No. of columns
aux <- (rho + 1/(nc - 1))/(1 - rho) # Modified Fisher's z: (10.35)
log(aux)                            # Second coefficient verified


###################################################
### code chunk: R14.4a
###################################################
pdSm0 <- pdSymm(mtxUN)
coef(pdSm0)                     # Unconstrained theta_D  
Dmtx <- pdMatrix(pdSm0)         # Matrix D
CholD <- chol(Dmtx)             # Cholesky factor U of D: D=U'U
vd <- svd(CholD, nu=0)          # SVD of U: (13.46)
vd$v %*% (log(vd$d) * t(vd$v))  # (13.47)


###################################################
### code chunk: R14.4b
###################################################
pdLCh <- pdLogChol(mtxUN)
coef(pdLCh)                     # Unconstrained coefficients theta_D  
LChD  <- CholD                  # U
diag(LChD) <- log(diag(LChD))   # \diag(U) log-transformed
LChD


###################################################
### code chunk: R14.4c
###################################################
pdNat <- pdNatural(mtxUN)
coef(pdNat)                     # Unconstrained theta_D 
log(sqrt(diag(Dmtx)))           # log(SDs)
corD  <- cov2cor(Dmtx)          # Corr(D)
rho   <- corD[upper.tri(corD)]  # rho_ij (for i<j) 
log((1+rho)/(1-rho))            # Fisher's z: (10.33)


###################################################
### code chunk: R14.5
###################################################
reSt <- reStruct(list(g1=pdSm,   # D_1
                      g2=pdCSf)) # D_2 
isInitialized(reSt)
names(reSt)                      # Note: order g1, g2 reversed 
formula(reSt)                    # Formulae for pdMat components
getGroupsFormula(reSt)           # Model hierarchy
Names(reSt)                      # Row/col names for pdMat components


###################################################
### code chunk: R14.6a
###################################################
as.matrix(reSt)            # D_1,D_2   
coef(reSt)                 # Unconstrained coefs for D_2,D_1   


###################################################
### code chunk: R14.6b
###################################################
reSt[["g1"]]                     # See pdSm in Panel R14.1b
g2.pdMat  <- reSt[["g2"]]        # See pdCSf in Panel R14.1c
all.equal(pdCSf, g2.pdMat)       # g2.pdMat and pdCSf are equal

###################################################
### code chunk: R14.7
###################################################
Zmtx1 <- model.matrix(formula(pdSm), dt1)                    
prmatrix(Zmtx1)                 # Design matrix Z_1 for pdSm
Zmtx2 <- model.matrix(formula(pdCSf),dt2)                    
prmatrix(Zmtx2)                 # Design matrix Z_2 for pdCSf 
dtz  <- data.frame(dt1,dt2)     # Data frame to evaluate reSt 
Zmtx <- model.matrix(reSt, dtz) # Design matrix Z for reSt 
prmatrix(Zmtx)                  # Matrix Z w/out attributes


###################################################
### code chunk: R14.8
###################################################
reSt  <- reStruct(list(g1 = pdSm, g2 = pdCSf))         # reStruct class
corSt <- corExp(c(0.3, 0.1), form=~tx, nugget = TRUE)  # corStruct class
vF    <- varExp(0.8, form = ~agex)                     # varFunc class
(lmeSt<- lmeStruct(reStruct = reSt, corStruct=corSt,   # lmeStruct class
          varStruct = vF))                             # ... created.
coefs <- coef(lmeSt,unconstrained=FALSE) # Constrained coefficients...
(as.matrix(coefs))                       # ...printed more compactly

#### sessionInfo ###
library(nlmeU)
sessionInfo()           # with packages attached
detach(package:nlmeU)
detach(package:nlme)
