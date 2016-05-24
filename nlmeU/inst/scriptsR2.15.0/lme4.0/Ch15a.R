##  NOTE: This code pertains to panels R15.1 - R15.4

###################################################
### code chunk Chap15init
###################################################
options(width=65, digits=5,show.signif.stars = FALSE)   
date()
packageVersion("Matrix")
sessionInfo()
SeedValue <- 17761
set.seed(SeedValue)

###################################################
### code chunk: R15.1
###################################################
n1 <- 2                 # Number of levels for the factor g1
n2 <- 3                 # Number of levels for the factor g2
i <- gl(n1, n2)         # i 
j <- gl(n2, 1, n1*n2)   # j
b1x <- rnorm(n1, 0, 1)  # b_i
b2x <- rnorm(n2, 0, 2)  # b_j
dt0 <- data.frame(i, j)
(dtc <- 
   within(dt0,
          {             # g1 and g2 are crossed
           eps <- rnorm(nrow(dt0), 0, 0.2)
           b1 <- b1x[i]
           b2 <- b2x[j]
           y <- 10 + b1 + b2 + eps
           g2 <- factor(j, labels = letters[1:n2])
           g1 <- factor(LETTERS[i])
           }))


###################################################
### code chunk: R15.2
###################################################
Zg1 <- model.matrix(~ 0 + g1, data = dtc)  # Z_1 for g1
Zg2 <- model.matrix(~ 0 + g2, data = dtc)  # Z_2 for g2
Z0 <- cbind(Zg1, Zg2)                      # Z for g1 and g2
A0 <- t(Z0)                                # A = Z' 
A0c <- tcrossprod(A0)                      # A*A' 
Dg <- diag(nrow(A0))
(A0q <- A0c + Dg)                          # A*A' + I


###################################################
### code chunk: R15.3a
###################################################
L0 <- t(chol(A0q))           # L such that L*L' = A*A' + I  
sum(L0 != 0.0)               # Count of non-zero elements 
max(abs(L0 %*% t(L0)- A0q))  # Verify L*L' = A*A' + I  


###################################################
### code chunk: R15.3b
###################################################
pvec <- c(3, 4, 5, 1, 2)     # Permutation vector
A1   <- A0[pvec, ]           # Rows permuted in A 
A1c  <- tcrossprod(A1)       
(A1q <- A1c + Dg)            # A*A' + I (permuted)


###################################################
### code chunk: R15.3b  (continued)
###################################################
A1q. <- A0q[pvec, pvec]      # Cols and rows permuted in A * A' 
identical(A1q, A1q.)
L1   <-  t(chol(A1q.))       # L*L' = A*A' + I (permuted) 
sum(L1 != 0.0)               # Count of non-zero elements 




###################################################
### code chunk: R15.4a  
###################################################
library(Matrix)
A0 <- as(A0, "dgCMatrix")           # A0 matrix coerced to sparse
A0c <- tcrossprod(A0)                   # A * A' 
L0 <- Cholesky(A0c, perm = FALSE, Imult = 1, LDL = FALSE)
nnzero(L0. <- as(L0, "sparseMatrix"))   # Coerced to verify 
Dg <- Diagonal(nrow(A0))
(A0q   <- A0c + Dg)
max(abs(L0. %*% t(L0.) - A0q))          # L*L' = A*A' + I 


###################################################
### code chunk: R15.4b  
###################################################
pvec <-  c(3, 4, 5, 1, 2)               # Permutation vector
P1   <-  as(pvec, "pMatrix")            # Permutation matrix
A1c  <-  P1 %*% A0c %*% t(P1)
L1   <-  Cholesky(A1c, perm = FALSE, Imult = 1, LDL = FALSE)
nnzero(as(L1, "sparseMatrix"))


###################################################
### code chunk: R15.4c 
###################################################
L2 <-  Cholesky(A0c, perm=TRUE, Imult =1, LDL = FALSE)
nnzero(as(L2, "sparseMatrix"))
slot(L2,"perm") + 1L                    # Permutation

### SessionInfo 
sessionInfo()            # before detaching package Matrix
detach(package:Matrix)
