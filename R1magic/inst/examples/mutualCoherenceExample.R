set.seed(42)
A <- DFTMatrix0(10) # Fourier Bases
B <- matrix(rnorm(100), 10, 10) # Gaussian Random Matrix
C <- A %*% B # A sensing matrix with A and B as above
mutualCoherence(A, 3)
mutualCoherence(B, 3)
mutualCoherence(C, 3)
aa<-mutualCoherence(A, 8)
bb<-mutualCoherence(A, 8)
bb<-mutualCoherence(A, 8)

# Below is to test correctness of mutualCoherence
# Manual Example to check mutualCoherence
set.seed(42)
A <- matrix(rnorm(16), 4, 4)
# k =1 
#  Qs =  1, 2, 3, 4
Q1p1 = .oneMutualcolIndex(A, 1, 2)
Q2p1 = .oneMutualcolIndex(A, 1, 3)
Q3p1 = .oneMutualcolIndex(A, 1, 4)
Qp1  = max(c(Q1p1, Q2p1, Q3p1))

Q1p2 = .oneMutualcolIndex(A, 2, 1)
Q2p2 = .oneMutualcolIndex(A, 2, 3)
Q3p2 = .oneMutualcolIndex(A, 2, 4)
Qp2  = max(c(Q1p2, Q2p2, Q3p2))

Q1p3 = .oneMutualcolIndex(A, 3, 1)
Q2p3 = .oneMutualcolIndex(A, 3, 2)
Q3p3 = .oneMutualcolIndex(A, 3, 4)
Qp3  = max(c(Q1p3, Q2p3, Q3p3))

Q1p4 = .oneMutualcolIndex(A, 4, 1)
Q2p4 = .oneMutualcolIndex(A, 4, 2)
Q3p4 = .oneMutualcolIndex(A, 4, 3)
Qp4  = max(c(Q1p4, Q2p4, Q3p4))

mcK1 = max(c(Qp1, Qp2, Qp3, Qp4))


# k = 2
#  Qs =  (1, 2), (1, 3), (1, 4), (2, 3), (3, 4)
Q1p1 = .oneMutualcolIndex(A, 1, 2)
Q2p1 = .oneMutualcolIndex(A, 1, 3)
Q3p1 = .oneMutualcolIndex(A, 1, 4)
Q4p1 = .oneMutualcolIndex(A, 1, 2) +  .oneMutualcolIndex(A, 1, 3)
Q5p1 = .oneMutualcolIndex(A, 1, 4) +  .oneMutualcolIndex(A, 1, 3)
Qp1  = max(c(Q1p1, Q2p1, Q3p1, Q4p1, Q5p1))

Q1p2 = .oneMutualcolIndex(A, 2, 1)
Q2p2 = .oneMutualcolIndex(A, 2, 3) + .oneMutualcolIndex(A, 2, 1)
Q3p2 = .oneMutualcolIndex(A, 2, 4) + .oneMutualcolIndex(A, 2, 1)
Q4p2 = .oneMutualcolIndex(A, 2, 3)
Q5p2 = .oneMutualcolIndex(A, 2, 4) + .oneMutualcolIndex(A, 2, 3)
Qp2  = max(c(Q1p2, Q2p2, Q3p2, Q4p2, Q5p2))

Q1p3 = .oneMutualcolIndex(A, 3, 1) + .oneMutualcolIndex(A, 3, 2)
Q2p3 = .oneMutualcolIndex(A, 3, 1) 
Q3p3 = .oneMutualcolIndex(A, 3, 4) + .oneMutualcolIndex(A, 3, 1)
Q4p3 = .oneMutualcolIndex(A, 3, 2)
Q5p3 = .oneMutualcolIndex(A, 3, 4)
Qp3  = max(c(Q1p3, Q2p3, Q3p3, Q4p3, Q5p3))

Q1p4 = .oneMutualcolIndex(A, 4, 1) + .oneMutualcolIndex(A, 4, 2)
Q2p4 = .oneMutualcolIndex(A, 4, 1) + .oneMutualcolIndex(A, 4, 3)
Q3p4 = .oneMutualcolIndex(A, 4, 1) 
Q4p4 = .oneMutualcolIndex(A, 4, 2) + .oneMutualcolIndex(A, 4, 3)
Q5p4 = .oneMutualcolIndex(A, 4, 3)
Qp4  = max(c(Q1p4, Q2p4, Q3p4, Q4p4, Q5p4))

mcK2 = max(c(Qp1, Qp2, Qp3, Qp4))
