# This file is intented to be called by calex_1d.R.  It creates a
# design matrix for the code, D1.1d, and the field observations,
# D2.1d. 

D1.1d <- latin.hypercube(n1,2)
rownames(D1.1d) <- paste("coderun",1:nrow(D1.1d),sep=".")
colnames(D1.1d) <- c("x",  "A")
head(D1.1d)

D2.1d <- latin.hypercube(n2,1)
rownames(D2.1d) <- paste("obs",1:nrow(D2.1d),sep=".")
colnames(D2.1d) <- c("x")
head(D2.1d)

