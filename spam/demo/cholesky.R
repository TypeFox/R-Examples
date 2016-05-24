# This is file ../spam/demo/cholesky.R
# This file is part of the spam package, 
#      http://www.math.uzh.ch/furrer/software/spam/
# by Reinhard Furrer [aut, cre], Florian Gerber [ctb]
     








# We illustrate the Cholesky decompostion approaches

set.seed(14)




# first start with a full matrix.
xn <- 750
fmat1 <- matrix(rnorm(xn*xn),xn,xn)
fmat1 <- t( fmat1) %*% fmat1
smat1 <- as.spam(fmat1)
smat2 <- smat1 + diag.spam(xn)

# Generic Cholesky
tmp <- gc(F);Rprof(memory.profiling=TRUE, interval = 0.01)
ch1 <- chol( fmat1)
Rprof(NULL);print( summaryRprof(memory="both")$by.total)


# Sparse Cholesky
tmp <- gc(F);Rprof(memory.profiling=TRUE, interval = 0.01)
ch1 <- chol( smat1)
Rprof(NULL);print( summaryRprof(memory="both")$by.total)

# Sparse Cholesky, direct call
tmp <- gc(F);Rprof(memory.profiling=TRUE, interval = 0.01)
ch1 <- chol.spam( smat1)
Rprof(NULL);print( summaryRprof(memory="both")$by.total)

# Sparse Cholesky, without symmetry check
spam.options(cholsymmetrycheck=FALSE)
tmp <- gc(F);Rprof(memory.profiling=TRUE, interval = 0.01)
ch1 <- chol.spam( smat1)
Rprof(NULL);print( summaryRprof(memory="both")$by.total)


# Sparse Cholesky, reusing pivoting
tmp <- gc(F);Rprof(memory.profiling=TRUE, interval = 0.01)
ch1 <- chol.spam( smat1,pivot=ch1@pivot)
Rprof(NULL);print( summaryRprof(memory="both")$by.total)


# Sparse Cholesky, updating
tmp <- gc(F);Rprof(memory.profiling=TRUE, interval = 0.01)
ch1 <- update.spam.chol.NgPeyton( ch1, smat2)
Rprof(NULL);print( summaryRprof(memory="both")$by.total)

# reset to default
spam.options(cholsymmetrycheck=TRUE)




# now create a sparse matrix.
fmat1[fmat1<3] <- 0
smat1 <- as.spam(fmat1)
smat2 <- smat1 + diag.spam(xn)

# Generic Cholesky
tmp <- gc(F);Rprof(memory.profiling=TRUE, interval = 0.01)
ch1 <- chol( fmat1)
Rprof(NULL);print( summaryRprof(memory="both")$by.total)


# Sparse Cholesky
tmp <- gc(F);Rprof(memory.profiling=TRUE, interval = 0.01)
ch1 <- chol( smat1)
Rprof(NULL);print( summaryRprof(memory="both")$by.total)

# Sparse Cholesky, direct call
tmp <- gc(F);Rprof(memory.profiling=TRUE, interval = 0.01)
ch1 <- chol.spam( smat1)
Rprof(NULL);print( summaryRprof(memory="both")$by.total)

# Sparse Cholesky, without symmetry check
spam.options(cholsymmetrycheck=FALSE)
tmp <- gc(F);Rprof(memory.profiling=TRUE, interval = 0.01)
ch1 <- chol.spam( smat1)
Rprof(NULL);print( summaryRprof(memory="both")$by.total)


# Sparse Cholesky, reusing pivoting
tmp <- gc(F);Rprof(memory.profiling=TRUE, interval = 0.01)
ch1 <- chol.spam( smat1,pivot=ch1@pivot)
Rprof(NULL);print( summaryRprof(memory="both")$by.total)


# Sparse Cholesky, updating
tmp <- gc(F);Rprof(memory.profiling=TRUE, interval = 0.01)
ch1 <- update.spam.chol.NgPeyton( ch1, smat2)
Rprof(NULL);print( summaryRprof(memory="both")$by.total)

# reset to default
spam.options(cholsymmetrycheck=TRUE)



# now create an even sparser matrix.
fmat1 <- fmat1+20*diag(xn)
fmat1[fmat1<32] <- 0
smat1 <- as.spam(fmat1)
smat2 <- smat1 + 1* diag.spam(xn)


# Generic Cholesky
tmp <- gc(F);Rprof(memory.profiling=TRUE, interval = 0.01)
ch1 <- chol( fmat1)
Rprof(NULL);print( summaryRprof(memory="both")$by.total)


# Sparse Cholesky
tmp <- gc(F);Rprof(memory.profiling=TRUE, interval = 0.01)
ch1 <- chol( smat1)
Rprof(NULL);print( summaryRprof(memory="both")$by.total)

# Sparse Cholesky, direct call
tmp <- gc(F);Rprof(memory.profiling=TRUE, interval = 0.01)
ch1 <- chol.spam( smat1)
Rprof(NULL);print( summaryRprof(memory="both")$by.total)

# Sparse Cholesky, without symmetry check
spam.options(cholsymmetrycheck=FALSE)
tmp <- gc(F);Rprof(memory.profiling=TRUE, interval = 0.01)
ch1 <- chol.spam( smat1)
Rprof(NULL);print( summaryRprof(memory="both")$by.total)


# Sparse Cholesky, reusing pivoting
spam.options(cholsymmetrycheck=FALSE)
tmp <- gc(F);Rprof(memory.profiling=TRUE, interval = 0.01)
ch1 <- chol.spam( smat1,pivot=ch1@pivot)
Rprof(NULL);print( summaryRprof(memory="both")$by.total)


# Sparse Cholesky, updating
spam.options(cholsymmetrycheck=FALSE)
tmp <- gc(F);Rprof(memory.profiling=TRUE, interval = 0.01)
ch1 <- update.spam.chol.NgPeyton( ch1, smat2)
Rprof(NULL);print( summaryRprof(memory="both")$by.total)


# reset to default
spam.options(cholsymmetrycheck=TRUE)
