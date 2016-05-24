library(abind)
options(error=function() NULL)
x <- array(0, dim=c(2,3,4),dimnames=list(letters[1:2],LETTERS[1:3],letters[23:26]))
# What we want to do here is get slices x[i,c("A","B"),c("w","x","y")] (for all i)
# to be matrix(1:6,ncol=3).
# If we assign in the standard way, and specify submatrices of the same shape in
# indices and as the value of the assign (but leave some indices on the LHS empty),
# then it is possible that submatrices do not end up being assigned to in the way
# we might expect, because the RHS value is flattened and replicated.
# This is only a problem when the repeating index has a lower dim number (here 1)
# than the specified indices (here 2 and 3).
x[,c("A","B"),c("w","x","y")] <- matrix(1:6,ncol=3)
x[1,c("A","B"),c("w","x","y")]
x[2,c("A","B"),c("w","x","y")]
# Assign in a way so that the RHS has its elements laid out in the same way
# as the LHS.
x[,c("A","B"),c("w","x","y")] <- rep(matrix(1:6,ncol=3), each=2)
# first slice
x[1,c("A","B"),c("w","x","y")]
# second slice
x[2,c("A","B"),c("w","x","y")]

# now do it with afill()
x[] <- 0
afill(x, TRUE, , ) <- matrix(1:6,ncol=3, dimnames=list(c("A","B"),c("w","x","y")))
x[1,c("A","B"),c("w","x","y")]
x[2,c("A","B"),c("w","x","y")]
# mix up the order of the RHS of the assignment, afill will sort it back to match the LHS
x[] <- 0
afill(x, T, , ) <- matrix(1:6,ncol=3, dimnames=list(c("A","B"),c("w","x","y")))[2:1,]
x[1,c("A","B"),c("w","x","y")]
x[2,c("A","B"),c("w","x","y")]
table(x==0)

# 4-d example
x <- array(0, dim=c(2,3,3,4),dimnames=list(letters[1:2],LETTERS[1:3],letters[24:26],LETTERS[23:26]))
x[1,c("A","B"),1,c("W","X","Y")] <- 1:6
x[1,c("A","B"),2,c("W","X","Y")] <- 1:6
x[1,c("A","B"),3,c("W","X","Y")] <- 1:6
x[2,c("A","B"),1,c("W","X","Y")] <- 1:6
x[2,c("A","B"),2,c("W","X","Y")] <- 1:6
x[2,c("A","B"),3,c("W","X","Y")] <- 1:6
c(x[1:2,c("A","B"),1:3,c("W","X","Y")])
c(matrix(1:6, ncol=3)[rep(1:2, each=2),rep(1:3,each=3)])

afill(x, T, , T, ) <- matrix(1:6,ncol=3, dimnames=list(c("A","B"),c("W","X","Y")))
x[1,c("A","B"),1,c("W","X","Y")]
x[1,c("A","B"),2,c("W","X","Y")]
x[1,c("A","B"),3,c("W","X","Y")]
x[2,c("A","B"),1,c("W","X","Y")]
x[2,c("A","B"),2,c("W","X","Y")]
x[2,c("A","B"),3,c("W","X","Y")]
table(x==0)

# 2-d example
x <- array(1:24, dim=c(6,4), dimnames=list(LETTERS[1:6], letters[23:26]))
x1 <- x
x1[2:4,2:3] <- -(1:6)
x1
x1 <- x
x1[LETTERS[2:4],letters[24:25]] <- -(1:6)
x1
x2 <- x
afill(x2) <- array(-(1:6),dim=c(3,2), dimnames=list(LETTERS[2:4],letters[24:25]))
x2
identical(x1, x2)
x2 <- x
afill(x2) <- array(-(1:6),dim=c(3,2), dimnames=list(LETTERS[5:7],letters[24:25]))
x2 <- x
afill(x2,excess.ok=T) <- array(-(1:6),dim=c(3,2), dimnames=list(LETTERS[5:7],letters[24:25]))
x2
x2 <- x
afill(x2, local=T) <- array(-(1:6),dim=c(3,2), dimnames=list(LETTERS[2:4],letters[24:25]))
x2

# 1-d named-vector example
x <- c(A=0,B=0,C=0,D=0)
afill(x) <- c(B=1,C=2)
x
# return value is the part of x that is assigned to
(afill(x) <- c(B=1,C=2))
(x[2:3] <- 0)
