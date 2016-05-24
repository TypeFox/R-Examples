# This is file ../spam/tests/demo_article-jss.R
# This file is part of the spam package, 
#      http://www.math.uzh.ch/furrer/software/spam/
# by Reinhard Furrer [aut, cre], Florian Gerber [ctb]
     








# This demo contains the R code to construct the figures and the table of the
# article:
#     "spam: A Sparse Matrix R Package with Emphasis on
#            MCMC Methods for Gaussian Markov Random Fields"
# submitted to JSS.


# The code presented here differs in the following points form the actually used
# one:
# - Very large grid sizes or very high order neighbor structures are not included
#   here;
# - Instead of (100+1) factorizations only (10+1) are performed here;
# - No figure fine-tuning is done here.
# - We had a few additional gc(), just to be sure.  



# The following are tests specific. Not all computers run with profiling. Instead
# of commenting, we define dummies.
options( echo=FALSE)
library( spam, warn.conflict=FALSE)


Rprof <- function(memory.profiling=TRUE, interval=0.1)
  return()
summaryRprof <- function(memory="both")
  return(list(by.total=rbind(1:4)))









# Figure 1:
i <- c( 2,4,4,5,5)
j <- c( 1,1,2,1,3)

A <- spam(0,5,5)
A[cbind(i,j)] <- rep(.5, length(i))
A <- t( A)+A+diag.spam(5)

U <- chol( A)
pivot <- U@pivot
B <- A[pivot,pivot]
R <- chol( B)

U@pivot
U@snmember
U@supernodes

U@entries
U@colindices
U@colpointers
U@rowpointers

if (F){ 
  display( A)
  display( as.spam( chol(as.matrix( A))))
  display( B)
  display( as.spam(R))
  abline( h=-U@supernodes+.5,col=3,lty=2)
}

# Figure 2:
theta1 <- .1
theta2 <- .01
n <- dim( UScounties.storder)[1]

USmat <- diag.spam(n) + theta1 *  UScounties.storder + theta2 *  UScounties.ndorder


U <- chol( USmat,memory=list(nnzR=146735))
if (F) {
  display( as.spam(U))
  text(400,-2200,"MMD\nz=146735\nw=30182\ns=1262",adj=0)
}

U <- chol( USmat, pivot="RCM",memory=list(nnzR=256198,nnzcolindices=140960))
if (F) {
  display( as.spam(U))
  text(400,-2200,"RCM\nz=256198\nw=140960\ns=1706",adj=0)
}

U <- chol( USmat, pivot=FALSE,memory=list(nnzR=689615,nnzcolindices=96463))
if (F) {
  display( as.spam(U))
  text(400,-2200,"no permutation\nz=689615\nw=96463\ns=711",adj=0)
}

# general parameters for the following
N <- 10         # would be 100 in the article 
stsel <- 1      # user.self
rPsx <- 1       # for function "system.time"
rPsy <- 3       # memory usage 
rPint <- .0001  # small interval


# Figure 3:
theta1 <- .1
theta2 <- .05

xseq <- ceiling(4 + exp(seq(0,to=4,by=1))/2)  # would be seq(0.5,to=6,by=.5) in the article
xseql <- length(xseq)

table <- array(NA,c(xseql,4))
for (ix in 1:xseql) {

  egdx <- expand.grid(1:xseq[ix],1:xseq[ix])
  Cspam <- nearest.dist( egdx, delta=1., upper=NULL)
  Dspam <- nearest.dist( egdx, delta=1.5,upper=NULL)

  mat <- diag.spam(xseq[ix]^2) + theta1 * Cspam + theta2 * Dspam

  Rprof( memory.profiling=TRUE, interval = rPint)
  table[ix,1] <- system.time( { ch1 <- chol(mat);
                                for (i in 1:N) ch1 <- chol(mat)}
                             )[stsel]
  Rprof( NULL)
  table[ix,2] <- summaryRprof( memory="both")$by.total[rPsx,rPsy]


  Rprof( memory.profiling=TRUE, interval = rPint)
  table[ix,3] <- system.time( { ch1 <- chol(mat);
                                for (i in 1:N) ch2 <- update(ch1,mat) }
                             )[stsel]
  Rprof( NULL)
  table[ix,4] <- summaryRprof( memory="both")$by.total[rPsx,rPsy]

}

if (F) {
  # Since we have a small N, elements in table might be zero.
  table <- pmax(table, 0.0001)

  par(mfcol=c(1,2))
  plot(xseq, table[,1], type='l', log='xy', ylim=range(table[,c(1,3)]),
       xlab="L (log scale)", ylab="seconds (log scale)")
  lines(xseq, table[,3], lty=2)
  
  plot(xseq, table[,2], type='l', log='xy', ylim=range(table[,c(2,4)]+0.01),
       xlab="L (log scale)", ylab="Mbytes (log scale)")
  lines(xseq, table[,4], lty=2)
}


# Figure 4:

x <- 20     # was 50 in article
maxnn <- 3  # was 6 in article

egdx <- expand.grid( 1:(maxnn+1), 1:(maxnn+1))
dval <- sort(unique(nearest.dist( egdx, delta=maxnn)@entries))
dvall <- length( dval)


egdx <- expand.grid( 1:x, 1:x)


table <- array(NA, c(dvall,5))

for (id in 1:dvall) {

  mat <- nearest.dist( egdx, delta=dval[id],upper=NULL)
  mat@entries <- exp(-2*mat@entries)         # arbitrary values to get a spd precision matrix

  table[id,5] <- length(Cspam)

  Rprof( memory.profiling=TRUE, interval = rPint)
  table[id,1] <- system.time( { ch1 <- chol(mat);
                                for (i in 1:N) ch1 <- chol(mat)}
                             )[stsel]
  Rprof( NULL)
  table[id,2] <- summaryRprof( memory="both")$by.total[rPsx,rPsy]
  
  Rprof( memory.profiling=TRUE, interval = rPint)
  table[id,3] <- system.time( { ch1 <- chol(mat);
                                for (i in 1:N) ch2 <- update(ch1,mat) }
                             )[stsel]
  Rprof( NULL)
  table[id,4] <- summaryRprof( memory="both")$by.total[rPsx,rPsy]

}


if (F) {
  # Since we have a small N, elements in table might be zero.
  table <- pmax(table, 0.0001)

  par(mfcol=c(1,2))
  plot( dval, table[,1], type='l', log='xy',ylim=range(table[,c(1,3)]),
       xlab="distance (log scale)", ylab="seconds (log scale)")
  lines( dval, table[,3],lty=2)
  
  plot( dval, table[,2], type='l', log='xy',ylim=range(table[,c(2,4)]),
       xlab="distance (log scale)", ylab="Mbytes (log scale)")
  lines( dval, table[,4],lty=2)
}



# Table 1:
table <- array(NA,c(9,4))

x <- 10    #  was 50 in article
egdx <- expand.grid(1:x,1:x)

# As above hence shortend
gridmat <- diag.spam(x^2) + .2 * nearest.dist( egdx, delta=1.,upper=NULL) +
  .1 * nearest.dist( egdx, delta=1.5,upper=NULL)
# USmat was constructed above.


# Generic call first:
Rprof( memory.profiling=TRUE, interval = rPint)
table[1,1] <- system.time(   for (i in 1:N) ch1 <- chol(gridmat) )[stsel]
Rprof( NULL)
table[1,2] <- summaryRprof(memory="both")$by.total[rPsx,rPsy]
Rprof( memory.profiling=TRUE, interval = rPint)
table[1,3] <- system.time(   for (i in 1:N) ch2 <- chol(USmat)   )[stsel]
Rprof( NULL)
table[1,4] <- summaryRprof(memory="both")$by.total[rPsx,rPsy]

# Call a chol.spam directly
Rprof( memory.profiling=TRUE, interval = rPint)
table[2,1] <- system.time(   for (i in 1:N) ch1 <- chol.spam(gridmat))[stsel]
Rprof( NULL)
table[2,2] <- summaryRprof(memory="both")$by.total[rPsx,rPsy]
Rprof( memory.profiling=TRUE, interval = rPint)
table[2,3] <- system.time(   for (i in 1:N) ch2 <- chol.spam(USmat)  )[stsel]
Rprof( NULL)
table[2,4] <- summaryRprof(memory="both")$by.total[rPsx,rPsy]

# Less checking:
spam.options( safemode=c(FALSE, FALSE, FALSE))
Rprof( memory.profiling=TRUE, interval = rPint)
table[3,1] <- system.time(   for (i in 1:N) ch1 <- chol( gridmat)    )[stsel]
Rprof( NULL)
table[3,2] <- summaryRprof(memory="both")$by.total[rPsx,rPsy]
Rprof( memory.profiling=TRUE, interval = rPint)
table[3,3] <- system.time(   for (i in 1:N) ch2 <- chol( USmat)    )[stsel]
Rprof( NULL)
table[3,4] <- summaryRprof(memory="both")$by.total[rPsx,rPsy]
spam.options( safemode=c(TRUE, TRUE, TRUE))


# lesser checking
spam.options( cholsymmetrycheck=FALSE)
Rprof( memory.profiling=TRUE, interval = rPint)
table[4,1] <- system.time(   for (i in 1:N) ch1 <- chol( gridmat)    )[stsel]
Rprof( NULL)
table[4,2] <- summaryRprof(memory="both")$by.total[rPsx,rPsy]

Rprof( memory.profiling=TRUE, interval = rPint)
table[4,3] <- system.time(   for (i in 1:N) ch2 <- chol( USmat)    )[stsel]
Rprof( NULL)
table[4,4] <- summaryRprof(memory="both")$by.total[rPsx,rPsy]
spam.options( cholsymmetrycheck=TRUE)

# Pass optimal memory parameters (from above
memory1 = summary(ch1)[1:2]
memory2 = summary(ch2)[1:2]
Rprof( memory.profiling=TRUE, interval = rPint)
table[5,1] <- system.time(   for (i in 1:N) ch1 <- chol( gridmat,memory=memory1)    )[stsel]
Rprof( NULL)
table[5,2] <- summaryRprof(memory="both")$by.total[rPsx,rPsy]
Rprof( memory.profiling=TRUE, interval = rPint)
table[5,3] <- system.time(   for (i in 1:N) ch2 <- chol( USmat,memory=memory2)    )[stsel]
Rprof( NULL)
table[5,4] <- summaryRprof(memory="both")$by.total[rPsx,rPsy]

# All of the above
spam.options( cholsymmetrycheck=FALSE, safemode=c(FALSE,FALSE,FALSE))
Rprof( memory.profiling=TRUE, interval = rPint)
table[6,1] <- system.time(   for (i in 1:N) ch1 <- chol.spam(gridmat,memory=memory1)    )[stsel]
Rprof( NULL)
table[6,2] <- summaryRprof(memory="both")$by.total[rPsx,rPsy]
Rprof( memory.profiling=TRUE, interval = rPint)
table[6,3] <- system.time(   for (i in 1:N) ch2 <- chol.spam(USmat,memory=memory2)    )[stsel]
Rprof( NULL)
table[6,4] <- summaryRprof(memory="both")$by.total[rPsx,rPsy]

# supply the permutation
pivot1 <- ch1@pivot
pivot2 <- ch2@pivot
Rprof( memory.profiling=TRUE, interval = rPint)
table[7,1] <- system.time(   for (i in 1:N) ch1 <- chol.spam(gridmat,pivot=pivot1,
                                                             memory=memory1)    )[stsel]
Rprof( NULL)
table[7,2] <- summaryRprof(memory="both")$by.total[rPsx,rPsy]
Rprof( memory.profiling=TRUE, interval = rPint)
table[7,3] <- system.time(   for (i in 1:N) ch1 <- chol.spam(USmat,pivot=pivot2,
                                                             memory=memory2)    )[stsel]
Rprof( NULL)
table[7,4] <- summaryRprof(memory="both")$by.total[rPsx,rPsy]

# Do not check the permutation
spam.options( cholpivotcheck=FALSE)
Rprof( memory.profiling=TRUE, interval = rPint)
table[8,1] <- system.time(   for (i in 1:N) ch1 <- chol.spam(gridmat,pivot=pivot1,
                                                              memory=memory1)    )[stsel]
Rprof( NULL)
table[8,2] <- summaryRprof(memory="both")$by.total[rPsx,rPsy]

Rprof( memory.profiling=TRUE, interval = rPint)
table[8,3] <- system.time(   for (i in 1:N) ch2 <- chol.spam(USmat,pivot=pivot2,
                                                              memory=memory2)    )[stsel]
Rprof( NULL)
table[8,4] <- summaryRprof(memory="both")$by.total[rPsx,rPsy]

# Update only
Rprof( memory.profiling=TRUE, interval = rPint)
table[9,1] <- system.time(   for (i in 1:N) ch1 <- update(ch1,gridmat)   )[stsel]
Rprof( NULL)
table[9,2] <- summaryRprof(memory="both")$by.total[rPsx,rPsy]

Rprof( memory.profiling=TRUE, interval = rPint)
table[9,3] <- system.time(   for (i in 1:N) ch2 <- update(ch2,USmat)   )[stsel]
Rprof( NULL)
table[9,4] <- summaryRprof(memory="both")$by.total[rPsx,rPsy]


# assemble the table
colnames(table) <- c("grid_time","grid_mem","US_time","US_mem")
rownames(table) <- c("Generic chol","chol.spam","safemode",
                     "symmetrycheck","memory","all","reusing pivot","best cast","update only")


normed.table <- t( round( t(table[-1,])/table[1,],3))

if (F) {
 print( t( round( t(table[-1,])/table[1,],3)))
}




# Figure 5
In <- diag.spam(nrow(UScounties.storder))
struct <- chol(In + .2 * UScounties.storder + .1 * UScounties.ndorder)

len.1 <- 10 # in the article, is set to 180
len.2 <- 5 # in the article, is set to 100
theta.1 <- seq(-.225, to=.515, len=len.1)
theta.2 <- seq(-.09, to=.235, len=len.2)

grid <- array(NA, c(len.1, len.2))
spam.options('cholupdatesingular'='null')

for (i in 1:len.1)
  for(j in 1:len.2) 
    grid[i,j] <- !is.null(update(struct, In + theta.1[i]*UScounties.storder
                       + theta.2[j]* UScounties.ndorder))


options( echo=TRUE)
