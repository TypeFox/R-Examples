generateLinCon <- function (softinfo, epiparameters, x, explicitBoundConstraints=FALSE, verbose = 0) 
{
#      softinfo: List of softinfo
# epiparameters: list of epiparameters
#             x: vector of data -- unused
# explicitBoundConstraints: if TRUE, add additional constraints by bounding parametets
# using the values in r_lower and r_upper computed in this function. This is for when we
# use solvers (like alabama:::auglag) that don't support passing bounds.
#       verbose: numeric, default 0. For debugging.
#
# Ensure that all required epiparameters are present
epi.required <- c("Ndiscr", "m0", "mN", "order")
epi.missing <- epi.required[!is.element (epi.required, names (epiparameters))]
if (length (epi.missing) > 0) {
    cat ("Missing required epiparameters ", paste (epi.missing, collapse=", "), "\n")
    return (Error = -1)
}
N <- epiparameters$Ndiscr
m0 <- epiparameters$m0
mN <- epiparameters$mN
p <- epiparameters$order
Delta <- (mN-m0)/N
r <- NULL # this is just to pass CMD check?
#
# Set up outputs.
#
A <- NULL   # leq constraint matrix 
Aeq <- NULL # = constraint matrix
b <- NULL   # leq right-hand side 
beq <- NULL # = right-hand side
LARGE_AND_NEGATIVE <- -1000
LARGE_AND_POSITIVE <- 1000
r_lower <- rep (LARGE_AND_NEGATIVE, N+1+(p+1)*N) # lower variable bounds, large but finite
r_upper <- rep (LARGE_AND_POSITIVE, N+1+(p+1)*N) # upper variable bounds

# Ensure that all required soft parameters are present
soft.required <- c("M")
soft.missing <- soft.required[!is.element (soft.required, names (softinfo))]
if (length (soft.missing) > 0) {
    cat ("Missing required softinfo", paste (soft.missing, collapse=", "), "\n")
    return (Error = -1)
}
M <- softinfo$M

#-------------------------------
# Lower semicontinuous
#-------------------------------
if (!is.null(softinfo$semicontinuous) && softinfo$semicontinuous == "lower") {
scoef <- matrix (0, N,N+1)
acoef <- matrix (0, N, (p+1)*N)
constants <- numeric (N)

for (k in 1:N) {
    scoef[k,k] = -1;      
    acoef[k,(p+1)*(k-1)+1] = 1;      
}
if (verbose > 0)
    cat ("Adding ", length (constants), " LEQ constraints in lower sc\n")
A <- rbind (A, cbind (scoef, acoef))
b  <- c(b, constants)

# sk geq sum i=0 to 2 aki(mk - mkminus1)^i for k=1,..., N
scoef <- matrix (0, N, N+1)
acoef <- matrix (0, N,(p+1)*N)
constants <- numeric (N)
for (k in 1:N) {
    scoef[k,k+1] <- -1
    for (i in 0:2) {
        acoef[k,(p+1)*(k-1)+1+i] <- Delta^i
    }
}
if (verbose > 0)
    cat ("Adding ", length (constants), " more LEQ constraints in lower sc\n")

A <- rbind (A, cbind (scoef, acoef))
b  <- c(b, constants)
}

#-------------------------------
# Upper semicontinuous
#-------------------------------

if (!is.null(softinfo$semicontinuous) && softinfo$semicontinuous == "upper") {
# sk-1 leq ak0 for k=1,..., N
scoef <- matrix (0, N, N+1)
acoef <- matrix (0, N, (p+1)*N)
constants <- numeric(N)
for (k in 1:N) {
    scoef[k,k] <- 1      
    acoef[k,(p+1)*(k-1)+1] = -1
}
if (verbose > 0)
    cat ("Adding ", length (constants), " LEQ constraints in upper sc\n")
A <- rbind (A, cbind (scoef, acoef))
b  <- c(b, constants)

# sk leq sum i=0 to 2 aki(mk - mkminus1)^i for k=1,..., N
scoef <- matrix (0, N, N+1)
acoef <- matrix (0, N,(p+1)*N)
constants <- numeric (N)
for (k in 1:N) {
    scoef[k,k+1] = 1
    for (i in 0:2) {
        acoef[k,(p+1)*(k-1)+1+i] = -Delta^i  
    }
}
if (verbose > 0)
    cat ("Adding ", length (constants), " LEQ constraints in upper sc\n")
A <- rbind (A, cbind (scoef, acoef))
b  <- c(b, constants)
}

#--------------------------------
# Pointwise Fisher upper bound
#--------------------------------
if (!is.null (softinfo$pointwiseFisherUpper)) {
# Pointwise 'Fisher' upper bound
# Imposing h'(xi)/h(\xi) leq softinfo.pointwiseFisherUpper at M points in each segment (mkminus1, mk)
# ------------------------------------------------------
# -sum i=1 to 2   i*aki( (j-1)*Delta/(M-1) )^(i-1) leq softinfo.pointwiseFisherUpper for k=1,..., N, j=1,..., M
scoef <- matrix (0, N*M, N+1)
acoef <- matrix (0, N*M,(p+1)*N)
constants <- rep (softinfo$pointwiseFisherUpper, N * M)
for (k in 1:N) {
    for (j in 1:M) {
    for (i in 1:2) {
        acoef[(k-1)*M + j,(p+1)*(k-1)+1+i] <- -i*(  (j-1)*Delta/(M-1)  )^(i-1)
    }
    }
}
A <- rbind (A, cbind (scoef, acoef))
b  <- c(b, constants)
}

#--------------------------------
# Pointwise Fisher lower bound
#--------------------------------
if (!is.null (softinfo$pointwiseFisherLower)) {
# Pointwise 'Fisher' lower bound
# Imposing h'(xi)/h(\xi) geq softinfo.pointwiseFisherLower at M points each segment (mkminus1, mk)
# ------------------------------------------------------
# -sum i=1 to 2 i*aki( (j-1)*Delta/(M-1) )^(i-1) geq softinfo.pointwiseFisherLower for k=1,..., N, j=1,..., M
scoef <- matrix (0, N*M, N+1)
acoef <- matrix (0, N*M,(p+1)*N)
constants <- rep (- softinfo$pointwiseFisherLower, N * M)

for (k in 1:N) {
    for (j in 1:M) {
    for (i in 1:2) {
        acoef[(k-1)*M + j,(p+1)*(k-1)+1+i] = i*( (j-1)*Delta/(M-1) )^(i-1)
    }
    }
}
A <- rbind (A, cbind (scoef, acoef))
b  <- c(b, constants)
}

if (!is.null (softinfo$monotone) && softinfo$monotone == "nonincreasing") {
# Nonincreasing
# ------------------------------------------------------
# first in intervals (mkminus1, mk)
# Imposing h'(xi)/h(\xi) leq 0 at M points in each segment (mkminus1, mk)
# ------------------------------------------------------
# -sum i=1 to 2 i*aki( (j-1)*Delta/(M-1) )^(i-1) leq 0 for k=1,..., N, j=1,..., M
scoef <- matrix (0, N*M, N+1) 
acoef <- matrix (0, N*M,(p+1)*N)
constants <- numeric(N*M)

for (k in 1:N)  {
    for (j in 1:M) {
    		for (i in 1:2) {
                acoef[(k-1)*M + j,(p+1)*(k-1)+1+i] <- 
                          -i*( (j-1)*Delta/(M-1) )^(i-1)
          }
    }
}
if (verbose > 0)
    cat ("Adding ", length (constants), " LEQ constraints in nonincr\n")
A <- rbind (A, cbind (scoef, acoef))
b  <- c(b, constants)

# second at the end points; same as usc epispline constraint
# sk-1 leq ak0 for k=1,..., N

scoef <- matrix (0,N,N+1)
acoef <-  matrix(0, N,(p+1)*N)
constants <-  numeric (N)
for (k in 1:N) {
    scoef[k,k] <- 1   
    acoef[k,(p+1)*(k-1)+1] <- -1
}    
if (verbose > 0)
    cat ("Adding ", length (constants), " more LEQ constraints in nonincr\n")
A <- rbind (A, cbind (scoef, acoef))
b <- c(b, constants)

# sk geq sum i=0 to 2 aki(mk - mkminus1)^i for k=1,..., N

scoef <- matrix (0,N,N+1)
acoef <- matrix (0, N,(p+1)*N)
constants <- numeric(N)
for (k in 1:N) {
    scoef[k,k+1] <- -1
    for (i in 0:2) {
        acoef[k,(p+1)*(k-1)+1+i] <- Delta^i;      
    }
}
if (verbose > 0)
    cat ("Adding ", length (constants), " yet more LEQ constraints in nonincr\n")

A <- rbind (A, cbind (scoef, acoef))
b <- c (b, constants)

} # end "nonincreasing"
else if (!is.null (softinfo$monotone) && softinfo$monotone == "nondecreasing") {
scoef <- matrix (0, N*M, N+1) 
acoef <- matrix (0, N*M,(p+1)*N)
constants <- numeric(N*M)

for (k in 1:N)  {
    for (j in 1:M) {
    		for (i in 1:2) {
                acoef[(k-1)*M + j,(p+1)*(k-1)+1+i] <- 
                          i*( (j-1)*Delta/(M-1) )^(i-1)
          }
    }
}
if (verbose > 0)
    cat ("Adding ", length (constants), " LEQ constraints in nondecr\n")

     
A <- rbind (A, cbind (scoef, acoef))
b  <- c(b, constants)

# second at the end points; same as lsc epispline constraint
# sk-1 geq ak0 for k=1,..., N

scoef <- matrix (0,N,N+1)
acoef <-  matrix(0, N,(p+1)*N)
constants <-  numeric (N)
for (k in 1:N) {
    scoef[k,k] <- -1   
    acoef[k,(p+1)*(k-1)+1] <- 1
}    
if (verbose > 0)
    cat ("Adding ", length (constants), " more LEQ constraints in nondecr\n")


A <- rbind (A, cbind (scoef, acoef))
b <- c(b, constants)

# sk leq sum i=0 to 2 aki(mk - mkminus1)^i for k=1,..., N

scoef <- matrix (0,N,N+1)
acoef <- matrix (0, N,(p+1)*N)
constants <- numeric(N)
for (k in 1:N) {
    scoef[k,k+1] <- 1
    for (i in 0:2) {
        acoef[k,(p+1)*(k-1)+1+i] <- -Delta^i;      
    }
}
if (verbose > 0)
    cat ("Adding ", length (constants), " yet more LEQ constraints in nondecr\n")


A <- rbind (A, cbind (scoef, acoef))
b <- c (b, constants)

} # end "nondecreasing"

if (!is.null (softinfo$unimodal) && softinfo$unimodal == TRUE)
{
# Unimodal
# ------------------------------------------------------
# continuous
#sk-1 = ak0 for k=1,..., N
scoef <- matrix (0, N,N+1)
acoef <- matrix (0, N,(p+1)*N)
constants <- numeric (N)
for (k in 1:N){
    scoef[k,k] <- 1
    acoef[k,(p+1)*(k-1)+1] <- -1
}
if (verbose > 0)
    cat ("Adding ", length (constants), "  EQ constraints in unimodal\n")


Aeq <- rbind (Aeq, cbind (scoef, acoef))
beq <- c (beq, constants)

#sk = sum i=0 to 2 aki(mk - mkminus1)^i for k=1,..., N
scoef <- matrix (0, N,N+1)
acoef <- matrix (0, N,(p+1)*N)
constants <- numeric(N)
for (k in 1:N) {
    scoef[k,k+1] = 1
    for (i in 0:2)
        acoef[k,(p+1)*(k-1)+1+i] = -Delta^i
}
if (verbose > 0)
    cat ("Adding ", length (constants), " more EQ constraints in unimodal\n")

Aeq <- rbind (Aeq, cbind (scoef, acoef))
beq <- c (beq, constants)

# left-derivative
#sum i=1 to 2 i*aki(mk - mkminus1)^(i-1) leq ak+1,1    for k=1,..., N-1
scoef <- matrix (0, N-1, N+1)
acoef <- matrix(0, N-1, (p+1)*N) 
constants <- numeric (N-1)
for (k in 1:(N-1)){
    for (i in 1:2){ 
        acoef[k,(p+1)*(k-1)+1+i] <- i * Delta^(i-1)
        acoef[k,(p+1)*k+2] <- -1
    }
}
if (verbose > 0)
    cat ("Adding ", length (constants), "  LEQ constraints in unimodal\n")

A <- rbind (A, cbind (scoef, acoef))
b <- c (b, constants)


# convexity
# for p=2 this implies just nonnegativity of ak2 for k=1, ..., N, which is
scoef_lower <- rep (LARGE_AND_NEGATIVE, N+1)     # was "-Inf *"
acoef_lower <- rep (LARGE_AND_NEGATIVE, (p+1)*N) # was "-Inf *"
for (k in 1:N)
        acoef_lower[(k-1)*(p+1)+p+1] <- 0
# r_lower = (max([r_lower'; [scoef_lower; acoef_lower]']))';

r_lower <- c(scoef_lower, acoef_lower)

}

#-------------------------------
# Unimodal upper tail
#-------------------------------
if (!is.null (softinfo$unimodaluppertail) && softinfo$unimodaluppertail > 0) {
# Upper unimodal
# Only applies to floor(N*softinfo.unimodaluppertail) segments
# ------------------------------------------------------
Nred <- floor( N * softinfo$unimodaluppertail)
if (Nred > 0) {
    # continuous
    #sk-1 = ak0 for k=1,..., N
    scoef <- matrix (0, Nred, N+1)
    acoef <- matrix (0, Nred, (p+1) * N)
    constants <- numeric (Nred)
    for (k in (N-Nred+1):N) {
        scoef[k-N+Nred,k] <- 1
        acoef[k-N+Nred,(p+1)*(k-1)+1] <- -1
    }
    Aeq <- rbind (Aeq, cbind (scoef, acoef))
    beq <- c (beq, constants)

    #sk = sum i=0 to 2 aki(mk - mkminus1)^i for k=1,..., N
    scoef <- matrix (0, Nred, N+1)
    acoef <- matrix (0, Nred, (p+1) * N)
    constants <- numeric (Nred)
    for (k in (N-Nred+1):N) {
        scoef[k-N+Nred,k+1] <- 1
        for (i in 0:2) {
            acoef[k-N+Nred,(p+1)*(k-1)+1+i] <- -Delta^i
        }
    }
    Aeq <- rbind (Aeq, cbind (scoef, acoef))
    beq <- c (beq, constants)

    # left-derivative
    #sum i=1 to 2 i*aki(mk - mkminus1)^(i-1) leq ak+1,1    for k=1,..., N-1
    scoef <- matrix (0, Nred-1, N+1)
    acoef <- matrix (0, Nred-1,(p+1) * N)
    constants <- numeric(Nred-1)
    for (k in (N-Nred+1):(N-1)) {
        for (i in 1:2) {
            acoef[k-N+Nred,(p+1)*(k-1)+1+i] <- i * Delta^(i-1)
            acoef[k-N+Nred,(p+1)*k+2] <- -1   
        }
    }
    A <- rbind (A, cbind (scoef, acoef))
    b <- c (b, constants)

    # convexity
    # for p=2 this implies just nonnegativity of ak2 for k=1, ..., N, which is
    scoef_lower <- rep (LARGE_AND_NEGATIVE, N+1)
    acoef_lower <- rep (LARGE_AND_NEGATIVE, (p+1)*N)
    for (k in (N-Nred+1):N) {
        acoef_lower[(k-1)*(p+1)+p+1,1] <- 0
    }
    r_lower <- pmax(r_lower, c(scoef_lower, acoef_lower))
}
}

#-----------------------------
# Unimodal lower tail
#-----------------------------
if (!is.null (softinfo$unimodallowertail) && softinfo$unimodallowertail > 0) {
# Lower unimodal
# Only applies to floor(N*softinfo.unimodallowertail) segments
# ------------------------------------------------------
Nred <- floor(N * softinfo$unimodallowertail)
if (Nred > 0) {
    # continuous
    #sk-1 = ak0 for k=1,..., N
    scoef <- matrix (0, Nred, N+1)
    acoef <- matrix (0, Nred, (p+1)*N)
    constants <- numeric (Nred)
    for (k in 1:Nred) {
        scoef[k,k] <- 1
        acoef[k,(p+1)*(k-1)+1] <- -1     
    }
    Aeq <- rbind (Aeq, cbind (scoef, acoef))
    beq <- c (beq, constants)

    #sk = sum i=0 to 2 aki(mk - mkminus1)^i for k=1,..., N
    scoef <- matrix (0, Nred, N+1)
    acoef <- matrix (0, Nred, (p+1) * N)
    constants <- numeric (Nred)
    for (k in 1:Nred) {
        scoef[k,k+1] <- 1
        for (i in 0:2) {
            acoef[k,(p+1)*(k-1)+1+i] <- -Delta^i
        }
    }
    Aeq <- rbind (Aeq, cbind (scoef, acoef))
    beq <- c (beq, constants)

    # left-derivative
    #sum i=1 to 2 i*aki(mk - mkminus1)^(i-1) leq ak+1,1    for k=1,..., N-1
    scoef <- matrix (0, Nred-1, N+1)
    acoef <- matrix (0, Nred-1,(p+1) * N)
    constants <- numeric (Nred-1)
    for (k in 1:(Nred-1)) {
        for (i in 1:2) {
            acoef[k,(p+1)*(k-1)+1+i] <- i * Delta^(i-1)
            acoef[k,(p+1)* k + 2] <- -1
        }
    }
    A <- rbind (A, cbind (scoef, acoef))
    b <- c (b, constants)


    # convexity
    # for p=2 this implies just nonnegativity of ak2 for k=1, ..., N, which is
    scoef_lower <- rep (LARGE_AND_NEGATIVE, N + 1)
    acoef_lower <- rep (LARGE_AND_NEGATIVE, (p+1) * N)
    for (k in 1:Nred) {
        acoef_lower[(k-1)*(p+1)+p+1, 1] <- 0
    }
    r_lower = pmax(r_lower, c(scoef_lower, acoef_lower))
}
}
#
# Don't impose continuity if unimodal is TRUE). But do impose it if
# continuous isn't set, if unimodal is FALSE and continuousDiff is TRUE.
#
if ((is.null (softinfo$unimodal) || softinfo$unimodal == FALSE)
&& (  !is.null (softinfo$continuous) && softinfo$continuous == TRUE
   || !is.null (softinfo$continuousDiff) && softinfo$continuousDiff == TRUE))
{
# continuity
# ------------------------------------------------------
# sk-1 = ak0 for k=1,..., N
scoef <- matrix (0, N, N+1)
acoef <- matrix (0, N, (p+1) * N)
constants <- numeric (N)

for (k in 1:N) {
    scoef[k,k] <- 1
    acoef[k,(p+1)*(k-1)+1] <- -1
}
if (verbose > 0)
    cat ("Adding ", length (constants), "  EQ constraints in continuous\n")


Aeq <- rbind (Aeq, cbind (scoef, acoef))
beq <- c(beq, constants)

#sk = sum i=0 to 2 aki(mk - mkminus1)^i for k=1,..., N
scoef <- matrix(0, N,N+1)
acoef <- matrix(0, N,(p+1)*N)
constants <- numeric(N)
for (k in 1:N) {
    scoef[k,k+1] <- 1
    for (i in 0:2)
        acoef[k,(p+1)*(k-1)+1+i] = -Delta^i
}
if (verbose > 0)
    cat ("Adding ", length (constants), "  EQ constraints in continuous\n")

Aeq <- rbind (Aeq, cbind (scoef, acoef))
beq <- c(beq, constants)
} # end continuity

if (!is.null(softinfo$continuousDiff) && softinfo$continuousDiff == TRUE)
{
#---------------------------------------------
# Continuously differentiable
#---------------------------------------------
# ak+1,1 = sum i=1 to 2 i*aki(mk - mkminus1)^(i-1) for k=1,..., N-1
scoef <- matrix (0, N-1,N+1)
acoef <- matrix (0, N-1,(p+1)*N)
constants <- numeric(N-1)
for (k in 1:(N-1)) {
    for (i in 1:2)
        acoef[k,(p+1)*(k-1)+1+i] <- i*Delta^(i-1)
    acoef[k,(p+1)*k+2] = -1
}
if (verbose > 0)
    cat ("Adding ", length (constants), "  EQ constraints in contdiff\n")


Aeq <- rbind (Aeq, cbind (scoef, acoef))
beq <- c(beq, constants)
}

#--------------------------------------------------
# Upper bounds on density values
#--------------------------------------------------
for (k in 1:N) {
if (!is.null(softinfo$upperdensityvalue) && softinfo$upperdensityvalue[k] < Inf) {
# Upper bounds on density values
# ------------------------------------------------------
# sum i=0 to 2 aki(  (j-1)*Delta/(M-1)   )^(i) geq -log softinfo$upperdensityvalue(k) for j=1,..., M
scoef <- matrix (0, M, N+1)
acoef <- matrix (0, M,(p+1)*N)
constants <- numeric (M)
for (j in 1:M) {
    for (i in 0:2) {
        acoef[j,(p+1)*(k-1)+1+i] <- -( (j-1)*Delta/(M-1) )^i
    }
    constants[j] <- log(softinfo$upperdensityvalue[k]);
}
A <- rbind (A, cbind (scoef, acoef))
b <- c (b, constants)
} # end "if"
} # end "for k"

# sk geq -log softinfo.upperdensityvalueEndpt(k) for k=0,1,..., N
#----------------------------------------------
# upperdensityvalueEndpt can come in as Inf; log hates that. Adjust. Where we have Inf,
# just use scoef_lower's pre-existing values. Otherwise find the max of that and -log(Endpt).
# Finally replace the first half of r_lower with these new values.
#

if (!is.null (softinfo$upperdensityvalueEndpt)) {
    scoef_lower <- numeric (length (softinfo$upperdensityvalueEndpt))
    infy <- is.infinite (softinfo$upperdensityvalueEndpt) | softinfo$upperdensityvalueEndpt == LARGE_AND_POSITIVE
    scoef_lower[infy] <- r_lower[1:(N+1)][infy]
    scoef_lower[!infy] <- pmax(r_lower[1:(N+1)][!infy], -log(softinfo$upperdensityvalueEndpt[!infy]))
    r_lower  <- c(scoef_lower, r_lower[(N+2):length(r_lower)])

}

#--------------------------------------------------
# Lower bounds on density values
#--------------------------------------------------
for (k in 1:N) {
if (!is.null(softinfo$lowerdensityvalue) && softinfo$lowerdensityvalue[k] > 0) {
# Lower bounds on density values
# ------------------------------------------------------
# sum i=0 to 2 aki(  (j-1)*Delta/(M-1)   )^(i) leq -log softinfo$lowerdensityvalue(k) for j=1,..., M
scoef <- matrix (0, M, N+1)
acoef <- matrix (0, M,(p+1)*N)
constants <- numeric (M)
for (j in 1:M) {
    for (i in 0:2) {
        acoef[j,(p+1)*(k-1)+1+i] <- ( (j-1)*Delta/(M-1) )^i
    }
    constants[j] <- -log(softinfo$lowerdensityvalue[k]);
}
A <- rbind (A, cbind (scoef, acoef))
b <- c (b, constants)
} # end "if"
} # end "for k"

# sk leq -log softinfo.upperdensityvalueEndpt(k) for k=0,1,..., N
#----------------------------------------------
if (any (names (softinfo) == "lowerdensityvalueEndpt")) { # never Inf
    scoef_upper <- pmin(r_upper[1:(N+1)], -log(softinfo$lowerdensityvalueEndpt))
    r_upper  <- c(scoef_upper, r_upper[(N+2):length(r_upper)])
}
#--------------------------------------------------
# Lower bounds on density values at specific points
#---------------------------------------------


if (!is.null(softinfo$lowerdensityvalueSpecific)) {
# Lower bounds on density values at specific points 
# ------------------------------------------------------
# lowerdensityvalueSpecific is a two-column matrix whose first column gives
# x values and whose second column gives density minima
npts <- nrow (softinfo$lowerdensityvalueSpecific) # number of points with density bounds
for (j in 1:npts) {
    scoef <- matrix (0, 1, N+1)
    acoef <- matrix (0, 1, (p+1)*N)
    for (k in 0:N) { # check if at mesh point
        if (softinfo$lowerdensityvalueSpecific[j,1] == m0+k*Delta) {
            scoef[k+1] <- 1
            break
        }
    }
    for (k in 1:N) { # check if in segments
        if (softinfo$lowerdensityvalueSpecific[j,1] > m0+(k-1)*Delta && softinfo$lowerdensityvalueSpecific[j,1] < m0+k*Delta) {
            for (i in 0:2) {
                acoef[1,(p+1)*(k-1)+1+i] = ( softinfo$lowerdensityvalueSpecific[j,1] - (m0 + (k-1)*Delta) )^i
            }
            break
        }
    }
    constants = -log(softinfo$lowerdensityvalueSpecific[j,2])
    A <- rbind (A, cbind (scoef, acoef))
    b <- c (b, constants)  
}
}


if (!is.null (softinfo$KLDivergenceUpper))
{
#---------------------------------------------
# Kullback Liebler divergence, upper
#---------------------------------------------
tol <- 1e-10
# compute constant term

constterm <- pracma::quadv (llfcn, m0, mN, tol, epiparameters = epiparameters, 
                                             softinfo = softinfo)
acoef <- pracma::quadv (llfcn.vec, m0, mN, tol, epiparameters = epiparameters, 
                                             softinfo = softinfo)

scoef <- numeric (N + 1)
if (verbose > 0)
    cat ("Adding 1 (?) LEQ in KL Upper\n")
A <- rbind (A, c (scoef, acoef$Q))
b <- c(b, softinfo$KLDivergenceUpper - constterm$Q)
}

if (!is.null (softinfo$KLDivergenceLower))
{
#---------------------------------------------
# Kullback Liebler divergence, lower
#---------------------------------------------
tol <- 1e-10
# compute constant term

constterm <- pracma::quadv (llfcn, m0, mN, tol, epiparameters = epiparameters, 
                                             softinfo = softinfo)
acoef <- pracma::quadv (llfcn.vec, m0, mN, tol, epiparameters = epiparameters, 
                                             softinfo = softinfo)

scoef <- numeric (N + 1)
if (verbose > 0)
    cat ("Adding 1 (?) LEQ in KL Upper\n")
A <- rbind (A, c (scoef, acoef$Q))
b <- c(b, softinfo$KLDivergenceLower + constterm$Q)
}

# upper bound on segment endpoint values sk 
# -------------------------------------------------------------------------
# Use the smaller of the existing r_upper and the upperboundsk vector
scoef_upper <- pmin(r_upper[1:(N+1)], softinfo$upperboundsk)
r_upper <- c(scoef_upper, r_upper[(N+2):length(r_upper)])

# lower bound on segment endpoint values sk 
# -------------------------------------------------------------------------
scoef_lower <- pmax(r_lower[1:(N+1)], softinfo$lowerboundsk)
r_lower  <- c(scoef_lower, r_lower[(N+2):length(r_lower)])


# upper bound on constant ak0 in segment k 
# -------------------------------------------------------------------------
for (k in 1:N) {
    if (r_upper[N+1 + (k-1)*(p+1)+1 ] > softinfo$upperboundak0[k])
        r_upper[N+1 + (k-1)*(p+1)+1 ] <- softinfo$upperboundak0[k]  
}


# lower bound on constant ak0 in segment k 
# -------------------------------------------------------------------------
for (k in 1:N) {
    if (r_lower[N+1 + (k-1)*(p+1)+1 ] < softinfo$lowerboundak0[k])
        r_lower[N+1 + (k-1)*(p+1)+1 ] <- softinfo$lowerboundak0[k]
}


# upper bound on polynomial coefficients akp in segment k, all p >0 treated
# equally
# -------------------------------------------------------------------------
for (k in 1:N) {
    for (j in 1:p) {
        if (r_upper[N+1 + (k-1)*(p+1)+1+j ] > softinfo$upperboundakp[k])
            r_upper[N+1 + (k-1)*(p+1)+1+j ] <- softinfo$upperboundakp[k]
    }
}


# lower bound on polynomial coefficients akp in segment k, all p >0 treated
# equally
# -------------------------------------------------------------------------
for (k in 1:N) {
    for (j in 1:p) {
        if (r_lower[N+1 + (k-1)*(p+1)+1+j ] < softinfo$lowerboundakp[k])
            r_lower[N+1 + (k-1)*(p+1)+1+j ] <- softinfo$lowerboundakp[k]
    }
}


#--------------------------------------------------------------------------
# Explicit Bound constraints
#--------------------------------------------------------------------------
#
# If we're using alabama/auglag, we have to add bounds here, as
# regular linear inequality constraints. Let's impose only the ones
# that are other than LARGE_AND_NEGATIVE, LARGE_AND_POSITIVE
#
if (explicitBoundConstraints) {
    lowers <- which (r_lower > LARGE_AND_NEGATIVE)
    if (any (lowers))
    {
        A.supp <- matrix (0, nrow = length (lowers), ncol = length (r))
        A.supp[cbind (1:length (lowers), lowers)] <- -1 
        A <- rbind (A, A.supp)
        b <- c(b, r_lower[lowers])
    }

    uppers <- which (r_upper < LARGE_AND_POSITIVE)
    if (any (uppers))
    {
        A.supp <- matrix (0, nrow = length (uppers), ncol = length (r))
        A.supp[cbind (1:length (uppers), uppers)] <- +1 
        A <- rbind (A, A.supp)
        b <- c(b, r_upper[uppers])
    }
#
# Key Point here: the A constraints are constructed to be "<=" ones, but
# alabama::auglag wants them to be ">=". So multiply by -1. Of course we don't
# have to do this with the Aeq ones.
#
    A <- -A
}


return (list (A = A, b = b, Aeq = Aeq, beq = beq, r_lower = r_lower, r_upper = r_upper))

}
