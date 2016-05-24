logLikHmm <- function(y,par) {
#
# Function logLikHmm.  To calculate the log likelihood of a sequence,
# or collection (list) of sequences, of observations which come
# from a hidden Markov model with discrete non-parametric observation
# distributions.  These distributions are specified by a matrix
# Rho = [rho_ij] where rho_ij = P(X = x_i | S = j), X being the
# observable random variable and S being the hidden state.

# If y is a matrix, change it to a list, and put out a
# snarky message to the user.
y <- charList(y)

# Get the parameters.
Rho  <- par$Rho
tpm  <- par$tpm
ispd <- par$ispd
if(is.null(ispd)) {
   ispd <- revise.ispd(tpm=tpm)
}

# Make sure that the entries of the vectors in y correspond
# to the row names of Rho.
Rho <- check.yval(y,Rho)

# If K=1 do the triv thing:
K <- length(ispd)
if(K==1) return(sum(log(ffun(y,Rho))))

lns <- sapply(y,length)
fy  <- ffun(y,Rho)
rp  <- recurse(fy,tpm,ispd,lns)
sum(log(rp$llc))
}
