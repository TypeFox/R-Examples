#<<BEGIN>>
lhs <- function(distr="runif",nsv=ndvar(),nsu=ndunc(),nvariates=1,...)
#TITLE Random Latin Hypercube Sampling
#DESCRIPTION
# Creates a Latin Hypercube Sample (LHS) of the specified distribution.
#KEYWORDS design
#INPUTS
#{distr}<<The function for generating random sample or its name. If \samp{distr} is "rdist",
#the function "qdist" must be the quantile function of this distribution with argument
#\samp{p} as a vector of probabilities, as all univariates distributions of the \samp{stat}
#library.>>
#{nsv}<<The number of raws of the final matrix.>>
#{nsu}<<The number of columns of the final matrix>>
#[INPUTS]
#{nvariates}<<The number of variates>>
#{\dots}<<All arguments to be passed to \samp{distr} except the size of the sample.>>
#VALUE
#A \samp{nsv x nsu} matrix of random variates.
#NOTE
#The resulting lhs is in fact a latin hypersquare sampling: the lhs is provided only in the first 2 dimensions.</>
#It is not possible to send truncated distribution with \code{\link{rtrunc}}. Use \code{\link{mcstoc}} for
#this purpose, with \samp{lhs=TRUE} and \samp{rtrunc=TRUE}.</>
#The \dots arguments will be recycled.
#SEE ALSO
#\code{\link{mcstoc}}
#EXAMPLE
#ceiling(lhs(runif,nsu=10,nsv=10)*10)
#AUTHOR adapted from a code of Rob Carnell (library \samp{lhs})
#CREATED 08-01-25
#--------------------------------------------
{
    nsv
    nsu
    if (length(nsv) != 1 | length(nsu) != 1) stop("nsv and nsu may not be vectors")
    if (any(is.na(c(nsv, nsu)))) stop("nsv and nsu may not be NA or NaN")
    if (any(is.infinite(c(nsv, nsu)))) stop("nsv and nsu may not be infinite")
    if (floor(nsv) != nsv | nsv < 1) stop("nsv must be a positive integer\n")
    if (floor(nsu) != nsu | nsu < 1) stop("nsu must be a positive integer\n")

    arg <- list(...)

    if(!is.character(distr)) distr <- as.character(match.call()$distr)                     #retrieve the name of the function
    distr <- substr(distr, 2, 1000)                             #remove the r
    distr <- paste("q",distr,sep="")                            # add the q

    ranperm <- function(X, N) order(runif(N))
    
    P <- array(dim=c(nsv, nsu, nvariates))
    for(i in 1:nvariates) {
      P[,,i] <- apply(P[,,i,drop=FALSE], 2, ranperm, N = nsv)
      eps <- matrix(runif(nsv * nsu), nrow = nsv, ncol = nsu)
      P[,,i] <- (P[,,i] - 1 + eps) / nsv
      }
    return(as.vector(do.call(distr,c(list(p=P),arg))))
}
