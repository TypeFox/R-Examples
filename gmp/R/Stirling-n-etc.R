###--------- Stirling numbers of 1st and 2nd kind ----
###          ====================================
### The "double prec" version of this is currently in package 'nacopula'
### (MM: >>> ../../nacopula/R/special-func.R )


##' Compute Stirling numbers of the 1st kind
##'
##' s(n,k) = (-1)^{n-k} times
##' the number of permutations of 1,2,...,n with exactly k cycles
##'
##' NIST DLMF 26.8 --> http://dlmf.nist.gov/26.8
##'
##' @title  Stirling Numbers of the 1st Kind
##' @param n
##' @param k
##' @return s(n,k)
##' @author Martin Maechler
Stirling1 <- function(n,k)
{
    ## NOTA BENE: There's no "direct" method available here
    stopifnot(length(n) == 1, length(k) == 1)
    if (k < 0 || n < k) stop("'k' must be in 0..n !")
    if(n == 0) return(as.bigz(1))
    if(k == 0) return(as.bigz(0))
    S1 <- function(n,k) {
        if(k == 0 || n < k) return(as.bigz(0))
        if(is.na(S <- St[[n]][k])) {
            ## s(n,k) = s(n-1,k-1) - (n-1) * s(n-1,k) for all n, k >= 0
            St[[n]][[k]] <<- S <- if(n1 <- n-1L)
                S1(n1, k-1) - n1* S1(n1, k) else as.bigz(1)
        }
        S
    }
    if(compute <- (nt <- length(St <- get("S1.tab", .Stirl..env))) < n) {
        ## extend the "table":
        length(St) <- n
        for(i in (nt+1L):n) St[[i]] <- rep(as.bigz(NA), i)
    }
    else compute <- is.na(S <- St[[n]][k])
    if(compute) {
        S <- S1(n,k)
        ## store it back:
        assign("S1.tab", St, envir = .Stirl..env)
    }
    S
}

##' Full Vector of Stirling Numbers of the 1st Kind
##'
##' @title  Stirling1(n,k) for all k = 1..n
##' @param n
##' @return the same as sapply(1:n, Stirling1, n=n)
##' @author Martin Maechler
Stirling1.all <- function(n)
{
    stopifnot(length(n) == 1)
    if(!n) return(as.bigz(numeric(0)))
    if(get("S1.full.n", .Stirl..env) < n) {
        assign("S1.full.n", n, envir = .Stirl..env)
        do.call(c, lapply(seq_len(n), Stirling1, n=n))# which fills "S1.tab"
    }
    else get("S1.tab", .Stirl..env)[[n]]
}

##' Compute Stirling numbers of the 2nd kind
##'
##' S^{(k)}_n = number of ways of partitioning a set of $n$ elements into $k$
##'	non-empty subsets
##' (Abramowitz/Stegun: 24,1,4 (p. 824-5 ; Table 24.4, p.835)
##'   Closed Form : p.824 "C."
##'
##' @title  Stirling Numbers of the 2nd Kind
##' @param n
##' @param k
##' @param method
##' @return S(n,k) = S^{(k)}_n
##' @author Martin Maechler, "direct": May 28 1992
Stirling2 <- function(n,k, method = c("lookup.or.store","direct"))
{
    stopifnot(length(n) == 1, length(k) == 1)
    if (k < 0 || n < k) stop("'k' must be in 0..n !")
    method <- match.arg(method)
    switch(method,
           "direct" = {
               sig <- rep(c(1,-1)*(-1)^k, length.out= k+1) # 1 for k=0; -1 1 (k=1)
               k <- 0:k # (!)
               ga <- factorialZ(k) # = gamma(k+1)
               sum( sig * k^n /(ga * rev(ga)))
           },
           "lookup.or.store" = {
               if(n == 0) return(as.bigz(1)) ## else:
               if(k == 0) return(as.bigz(0))
               S2 <- function(n,k) {
                   if(k == 0 || n < k) return(as.bigz(0))
                   if(is.na(S <- St[[n]][k]))
                       ## S(n,k) = S(n-1,k-1) + k * S(n-1,k)   for all n, k >= 0
                       St[[n]][[k]] <<- S <- if(n1 <- n-1L)
                           S2(n1, k-1) + k* S2(n1, k) else as.bigz(1) ## n = k = 1
                   S
               }
               if(compute <- (nt <- length(St <- get("S2.tab", .Stirl..env))) < n) {
                   ## extend the "table":
                   length(St) <- n
                   for(i in (nt+1L):n) St[[i]] <- rep(as.bigz(NA), i)
               }
               else compute <- is.na(S <- St[[n]][k])
               if(compute) {
                   S <- S2(n,k)
                   ## store it back:
                   assign("S2.tab", St, envir = .Stirl..env)
               }
               S
           })
}

##' Full Vector of Stirling Numbers of the 2nd Kind
##'
##' @title  Stirling2(n,k) for all k = 1..n
##' @param n
##' @return the same as sapply(1:n, Stirling2, n=n)
##' @author Martin Maechler
Stirling2.all <- function(n)
{
    stopifnot(length(n) == 1)
    if(!n) return(as.bigz(numeric(0)))
    if(get("S2.full.n", .Stirl..env) < n) {
        assign("S2.full.n", n, envir = .Stirl..env)
        do.call(c, lapply(seq_len(n), Stirling2, n=n))# which fills "S2.tab"
    }
    else get("S2.tab", .Stirl..env)[[n]]
}


##' Compute Eulerian numbers (German "Euler Zahlen")  A(n,k)
##'
##' A(n,k) = number of permutations of n with exactly  k ascents (or k descents)
##' --> http://dlmf.nist.gov/26.14
##'
##' @title Eulerian Numbers
##' @param n
##' @param k
##' @param method
##' @return A(n,k) = < n \\ k >
##' @author Martin Maechler, April 2011
Eulerian <- function(n,k, method = c("lookup.or.store","direct"))
{
    stopifnot(length(n) == 1, length(k) == 1)
    if(k < 0 || n < k) stop("'k' must be in 0..n !")
    if(n && k == n) return(as.bigz(0))
    ## have  __ 0 <= k < n __
    method <- match.arg(method)
    switch(method,
	   "direct" = {
	       if(k == 0) return(as.bigz(1))
	       if(k == n) return(as.bigz(0))
	       ## else 0 <= k < n >= 2
	       ## http://dlmf.nist.gov/26.14.E9 :  A(n,k) = A(n, n-1-k),  n >= 1
	       if(k >= (n+1)%/% 2) k <- n-(k+1L)
	       k1 <- k+1L
	       sig <- rep(c(1,-1), length.out = k1) # 1 for k=0; 1 -1 (k=1), ...
	       sum( sig * chooseZ(n+1, 0:k) * (k1:1L)^n )
	   },
	   "lookup.or.store" = {
	       Eul <- function(n,k) {
		   ## Quick return for those that are *not* stored:
		   if(k < 0 || k > n || (0 < n && n == k)) return(as.bigz(0))
		   if(n == 0) return(as.bigz(1))
		   ## now -- 0 <= k < n -- are stored
		   if(is.na(r <- E.[[n]][k1 <- k+1L])) ## compute it (via recursion)
		       ## A(n,k) = (k+1)* A(n-1,k) + (n-k)*A(n-1,k-1)	for n >= 2, k >= 0
		       ## http://dlmf.nist.gov/26.14.E8
		       E.[[n]][[k1]] <<- r <- if(n1 <- n-1L)
			   k1*Eul(n1, k)+ (n-k)* Eul(n1, k-1) else as.bigz(1) ## n=1, k=0
		   r
	       }
	       if(compute <- (nt <- length(E. <- get("Eul.tab", .Stirl..env))) < n) {
		   ## extend the "table":
		   length(E.) <- n
                   for(i in (nt+1L):n) E.[[i]] <- rep(as.bigz(NA), i)
	       }
	       else compute <- is.na(E <- E.[[n]][k+1L])
	       if(compute) {
		   E <- Eul(n,k)
		   ## store it back:
		   assign("Eul.tab", E., envir = .Stirl..env)
	       }
	       E
	   })
}

##' Full Vector of Eulerian Numbers == row of Euler triangle
##'
##' @title Eulerian(n,k) for all k = 0..n-1
##' @param n
##' @return (for n >= 1), the same as sapply(0:(n-1), Eulerian, n=n)
##' @author Martin Maechler, April 2011
Eulerian.all <- function(n)
{
    stopifnot(length(n) == 1, n >= 0)
    if(!n) return(as.bigz(1))
    if(get("Eul.full.n", .Stirl..env) < n) {
        ## FIXME: do the assign() only when the lapply() below does *not* fail
        ##  on.exit( ... ) ?
	assign("Eul.full.n", n, envir = .Stirl..env)
	do.call(c, lapply(0:(n-1L), Eulerian, n=n))# which fills "Eul.tab"
    }
    else get("Eul.tab", .Stirl..env)[[n]]
}

## Our environment for tables etc:  no hash, as it will contain *few* objects:
.Stirl..env <- new.env(parent=emptyenv(), hash=FALSE)
assign("S2.tab", list(), envir = .Stirl..env) ## S2.tab[[n]][k] == S(n, k)
assign("S1.tab", list(), envir = .Stirl..env) ## S1.tab[[n]][k] == s(n, k)
assign("S2.full.n", 0  , envir = .Stirl..env)
assign("S1.full.n", 0  , envir = .Stirl..env)
assign("Eul.tab", list(), envir = .Stirl..env) ## Eul.tab[[n]][k] == A(n, k) == < n \\ k > (Eulerian)
assign("Eul.full.n", 0	, envir = .Stirl..env)

