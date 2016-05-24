
# matrix powers, with flexible options and some callback apps

# matpow() arguments:
#     m:  the (square) input matrix
#     k:  the desired matrix power; if NULL, must be set by callback
#     squaring:  if TRUE, use the squaring method, i.e. form the 
#                square of m, then square that to get the fourth 
#                power of m, etc., so that only about log k 
#                iterations are needed; note: if k is not a power of 2,
#                a higher power than k will be computed
#     genmulcmd:  an R function that generates a multiplication 
#                 command in string form; formal arguments are
#                 a, b and c, quoted strings; places product of a and b
#                 into c; see genmulcmd.vanilla()
#                 for the ordinary matrix case
#     dup:  a function that will form a deep copy of a matrix of 
#           the given type
#     callback:  optional function to be called at the end of 
#                each iteration

# matpow() value:
#     ev:  see below

# matpow() maintains an R environment ev, accessible to the callback
# function (if any), with these components:
#    stop:  a boolean indicating whether to end the iteration process,
#           e.g. to stop when convergence is reached
#    m:  the argument m in matpow()
#    k:  the argument k in matpow()
#    squaring: the argument squaring in matpow()
#    i:  current iteration number 
#    finaliter:  TRUE value means no more iterations
#    prod1: most recent odd-iteration power of m
#    prod2: most recent even-iteration power of m

# callback call form is callbackname(ev,init,...), with the arguments
# being:
#    ev: the ev maintained in matpow()
#    init:  the first time matpow() calls the callback, it sets
#           init to TRUE, to enable application-specific initialization;
#           at that time, the callback has the opportunity to add
#           new variables to ev

matpow <- function(m,k=NULL,squaring=FALSE,
      genmulcmd=NULL,dup=NULL,callback=NULL,...) {
   # general initialization of ev (number of iterations etc.)
   ev <- new.env()
   ev$stop <- FALSE
   ev$finaliter <- FALSE
   ev$squaring <- squaring
   ev$m <- m
   ev$k <- k
   if (!is.null(genmulcmd)) {
      ev$genmulcmd <- genmulcmd
   } else setgenmulcmd(ev)
   if (!is.null(dup)) {
      ev$dup <- dup
   } else setdup(ev)   
   ev$prod1 <- ev$dup(m)
   ev$prod2 <- ev$dup(m)  # dummy, just to set up correct class etc.
   if (!is.null(callback)) callback(ev,cbinit=TRUE)
   if (!squaring) {
      niters <- ev$k - 1
   } else  {
      niters <- ceiling(log2(ev$k))
   }
   
   # OK, start iterations
   if (ev$k == 1) return(ev)
   for (i in 1:niters) {
      ev$i <- i
      doiter(ev)  # perform this iteration
      if (i == niters) {
         ev$finaliter <- TRUE
         if (niters %% 2 == 1) ev$prod1 <- ev$prod2
      }
      if (!is.null(callback)) {
        callback(ev,...)
      	if (ev$stop) {
           return(ev)
        }
      }
   }
   ev
}

setdup <- function(ev) {
   if (class(ev$m) == "matrix") {
      ev$dup <- dup.vanilla 
      return()
   } else if (class(ev$m) == "big.matrix")  {
      ev$dup <- dup.bigmemory 
      return()
   }
   stop("no dup() available")
}

setgenmulcmd <- function(ev) {
   if (class(ev$m) == "matrix") {
      ev$genmulcmd <- genmulcmd.vanilla 
      return()
   } else if (class(ev$m) == "big.matrix")  {
      ev$genmulcmd <- genmulcmd.bigmemory 
      return()
   }
   stop("no genmulcmd() available")
}

dup.vanilla <- function(mat) mat

dup.bigmemory <- function(mat) {
   require(bigmemory)
   tmp <- bigmemory::big.matrix(nrow=nrow(mat), ncol=ncol(mat))
   tmp[,] <- mat[,]
}

doiter <- function(ev) {
   p1 <- "ev$prod1"
   p2 <- "ev$prod2"
   m <- "ev$m"
   if (!ev$squaring) {
      if (ev$i %% 2 == 1)
         eval(parse(text=ev$genmulcmd(m,p1,p2))) else
         eval(parse(text=ev$genmulcmd(m,p2,p1))) 
   } else {
      if (ev$i %% 2 == 1)
         eval(parse(text=ev$genmulcmd(p1,p1,p2))) else
         eval(parse(text=ev$genmulcmd(p2,p2,p1))) 
   }
}

genmulcmd.vanilla <- function(a,b,c) {
   paste(c," <- ",a," %*% ",b)  
}

genmulcmd.bigmemory <- function(a,b,c) {
   paste(c,"[,] <- ",a,"[,] %*% ",b,"[,]")  
}

genmulcmd.gputools <- function(a,b,c) {
   paste(c," <- gpuMatMult(",a,",",b,")")
}

# norm of vector x
normvec <- function(x) sqrt(sum(x^2))

