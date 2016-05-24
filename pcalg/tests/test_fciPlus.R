library(pcalg)

source(system.file(package="Matrix", "test-tools-1.R", mustWork=TRUE))
##--> showProc.time(), assertError(), relErrV(), ...
cat("doExtras:", (doExtras <- pcalg:::doExtras()), "\n")

##################################################
## simulate graph; remove confounders
##################################################

##' (utility)
##' for (i in 1:p),
##'   vect[i] := new name of node i in the graph when series nodes are removed
conf2vec <- function(confs, p) {
  stopifnot(p >= 1) # confs may be *empty*
  vect <- sapply(1:p, function(r) r - sum(1 <= confs & confs <= r))
  vect[confs] <- NA_integer_
  vect
}


##' random DAG  with confounders removed
rDAG_noC <- function(p, prob, max.no.conf = 3,
                     V = c(LETTERS,letters)[1:p], verbose=TRUE) {
  stopifnot(p == as.integer(p), length(p) >= 1, p >= 2, prob > 0,
            !anyDuplicated(V))

  g <- randomDAG(p, prob, V=V)  ## true DAG
  ## get the adjacency matrix
  A <- as(g, "matrix")
  A[A != 0] <- 1

  ## first three confounders should be removed
  ## confounders - vector of confounding nodes to be removed
  n <- ncol(A) ## number of nodes
  confounders <- numeric(max.no.conf)
  counter <- 1L
  for(k in 1:n) {
    if (sum(A[k,]) > 1) { ## k has more than one child
      confounders[counter] <- k ## mark node k; leave out later
      counter <- counter +1L
      if(counter > max.no.conf)
        break
    }
  }
  confounders <- confounders[is.c <- confounders > 0]
  if(verbose) {
    cat("confounders:\n"); print(confounders)
    vect <- conf2vec(confounders, p)
    cat("-> nodes after removing confounders:\n"); print(vect)
  }

  ## define sufficient statistics (d-separation oracle)
  trueC <- trueCov(g)
  ## delete rows and columns belonging to confounder variables in confounders
  C.1 <- if(any(is.c)) trueC[-confounders, -confounders] else trueC
  list(g = g, trueC = trueC, trueCov.1 = C.1, confounders = confounders, p. = p - sum(is.c))
} ## {rDAG_noC}

chk_rDAG <- function(rD) {
  stopifnot(is.list(rD), c("trueC", "trueCov.1", "p.","confounders") %in% names(rD))
  p <- nrow(rD$trueC)
  p. <- rD$p.
  conf <- rD$confounders
  no.c <- sum(conf > 0)
  stopifnot(is(rD$g, "graph"))
  V <- rownames(rD$trueC)
  stopifnot(p. == p - no.c, dim(rD$trueCov.1) == c(p., p.),
            isSymmetric(rD$trueC),
            isSymmetric(rD$trueCov.1),
            identical(V, colnames(as(rD$g, "matrix"))),
            identical(if(no.c) V[-conf] else V, rownames(rD$trueCov.1)))
}

seed <- 547 ## Seed 547: Dsep Link wird benoetigt
set.seed(seed)

## will get slow if p > 14; these settings take already a few minutes
p <- 14
r1 <- rDAG_noC(p, prob = 0.2)
c(p = p, p. = r1$p.)
chk_rDAG(r1)
###

indepTest1 <- gaussCItest

## transform covariance matrix into a correlation matrix
p. <- r1$p.
suffStat1 <- list(C = cov2cor(r1$trueCov.1), n = 10^9)

## FCI+: neue Variante
showSys.time(
  fciplus.fit <- fciPlus(suffStat1, indepTest = indepTest1, alpha = 0.99, p = p.)
)

## Fit FCI
showSys.time(
  fci.fit <- fci(suffStat1, indepTest = indepTest1, alpha = 0.99, p = p.,
                 rules = rep(TRUE, 10), verbose = FALSE)
)

print.table(fciplus.fit@amat, zero.print = ".")
## MM{FIXME}: want to see node correct names !!

if(!all(fciplus.fit@amat == fci.fit@amat))
    stop("Test fciPlus() differing from fci(): PAG not found correctly?")
showProc.time()

if(doExtras) { ## more tests
 for(p in 2:20) { ##     p = 10
   cat("\n\np = ", p, ":\n~~~~~~\n")
   set.seed(p)
   for(pr in c(0.75, 1.5) / p) {
     seed <- 10*p + ceiling(2*p*pr)
     cat("\nprob = ", pr, ", seed = ", seed,":\n")
     set.seed(seed)
     rr <- rDAG_noC(p, prob = pr) # quite sparse
     chk_rDAG(rr)

     ## transform covariance matrix into a correlation matrix
     p. <- rr$p.
     suffStat1 <- list(C = cov2cor(rr$trueCov.1), n = 10^9)

     ## FCI+: neue Variante
     showSys.time(
       fciplus.fit <- fciPlus(suffStat1, indepTest = indepTest1, alpha = 0.99, p = p.)
     )

     ## Fit FCI
     showSys.time(
       fci.fit <- fci(suffStat1, indepTest = indepTest1, alpha = 0.99, p = p.,
                      rules = rep(TRUE, 10), verbose = FALSE)
     )

     cat("Adjacency matrix:\n")
     print.table(fciplus.fit@amat, zero.print = ".")
     ## MM{FIXME}: want to see node correct names !!

     if(!all(fciplus.fit@amat == fci.fit@amat))
       stop("Test fciPlus() differing from fci(): PAG not found correctly?")
     else cat("=================================================\n")
   }# for different 'pr'

 }# for a set of 'p'

} ## if (doExtras)
