#' Update the Grouping Matrix
#' 
#' Updates the grouping matrix. 
#' 
#' @param G Grouping matrix.
#' @param a Current value of the row scores.
#' @param bwts2 Squared column weights.
#' @param Fr.bk Product of Fr and bkmat.
#' @param const Constant part of the loss function.
#' @param n The number of observations.
#' @param m The number of items.
#' @param q The number of rating categories.
#' @param random Logical indicating whether to randomize the observations.
#' @param info.level Integer controlling the amount of printed.
#' @author Pieter Schoonees
#' @keywords multivariate
#' @export updateG
updateG <- function(G, a, bwts2, Fr.bk, const, n, m, q, random = FALSE, info.level = 3) {
  time.start <- proc.time()[3]
  if(nrow(G) != n) stop("G: incorrect number of rows")
  G.upd <- G
  K <- ncol(G)
ind <- if(random) (1:n)[sample(1:n, n)] else 1:n
  crit.start <- Lfun.G.upd(G = G, a.cur = a, bwts2 = bwts2, Fr.bk = Fr.bk, 
                            n = n, m = m, q = q, K = K)$out
  nr.moves <- 0
  
  for(i in ind) {
      crit <- 0.5*(m+q-2)*bwts2*a[i]^2 - (a[i]*Fr.bk[i,] + a[i+n]*Fr.bk[i+n,])
      which.G <- which.min(crit)
      if(which.G != (1:K)[G[i,]==1]) {
          nr.moves <- nr.moves + 1
          G.upd[i,] <- 0
        	G.upd[i,which.G] <- 1
          }
    }
  
  nr.grp <- apply(G.upd,2,sum)
  K.new <- sum(nr.grp>0)
  which.pos <- nr.grp>0
  G.upd <- G.upd[,which.pos]
  if(K.new==1) G.upd <- matrix(G.upd, ncol = 1)
  time.end <- proc.time()[3]
  crit.end <- Lfun.G.upd(G = G.upd, a.cur = a, bwts2 = bwts2[which.pos], Fr.bk = Fr.bk[,which.pos], 
                               n = n, m = m, q = q, K = K.new)
  if(crit.start < crit.end$out) warning("G: Update not monotone decreasing")

  if(info.level == 3 || info.level == 4) 
    cat("\n\t\t G loss =", (m+q-2)*crit.end$out + const, "\t | \t",
          round(time.end - time.start,2),"sec's \t | \t moves", 
          nr.moves, "\t | \t improve", round((m+q-2)*crit.start - (m+q-2)*crit.end$out, 1), "\n \n")
  list(G = G.upd, K = K.new, classes.left = which.pos, moves = nr.moves, kloss = crit.end$kloss)
  }
