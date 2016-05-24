#' Calculates Barber's Bipartite Modularity
#'
#' Modularity, community formation, boundary clumping. Call it what you want,
#' there is a fair amount of overlap in definition here. As such, we posit that
#' modularity statistics may be able to detect boundary clumping better than
#' Morisita's index. We offer a function here to calculate modularity.
#' Specifically, this function calculates Barber's Q statistic, and compares it
#' relative to null model randomizations.
#'
#' \code{method} is an argument handed to functions in the \code{vegan}
#' package. Leibold & Mikkelson advocated the use of equiprobable rows and
#' columns (provided that rows and columns had at least one entry). This method
#' is called \code{r00}. Methods maintaining row (site) frequencies include
#' \code{r0},\code{r1}, and \code{r2}. The default method argument is
#' \code{r1}, which maintains the species richness of a site (row totals) and
#' fills species ranges (columns) based on their marginal probabilities.
#' Arguably the most conservative null algorithm is the fixed row - fixed
#' column total null, which can be attained using many of swap algorithms
#' described in the vegan package (sequential methods like \code{tswap},
#' \code{swap}, and non-sequential \code{quasiswap} and \code{backtracking}).
#' See the help file for \code{commsim} or Wright et al. 1998 for more
#' information.
#'
#' If \code{order} is FALSE, the interaction matrix is not ordinated, allowing
#' the user to order the matrix based on site characteristics or other
#' biologically relevant characteristics.
#'
#' @param comm community data in the form of a presence absence matrix
#' @param method null model randomization method used by \code{NullMaker}. See
#' details below (and the help file of fucntion \code{NullMaker}) for more
#' information.
#' @param sims number of simulated null matrices to use in analysis
#' @param scores axis scores to ordinate matrix. 1: primary axis scores
#' (default) 2: secondary axis scores
#' @param order logical argument indicating whether to ordinate the interaction
#' matrix or not. See details.
#' @param c starting guess for the number of modules present. Defaults to the
#' maximum number of modules possible.
#' @param nstarts number of starts. Default is 100. More will both slow the
#' function down, and increase the likelihood of converging on the true
#' modularity value.
#' @return  A vector containing Barber's modularity statistic (Q), the z statistic comparing observed modularity against null matrices (z), p-value (pval), and mean (simulatedMean) and variance (simulatedVariance) from null model simulations
#'
#' @author Tomlin Pulliam and Tad Dallas
#' @export
#' @references
#'
#' Barber, M. J. 2007. Modularity and community detection in bipartite
#' networks. Physical Review E, 76(6), 066102.
#'
#' Leibold, M.A. and G.M. Mikkelson. 2002. Coherence, species turnover, and
#' boundary clumping: elements of meta-community structure. Oikos 97: 237 -
#' 250.
#'
#' Wright, D.H., Patterson, B.D., Mikkelson, G.M., Cutler, A. & Atmar, W.
#' (1998). A comparative analysis of nested subset patterns of species
#' composition. Oecologia 113, 1-20.
#' @examples
#'
#' #define an interaction matrix
#' data(TestMatrices)
#' intmat <- TestMatrices[[3]]
#'
#' #determine Barber"s modularity
#' mod.intmat <- Modularity(intmat, method="r1", sims=5,
#'  scores=1, order=TRUE, nstarts=10)
#'
#' #return results
#' mod.intmat
#'
Modularity <- function(comm, method='tswap', sims=1000, scores=1,  order=TRUE, c = length(comm), nstarts=100){

 #suport functions
 reduce <- function(v){
   n <- length(unique(v))
   while(any(sort(unique(v)) != 1:n)){
     su <- sort(unique(v))
     missing <- which(su != 1:n)[1]
     high <- which( v > missing)
     v[high] <- v[high] - 1
   }
   ret <- rep(0, n)
   for(i in 1:n){
     swap <- which(unique(v)==i)[1]
     ret[which(v==i)] <- swap
   }
   return(ret)
 }

 BarberQ <- function(comm, modules){
   modules <- reduce(modules)
   p <- nrow(comm)	# Number of red nodes
   q <- ncol(comm)	# Number of blue nodes
   n <- p+q		# Number of nodes
   if(n != length(modules)){
     warning("Please ensure that every node is assigned to exactly one module.")
   }
   m <- sum(comm != 0)	# Number of links
   c <- max(modules)	# Number of modules

   # Adjacency Matrix
   A <- rbind(cbind(matrix(0, p, p), comm), cbind(t(comm), matrix(0, q, q)))
   # Probabilities of potential links
   Ptilde <- (rowSums(comm) %*% t(colSums(comm)))/m
   # Probabilities of all links
   P <- rbind(cbind(matrix(0, p, p), Ptilde), cbind(t(Ptilde), matrix(0,q,q)))
   # Modularity Matrix
   B <- A - P
   # Index matrix
   S <- matrix(0, n, c)
   for(i in 1:n){
     S[i,modules[i]] <- 1
   }
   # Barber's Bipartite Modularity
   Q <- sum(diag(t(S) %*% B %*% S))/(2*m)
   return(Q)
 }


 maximizeQ <- function(comm, modules, c=sum(dim(comm))){
   modules <- reduce(modules)
   p <- nrow(comm)
   q <- ncol(comm)
   n <- p+q
   m <- sum(comm != 0)

   Ptilde <- (rowSums(comm) %*% t(colSums(comm)))/m
   Btilde <- comm - Ptilde

   S <- matrix(0, n, c)
   for(i in 1:n){
     S[i,modules[i]] <- 1
   }

   R <- S[1:p,]
   T <- S[(p+1):n,]

   lastMods <- rep(0, n)

   while(!all(lastMods == modules)){
     lastMods <- modules
     Ttilde <- Btilde %*% T
     redMods <- apply(Ttilde, 1, which.max)
     for(i in 1:p){
       tmp <- rep(0, c)
       tmp[redMods[i]] <- 1
       R[i,] <- tmp
     }
     Rtilde <- t(Btilde) %*% R
     blueMods <- apply(Rtilde, 1, which.max)
     for(i in 1:q){
       tmp <- rep(0, c)
       tmp[blueMods[i]] <- 1
       T[i,] <- tmp
     }
     modules <- c(redMods, blueMods)
   }
   modules <- reduce(modules)
   return(modules)
 }

 #actual modularity calculation
 getQ <- function(comm){
   n <- sum(dim(comm))
   best.modules <- rep(1,n)
   best.modules <- maximizeQ(comm, best.modules, c)
   best.Q <- BarberQ(comm, best.modules)

   for(i in 1:nstarts){
     modules <- reduce(sample(c, n, replace=TRUE))
     modules <- maximizeQ(comm, modules, c)
     Q <- BarberQ(comm, modules)
     if(Q > best.Q){
       best.modules <- modules
       best.Q <- Q
     }
   }
   ret <- list(Q=best.Q)
   return(ret)
 }

	nulls <- NullMaker(comm=comm, sims=sims, method=method, scores=scores, ordinate=order)
  Qemp <- getQ(comm)$Q
	Qnull <- as.numeric(unlist(lapply(nulls, getQ)))
  varstat <- sd(Qnull)
  z <- (mean(Qnull) - Qemp) / (varstat)
  pval <- 2 * pnorm(-abs(z))
 return(c(Q = Qemp, z = z, pval = pval, simulatedMean = mean(Qnull), simulatedVariance = varstat))
}
