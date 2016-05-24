#' Create the expected higher-order statistics calculating matrix for approval data
#'
#' Create the expected higher-order statistics calculating matrix for approval data
#' 
#' @param m the number of objects
#' @param vin the (lower order) grouping level of the data
#' @param vout the desired higher order grouping level
#' @return ...
#' @seealso \code{\link{Tmaker}}, \code{\link{Amaker}}, \code{\link{Mmaker}}, \code{\link{Pmaker}}, \code{\link{Smaker}}
#' @export Emaker
#' @examples
#' 
#' 
#' Emaker(6, 0, 1)
#' Emaker(6, 0, 2)
#' Emaker(6, 0, 3)
#' Emaker(6, 0, 4)
#' 
#' Emaker(6, 1, 1)
#' Emaker(6, 1, 2)
#' Emaker(6, 1, 3)
#' Emaker(6, 1, 4)
#' Emaker(6, 1, 5)
#' Emaker(6, 1, 6)
#'
#' # compare to Tmaker
#' Emaker(6, 1, 3) # contributors when bumping up from 1-groups to 3-groups
#' Tmaker(6, 3, 1)
#' 
Emaker <- function(m, vin, vout){
  
  ssetsNK <- subsets(m, vin)
  ssetsNKpd <- subsets(m, vout)
  
  ul <- unlist(
    lapply(
      ssetsNKpd,
      function(y) sapply(ssetsNK, function(x) containsQ(y,x))
    )  
  )
  
  mat <- matrix(ul, nrow = choose(m, vout), ncol = choose(m, vin), byrow = TRUE)
  mat <- (mat+0L) 
  
  colnames(mat) <- sapply(ssetsNK, paste, collapse = "")
  row.names(mat) <- sapply(ssetsNKpd, paste, collapse = "")

  mat

}























#' Create the sufficient statistics calculating matrix for approval data
#'
#' Create the sufficient statistics calculating matrix for approval data
#' 
#' @param m the number of objects
#' @param k the number of objects selected
#' @param d the order-effect for the desired matrix (0 to k)
#' @return ...
#' @seealso \code{\link{Emaker}}, \code{\link{Amaker}}, \code{\link{Mmaker}}, \code{\link{Pmaker}}, \code{\link{Smaker}}
#' @export Tmaker
#' @examples
#' 
#' 
#' Tmaker(4, 2, 0) # m
#' Tmaker(4, 2, 1) # generates how many of each
#' Tmaker(4, 2, 2) # gives data (order = subsets(1:4, 2))
#' 
#' Tmaker(5, 2, 0)
#' Tmaker(5, 2, 1)
#' Tmaker(5, 2, 2)
#' 
#' Tmaker(4, 3, 0) # 
#' Tmaker(4, 3, 1) # subsets(1:4, 3), 1 is in 1, 2, and 3
#' Tmaker(4, 3, 2) # subsets(1:4, 2)
#' Tmaker(4, 3, 3)
#'
#'
#'
#'
#'
#'
#'
#'
#'
#'
#'
#'
#' data(cookie)
#'
#'
#' ## voting statistics at different levels
#' ############################################################
#' 
#' # projection onto V0: the number of people in survey
#' effectsOnV0 <- Tmaker(6, 3, 0) %*% cookie$freq 
#' colnames(effectsOnV0) <- "Total Votes"
#' effectsOnV0 # = sum(cookie$freq)
#' 
#'
#' # projection onto V1: the number of people voting for each cookie
#' effectsOnV1 <- Tmaker(6, 3, 1) %*% cookie$freq  
#' row.names(effectsOnV1) <- cookie$cookies
#' colnames(effectsOnV1) <- "Total Votes" 
#' effectsOnV1
#'
#'
#' # projection onto V2: the number of people voting for each cookie-pair
#' effectsOnV2 <- Tmaker(6, 3, 2) %*% cookie$freq   
#' row.names(effectsOnV2) <- sapply(subsets(cookie$cookies, 2), paste, collapse = ", ")
#' colnames(effectsOnV2) <- "Total Votes"
#' effectsOnV2 
#' 
#'
#' # projection onto V3: the number of people voting for each cookie-triple
#' effectsOnV3 <- Tmaker(6, 3, 3) %*% cookie$freq    
#' row.names(effectsOnV3) <- sapply(subsets(cookie$cookies, 3), paste, collapse = ", ")    
#' colnames(effectsOnV3) <- "Total Votes" 
#' effectsOnV3 # = t(t(cookie$freq)) = the (freq) data
#' 
#'
#'
#'
#'
#'
#'
#'
#'
#'
#'
#'
#'
Tmaker <- function(m, k, d){
  
  kSubs <- subsets(1:m, k)
  dSubs <- subsets(1:m, d)

  inclusion_list <- lapply(kSubs, function(kv){
    lapply(dSubs, function(dv){
      as.integer( all(dv %in% kv) )
    })
  })  
  
  matrix(unlist(inclusion_list), nrow = choose(m, d), ncol = choose(m, k))
  
}













#' Distance transitive matrix
#'
#' Compute the distance transitive matrix for a survey in which k objects are selected from a group of m
#' 
#' @param m the number of objects
#' @param k the number of objects selected
#' @return ...
#' @seealso \code{\link{Tmaker}}, \code{\link{Emaker}}, \code{\link{Mmaker}}, \code{\link{Pmaker}}, \code{\link{Smaker}}
#' @export Amaker
#' @examples
#' 
#' 
#' Amaker(4, 2)
#' 
Amaker <- function(m, k){
  mat <- Tmaker(m, k, k-1)
  mat <- t(mat) %*% mat - k * diag(rep(1, choose(m, k))) 
  matrix(as.integer(mat), nrow = nrow(mat))
}









#' Marginals matrix
#'
#' Compute the marginals matrix for a full ranking of m objects
#'
#' This is the transpose of the marginals matrix presented in Marden (1995).
#' 
#' @param m the number of objects
#' @return ...
#' @seealso \code{\link{Tmaker}}, \code{\link{Amaker}}, \code{\link{Emaker}}, \code{\link{Pmaker}}, \code{\link{Smaker}}
#' @export Mmaker
#' @references Marden, J. I. (1995). \emph{Analyzing and Modeling Rank Data}, London: Chapman & Hall. p.42.
#' @examples
#' 
#'
#' data(city)
#'
#' Mmaker(3)
#' Mmaker(3) %*% city
#' 
Mmaker <- function(m){
  perms <- permutations(1:m)
  pairs <- as.matrix(expand.grid(1:m, 1:m))[,2:1]
  
  mat <- apply(perms, 1, function(s){
    apply(pairs, 1, function(ij){
      s[ij[1]] == ij[2]
    }) + 0L
  })
  
  colnames(mat) <- apply(perms, 1, paste, collapse = "")
  
  rn <- str_dup("+",m)
  str_sub(rn, rep(1:m, each = m), rep(1:m, each = m)) <- 
    rep(1:m, m)
  row.names(mat) <- rn    

  mat
}









#' Pairs matrix
#'
#' Compute the pairs matrix for a full ranking of m objects
#'
#' This is the transpose of the pairs matrix presented in Marden (1995).
#' 
#' @param m the number of objects
#' @return ...
#' @seealso \code{\link{Tmaker}}, \code{\link{Amaker}}, \code{\link{Emaker}}, \code{\link{Mmaker}}, \code{\link{Smaker}}
#' @export Pmaker
#' @references Marden, J. I. (1995). \emph{Analyzing and Modeling Rank Data}, London: Chapman & Hall. p.42.
#' @examples
#' 
#' data(city)
#'
#' Pmaker(3)
#' Pmaker(3) %*% city
#' # 1 = city, 2 = suburb, 3 = country
#'
#' # looking just among city folk, generate the pairs matrix
#' city[,"city",drop=FALSE] # the data
#' m <- sum(city[,"city"])
#' k <- (Pmaker(3) %*% city)[,1]
#' Khat <- upper(k) + lower(m-k)
#' colnames(Khat) <- row.names(Khat) <- colnames(city)
#' Khat
#' round(Khat / m, 2) # % times row is rated over column
#'
#'
#' # worked out: city is voted over suburb in 123 , 132, and 231, equaling
#' 210 + 23 + 8   # = Khat[1,2]
#' # whereas suburb is rated over city in 213, 312, 321, equaling
#' 111 + 204 + 81 # = Khat[2,1]
#'
#'
#' # is there a condorcet choice?
#' 
#' p <- ncol(Khat)
#' Khat[which(diag(p) == 1)] <- NA
#' K2 <- t(apply(Khat, 1, function(v) v[!is.na(v)])) # remove diag elts
#' boole <- apply(K2/m, 1, function(x) all(x > .5))
#' if(any(boole)) names(boole)[which(boole)]
#' # suburb is a condorcet choice
#' 
Pmaker <- function(m){

  perms <- permutations(1:m)	
  combos <- matrix(unlist(subsets(m, 2)), ncol = 2, byrow = TRUE)
  
  mat <- apply(perms, 1, function(s){
    apply(combos, 1, function(ij){
      s[ij[1]] < s[ij[2]]
    }) + 0L
  })
  
  colnames(mat) <- apply(perms, 1, paste, collapse = "")
  row.names(mat) <- apply(combos, 1, paste, collapse = ">")

  mat
}
























#' Means matrix (rank data)
#'
#' Compute the means matrix for a full ranking of m objects
#'
#' This is the transpose of the means matrix presented in Marden (1995); it projects onto the means subspace of a collection of ranked data.  See the examples for how to compute the average rank.
#' 
#' @param m the number of objects
#' @return ...
#' @seealso \code{\link{Tmaker}}, \code{\link{Amaker}}, \code{\link{Emaker}}, \code{\link{Mmaker}}, \code{\link{Pmaker}}
#' @export Smaker
#' @references Marden, J. I. (1995). \emph{Analyzing and Modeling Rank Data}, London: Chapman & Hall. p.41.
#' @examples
#' 
#' data(city)
#' 
#' X <- permutations(3)
#' 
#' # the average rank can be computed without this function
#' normalize <- function(x) x / sum(x)
#' factorial(3) * apply(t(X) %*% city, 2, normalize)
#' # the dataset city is really like three datasets; they can be pooled back
#' # into one via:
#' rowSums(city)
#' factorial(3) * apply(t(X) %*% rowSums(city), 2, normalize)
#'
#'
#' # the means matrix is used to summarize the data to the means subspace
#' # which is the subspace of m! spanned by the columns of permutations(m)
#' # note that when we project onto that subspace, the projection has the
#' # same average rank vector :
#' Smaker(3) %*% city # the projections, table 2.8
#' factorial(3) * apply(t(X) %*% Smaker(3) %*% city, 2, normalize)
#' 
#' # the residuals can be computed by projecting onto the orthogonal complement
#' (diag(6) - Smaker(3)) %*% city # residuals
#'
#' 
#' apply(t(X) %*% city, 2, function(x) x / sum(x) * factorial(3)) # average ranks by group
#' 
#' apply(t(X) %*% rowSums(city), 2, function(x) x / sum(x) * factorial(3)) # average ranks pooled
#' 
Smaker <- function(m){
  X <- permutations(1:m)	  
  zapsmall(X %*% solve(t(X) %*% X) %*% t(X))
}







































#' U matrix (rank data)
#'
#' Compute the generalized marginals matrix for a full ranking of m objects.  Umaker generalized Mmaker.
#'
#' This is the transpose of the generalized marginals matrix presented in Marden (1995).
#' 
#' @param m the number of objects
#' @return ...
#' @seealso \code{\link{Mmaker}}, \code{\link{Pmaker}}, \code{\link{Smaker}}
#' @export Smaker
#' @references Marden, J. I. (1995). \emph{Analyzing and Modeling Rank Data}, London: Chapman & Hall. pp.47--48.
#' @examples
#' 
#'
#' data(politicalGoals)
#' 
#' lambdas <- apply(partitions(4), 1, function(v) v[v != 0])
#' 
#' 
#' 
#' 
Umaker <- function(m){
  "this isn't finished"
}


