#' @title beta-diversity of two networks
#' @description
#' measures the beta-diversity between two networks
#' @param n1 network 1 (as an igraph object)
#' @param n2 network 2 (as an igraph object)
#' @param bf any function to measure beta-diversity between two sets
#'
#' @return a list with components S, OS, WN, and ST. While interpreting
#' the output, it is important to consider that ST is strongly constrained by
#' the values of S (the species composition dissimilarity). ST is only really
#' meaningful when the values of S are "intermediate"; a good example is when
#' the networks have been sampled along a gradient, and a more or less equal
#' proportion of the species show turnover from one step to the next. In the
#' situations where S is either really high or really low, the values of ST
#' are constrained and should no be given importance. The values of OS and WN,
#' and how they relate to S, have more informative value.
#' @export
betalink <- function(n1,n2,bf=B01){
   # Vertices in the two networks
   v1 <- igraph::V(n1)$name
   v2 <- igraph::V(n2)$name
   vs <- v1[v1 %in% v2] # Shared vertices
   beta_S <- bf(betapart(v1, v2))
   # Why can't igraph just expose the name of edges? WHY?
   # This is fugly
   # I hate this bullshit
   e1 <- plyr::aaply(igraph::get.edgelist(n1), 1, function(x) stringr::str_c(x, collapse='--', paste='_'))
   e2 <- plyr::aaply(igraph::get.edgelist(n2), 1, function(x) stringr::str_c(x, collapse='--', paste='_'))
   beta_WN <- bf(betapart(e1, e2))
   if(length(vs)>=2)
   {
      sn1 <- igraph::induced.subgraph(n1, which(igraph::V(n1)$name %in% vs))
      sn2 <- igraph::induced.subgraph(n2, which(igraph::V(n2)$name %in% vs))
      se1 <- plyr::aaply(igraph::get.edgelist(sn1), 1, function(x) stringr::str_c(x, collapse='--', paste='_'))
      se2 <- plyr::aaply(igraph::get.edgelist(sn2), 1, function(x) stringr::str_c(x, collapse='--', paste='_'))
      beta_OS <- bf(betapart(se1, se2))
      beta_ST <- beta_WN - beta_OS
   } else {
      beta_OS <- NaN
      beta_ST <- NaN
   }
	return(list(S = beta_S, OS = beta_OS, WN = beta_WN, ST = beta_ST))
}

#' @title Partition sets A and B
#' @description
#' given any two sets (arrays) A and B, return the size of components
#' a, b, and c, used in functions to measure beta-diversity
#' @param A any array
#' @param B any array
#' @export
#' @examples
#' A = c(1,2,3)
#' B = c(2,3,4)
#' betapart(A, B)
betapart <- function(A,B) list(b=sum(!(A %in% B)), c=sum(!(B %in% A)), a=sum(B %in% A))


#' @title Give names to networks
#' @description
#' If the networks (in a list) have no names, give them names
#' @param w A list (of networks, but who am I to judge?)
#' @export
name_networks <- function(w){
   if(is.null(names(w))){
      warning("It is recommended to give names to your networks. I've done it for you.")
      names(w) <- paste("network", c(1:length(w)), sep='_')
   }
   return(w)
}

#' @title Anemone/fish interaction data
#' @docType data
#' @keywords dataset
#' @name anemonefish
#' @format 16 adjacency matrices with species names
#' @description From http://mangal.io/data/dataset/2/
NULL
