#' Find a Condorcet Choice.
#'
#' Try to compute find a Condorcet choice given a full ranking of m objects.
#'
#' In a ranking of m objects, the Condorcet choice is the choice that wins over every other choice in pairwise comparisons.  See Marden (1995), p.20 for details.
#' 
#' @param data the data, a vector of counts of each permutation of the m objects (m is the length of data)
#' @param names character vector of the names of the m objects
#' @return ...
#' @seealso \code{\link{Pmaker}}
#' @export condorcet
#' @references Marden, J. I. (1995). \emph{Analyzing and Modeling Rank Data}, London: Chapman & Hall. p.20.
#' @examples
#' 
#' data(city)
#'
#' condorcet(city[,"city"], colnames(city))    # among city-dwellers
#' condorcet(city[,"suburb"], colnames(city))  # among suburb-dwellers
#' condorcet(city[,"country"], colnames(city)) # among country-dwellers
#' condorcet(rowSums(city), colnames(city))    # overall winner
#'  
#'
condorcet <- function(data, names){

  # determine m
  mFac <- length(data)
  boole <- TRUE
  m <- 0
  while(boole){
    m <- m + 1  	
    if(mFac == factorial(m)) boole <- FALSE
    if(m == 50) stop("condorcet size failsafe; the length of data is not factorial.", call. = FALSE)
  }
  
  # compute Khat (percentage scale)
  n <- sum(data)
  KhatVec <- Pmaker(m) %*% data
  pctKhat <- (upper(KhatVec) + lower(n-KhatVec))/n
  pctKhat[which(diag(m) == 1)] <- NA
  K2 <- t(apply(pctKhat, 1, function(v) v[!is.na(v)])) # remove diag elts
  winner <- apply(K2, 1, function(x) all(x > .5))
  
  if(any(winner)){
    if(!missing(names)) return(names[which(winner)])
    return(which(winner))
  }
  
  print("No Condorcet choice found.")
  invisible(pctKhat)
}


