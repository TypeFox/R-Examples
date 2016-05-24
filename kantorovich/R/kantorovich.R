#' @importFrom gmp as.bigq
#'
asab <- function(x) as.character(gmp::as.bigq(x))

#' Names for bigq vectors
#'
#' @param x a \code{bigq} vector
#' @return the names of \code{x}
#'
#' @export
names.bigq <- function(x){
    attr(x, "names")[1:length(x)]
}

# #' Extract vector/matrix elements by names
# #'
# #'
# `[.bigq` <- function (x, i = NULL, j = NULL, drop = TRUE)
# {
#   if(is.character(i)){ i <- sapply(i, function(k) which(names(x)==k), USE.NAMES=FALSE) }
#   .Call(gmp:::matrix_get_at_q, x, i, j)
# }


#' @importFrom gmp as.bigq apply
#' @importFrom methods formalArgs
#'
Vectorize_bigq <- function(f){
  if(length(formalArgs(f)) != 2) stop("Intended only for two variables function")
    return(function(x,y) gmp::as.bigq(gmp::apply(cbind(x,y), 1,
                  FUN=function(row) as.character(f(row[1], row[2])))))
}


#' @importFrom gmp as.bigq
#' @importFrom stats setNames
#'
arrange_names <- function(mu, nu){
  if(is.character(mu)) mu <- setNames(as.bigq(mu), names(mu))
  if(is.character(nu)) nu <- setNames(as.bigq(nu), names(nu))
  # check sum ==1
  if(sum(mu) != 1) stop("sum(mu) != 1")
  if(sum(nu) != 1) stop("sum(nu) != 1")
  # check class
  if(!all(c(class(mu), class(nu)) %in% c("integer", "numeric", "bigq"))) stop("mu and nu must be vectors in numeric or bigq/character class")
  if(class(mu) != class(nu)) stop("mu and nu have not the same class")
  #
  if(is.null(names(mu)) && is.null(names(nu)) && length(mu)==length(nu)){
    names(mu) <- seq_along(mu); names(nu) <- seq_along(nu)
    return(list(mu=mu, nu=nu))
  }
  #
  names_mu <- names(mu); names_nu <- names(nu)
  has_changed <- function(x, y) !identical(names(x), names_mu) || !identical(names(y), names_nu)
  if(is.null(names(mu))) names(mu) <- seq_along(mu)
  if(is.null(names(nu))) names(nu) <- seq_along(nu)
  if(setequal(names(mu), names(nu))){
    if(has_changed(mu,nu)) message("Caution: some names of mu and/or nu were missing or not compatible - automatic change")
    return(list(mu=mu, nu=nu))
  } else if(length(mu) == length(nu)){
    stop("Cannot deal with the names of mu and nu")
  } else if(length(mu) < length(nu)){
    if(all(names(mu) %in% names(nu))){
      message("Caution: some names of mu and/or nu were missing or not compatible - automatic change")
      if(class(mu)=="bigq"){
        mu_ch <- setNames(as.character(mu), names(mu))
        mu_ch[setdiff(names(nu), names(mu))] <- "0"
        mu <- setNames(gmp::as.bigq(mu_ch), names(mu_ch))
      } else {
        mu[setdiff(names(nu), names(mu))] <- 0
      }
      return(arrange_names(mu,nu))
    }else{
      stop("Cannot deal with the names of mu and nu")
    }
  }else if(length(mu) > length(nu)){
    return(setNames(arrange_names(nu,mu)[2:1], c("mu", "nu")))
  }else{
    stop("Cannot deal with the names of mu and nu")
  }
}


#' @importFrom gmp as.bigq
discrete <- function(x, y, gmp=FALSE){
  out <- if(gmp) gmp::as.bigq(x != y) else as.integer(x != y)
  return(out)
}

#' Extreme joinings
#'
#' Return extreme joinings between \code{mu} and \code{nu}.
#'
#' @param mu (row margins) probability measure in numeric or bigq/character mode
#' @param nu (column margins) probability measure in numeric or bigq/character mode
#' @param zeros logical; in case when \code{mu} and \code{nu} have differente lengths, set \code{FALSE} to remove lines or columns full of zeros
#'
#' @return A list containing the extreme joinings (matrices).
#'@examples
#' mu <- nu <- c(0.5, 0.5)
#' ejoinings(mu, nu)
#' # use exact arithmetic
#' library(gmp)
#' mu <- nu <- as.bigq(c(0.5,0.5))
#' ejoinings(mu, nu)
#' # different lengths example
#' mu <- setNames(as.bigq(c(1,2,4), 7), c("a", "b", "c"))
#' nu <- setNames(as.bigq(c(3,1), 4), c("b", "c"))
#' ejoinings(mu, nu)
#'
#' @importFrom stats model.matrix
#' @export
ejoinings <- function(mu, nu, zeros=FALSE){
  mu0 <- mu; nu0 <- nu
  munu <- arrange_names(mu, nu)
  mu <- munu$mu; nu <- munu$nu
  if(class(mu) != class(nu)) stop("Enter mu and nu in numeric or (preferably) in rational with the gmp package.")
  if(class(mu) != "bigq") message("Message: You should enter mu and nu in rational with the gmp package.")
  if(length(mu)>1){
    if(length(nu)>1){
      B <- c(mu,nu)
    }else{ B <- mu }
  }else{ B <- nu }
  m <- length(mu)
  n <- length(nu)
  if(m>1){ M1 <- t(model.matrix(~0+gl(m,n)))[,] }else{ M1 <- NULL }
  if(n>1){ M2 <- t(model.matrix(~0+factor(rep(1:n,m))))[,] }else{ M2 <- NULL }
  M <- rbind(M1,M2)
  if(class(mu)=="bigq"){
    mH0 <- rcdd::makeH(a1=asab(-diag(m*n)), b1=asab(rep(0,m*n)), a2=asab(M), b2=asab(B))
  }else{
    mH0 <- rcdd::makeH(a1=-diag(m*n), b1=rep(0,m*n), a2=M, b2=B)
  }
  extremals <- rcdd::scdd(mH0)$output[,-c(1,2)]
  if(is.null(dim(extremals))) extremals <- matrix(extremals, nrow=1)
  extremals <- lapply(1:nrow(extremals), function(i) matrix(extremals[i,], ncol=n, byrow=TRUE) )
  out <- lapply(extremals,
         function(M){
           rownames(M) <- names(mu)
           colnames(M) <- names(nu)
           return(M[, sapply(names(mu), function(k) which(names(nu)==k), USE.NAMES=FALSE)])
         })
  if(!zeros && length(mu0) != length(nu0)){
    if(length(mu0) < length(nu0)){
      out <- lapply(out, function(joining){
        which.zeros <- sapply(seq_along(mu), function(i) all(joining[i,]==0) || all(joining[i,]=="0"))
        return(joining[!which.zeros,])
      })
    } else {
      out <- lapply(out, function(joining){
        which.zeros <- sapply(seq_along(nu), function(i) all(joining[,i]==0) || all(joining[,i]=="0"))
        return(joining[,!which.zeros])
      })
    }
  }
  return(out)
}

#' Extremal distances
#'
#' Compute the distances at the extreme joinings.
#'
#' @param mu (row margins) probability measure in numeric or bigq/character mode
#' @param nu (column margins) probability measure in numeric or bigq/character mode
#' @param dist function or matrix, the distance to be minimized on average. If \code{NULL}, the 0-1 distance is used.
#' @param ... arguments passed to \code{dist}
#'
#' @return A list with two components: the extreme joinings in a list and the distances in a vector.
#'
#' @note This function, called by \code{\link{kantorovich}}, is rather for internal purpose.
#'
#' @importFrom gmp as.bigq
#' @export
edistances <- function(mu, nu, dist=NULL, ...){
  joinings <- ejoinings(mu, nu, zeros=TRUE)
  n.joinings <- length(joinings)
  j1 <- joinings[[1]]
  use_gmp <- class(mu) %in% c("bigq", "character")
  if(is.null(dist)){
    rho <- function(x, y) discrete(x, y, gmp=use_gmp)
  } else if(class(dist) == "function") {
    rho <- function(x, y) dist(x, y, ...)
  } else if(class(dist)=="matrix"){
    if(!use_gmp && mode(dist) != "numeric") stop("The dist matrix must be numeric if mu and nu are numeric")
    if(nrow(dist) != length(mu) || ncol(dist) != length(nu)) stop("Invalid dimension of the dist matrix")
    if(is.null(rownames(dist))) rownames(dist) <- 1:nrow(dist)
    if(is.null(colnames(dist))) colnames(dist) <- 1:ncol(dist)
    if(!setequal(rownames(j1), rownames(dist)) || !setequal(colnames(j1), colnames(dist))) stop("Invalid dimension names of the dist matrix")
  } else {
    if(!use_gmp) stop("dist must be a function or a numeric matrix")
    if(use_gmp) stop("dist must be a function or a numeric/character matrix")
  }
  if(class(dist) == "matrix") {
    Rho <- dist[rownames(j1), colnames(j1)]
  } else {
    if(use_gmp){
      Rho <- outer(rownames(j1), colnames(j1), FUN=Vectorize_bigq(rho))
    } else {
      Rho <- outer(rownames(j1), colnames(j1), FUN=rho)
    }
  }
  distances <- if(use_gmp) gmp::as.bigq(numeric(n.joinings)) else numeric(n.joinings)
  for(k in 1:n.joinings){
    joining <- joinings[[k]]
    if(use_gmp){
      distances[k] <- sum(Rho * as.bigq(joining))
    } else {
      distances[k] <- sum(Rho * joining)
    }
  }
  out <- list(joinings=joinings, distances=distances)
  return(out)
}

#' Kantorovich distance
#'
#' Compute the Kantorovich distance between two probability measures on a finite set.
#'
#' @param mu (row margins) probability measure in numeric or bigq/character mode
#' @param nu (column margins) probability measure in numeric or bigq/character mode
#' @param dist function or matrix, the distance to be minimized on average; if \code{NULL}, the 0-1 distance is used.
#' @param details prints the joinings achieving the Kantorovich distance and returns them in the \code{"joinings"} attribute of the output
#' @param ... arguments passed to \code{dist} (only if it is a function)
#'
#' @return The Kantorovich distance between \code{mu} and \code{nu}.
#'
#' @examples
#' mu <- c(1/7, 2/7, 4/7)
#' nu <- c(1/4, 1/4, 1/2)
#' kantorovich(mu, nu)
#' library(gmp)
#' mu <- as.bigq(c(1,2,4), 7)
#' nu <- as.bigq(c(1,1,1), c(4,4,2))
#' kantorovich(mu, nu)
#' mu <- c("1/7", "2/7", "4/7")
#' nu <- c("1/4", "1/4", "1/2")
#' kantorovich(mu, nu, details=TRUE)
#'
#' @details The function firstly computes all the extreme joinings of \code{mu} and \code{nu}, then evaluates the average distance for each of them, and then returns the minimal one.
#'
#' @export
kantorovich <- function(mu, nu, dist=NULL, details=FALSE, ...){
  distances <- edistances(mu=mu, nu=nu, dist=dist, ...)
  best <- which(distances$distances==min(distances$distances))
  kanto <- distances$distances[[best[1]]]
  if(details){
    joinings <- distances$joinings
    njoinings <- length(joinings)
    bestjoinings <- joinings[best]
    message1 <- sprintf("The Kantorovich distance is achieved for %s joining(s) among the %s extreme joining(s), given in the 'joinings' attribute of the output.\n", length(best), njoinings)
    cat(message1)
    attr(kanto, "joinings") <- bestjoinings
  }
  return(kanto)
}
