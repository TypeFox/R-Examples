#' ACSS complexity
#' 
#' Functions to obtain the algorithmic complexity for short strings, an approximation of the Kolmogorov Complexity of a short string using the coding theorem method.
#' 
#' @usage acss(string, alphabet = 9)
#' 
#' local_complexity(string, alphabet = 9, span = 5)
#' 
#' likelihood_d(string, alphabet = 9)
#' 
#' likelihood_ratio(string, alphabet = 9)
#' 
#' prob_random(string, alphabet = 9, prior= 0.5)
#' 
#' @param string \code{character} vector containing the to be analyzed strings (can contain multiple strings).
#' @param alphabet \code{numeric}, the number of possible symbols (not necessarily actually appearing in str). Must be one of \code{c(2, 4, 5, 6, 9)} (can also be \code{NULL} or contain multiple values for \code{acss()}). Default is 9.
#' @param prior \code{numeric},  the prior probability that the underlying process is random.
#' @param span size of substrings to be created from \code{string}.
#' 
#' @return
#' \describe{
#'   \item{"acss"}{A matrix in which the rows correspond to the strings entered and the columns to the algorithmic complexity K and the algorithmic probability D of the string (see \url{http://complexitycalculator.com/methodology.html}).}
#'   \item{"local_complexity"}{A list with elements corresponding to the strings. Each list containes a named vector of algorithmic complexities (K) of all substrings in each string with length span.}
#'   \item{"likelihood_d"}{A named vector with the likelihoods for \code{string} given a detreministic process.}
#'   \item{"likelihood_ratio"}{A named vector with the likelihood ratios (or Bayes factors) for \code{string} given a random rather than detreministic process.}
#'   \item{"prob_random"}{A named vector with the posterior probabilities that for a random process given the strings and the provided prior for being produced by a random process (default is 0.5, which correspond to a prior of 1 - 0.5 = 0.5 for a detereministic process).}
#'   }
#' 
#' @details The algorithmic complexity is computed using the coding theorem method: For a given alphabet size (number of different symbols in a string), all possible or a large number of random samples of Turing machines (TM) with a given number of states (e.g., 5) and number of symbols corresponding to the alphabet size were simulated until they reached a halting state or failed to end.
#' The outputs of the TMs at the halting states produces a distribution of strings known as the algorithmic probability of the strings. The algorithmic coding theorem (Levin, 1974) establishes the connection between the complexity of a string \eqn{K(s)} and its algorithmic probability \eqn{D(s)} as:
#' \deqn{K(s)\approx -\log_{2}(D(s))}{K(s) = log(D(s))}
#' 
#' This package accesses a database containing data on 4.5 million strings from length 1 to 12 simulated on TMs with 2, 4, 5, 6, and 9 symbols. 
#' 
#' For a more detailed discussion see Gauvrit, Singmann, Soler-Toscano, and Zenil (2014), \url{http://complexitycalculator.com/methodology.html}, or references below.
#' 
#' @references Delahaye, J.-P., & Zenil, H. (2012). Numerical evaluation of algorithmic complexity for short strings: A glance into the innermost structure of randomness. \emph{Applied Mathematics and Computation}, 219(1), 63-77. doi:10.1016/j.amc.2011.10.006 
#' 
#' Gauvrit, N., Singmann, H., Soler-Toscano, F., & Zenil, H. (2014). Algorithmic complexity for psychology: A user-friendly implementation of the coding theorem method. arXiv:1409.4080 [cs, stat]. \url{http://arxiv.org/abs/1409.4080}.
#' 
#' Gauvrit, N., Zenil, H., Delahaye, J.-P., & Soler-Toscano, F. (2014). Algorithmic complexity for short binary strings applied to psychology: a primer. \emph{Behavior Research Methods}, 46(3), 732-744. doi:10.3758/s13428-013-0416-0
#' 
#' Levin, L. A. (1974). Laws of information conservation (nongrowth) and aspects of the foundation of probability theory. \emph{Problemy Peredachi Informatsii}, 10(3), 30-35.
#' 
#' Soler-Toscano, F., Zenil, H., Delahaye, J.-P., & Gauvrit, N. (2012). Calculating Kolmogorov Complexity from the Output Frequency Distributions of Small Turing Machines. \emph{PLoS ONE}, 9(5): e96223.
#' 
#' @note The first time per session one of the functions described here is used, a relatively large dataset is loaded into memory which can take a considerable amount of time (> 10 seconds).
#' 
#' @example examples/examples.acss.R
#' 
#' @name acss
#' @aliases acss prob_random local_complexity likelihood_d likelihood_ratio
#' @export acss prob_random local_complexity likelihood_d likelihood_ratio
#' @importFrom zoo rollapply
#' @import acss.data
#' 


acss <- function(string, alphabet = 9) { #, return = "matrix") {
  check_string(string)
#  return <- match.arg(return, c("matrix", "data.frame"))
  names <- string
  string <- normalize_string(string)
  if (is.null(alphabet)) tmp <- acss_data[string,]  
  else {
    alphabet <- as.numeric(alphabet)
    if (any(!(alphabet %in% c(2, 4, 5, 6, 9)))) stop("alphabet must be in c(2, 4, 5, 6, 9)")
    tmp <- acss_data[string, paste("K", alphabet , sep = "."), drop = FALSE]
  }
  D <- apply(tmp, c(1,2), function(x) 2^(-x))
  colnames(D) <- paste0("D.", substr(colnames(D), 3, 3))  
#  if (return == "matrix") {
    tmp <- as.matrix(cbind(tmp, D))  
    rownames(tmp) <- names
    return(tmp)
#   } else if (return == "data.frame") {
#     tmp <- cbind(tmp, D)
#     rownames(tmp) <- make.unique(names)
#     return(tmp)
#   }
}


likelihood_d <- function(string, alphabet = 9) {
  if (length(alphabet) > 1) stop("'alphabet' needs to be of length 1.")
  alphabet <- as.numeric(alphabet)
  if (!(alphabet %in% c(2, 4, 5, 6, 9))) stop("alphabet must be in c(2, 4, 5, 6, 9)")
  check_string(string)
  l <- nchar(string)
  lu <- unique(l)
  rn <- nchar(rownames(acss_data))
  subtables <- lapply(lu, function(x) {
    tmp <- acss_data[rn == x, paste0("K.", alphabet), drop = FALSE]
    tmp <- tmp[!is.na(tmp[,paste0("K.", alphabet)]),,drop = FALSE]
    tmp$count <- count_class(rownames(tmp), alphabet = alphabet)
    tmp$D <- 2^(-tmp[,paste0("K.", alphabet)])
    tmp
  })
  #browser()
  ptot <- vapply(subtables, function(x) sum(x$count*x$D), 0)
  psgiventm <- acss(string, alphabet = alphabet)[,paste0("D.", alphabet)]/ptot[match(l, lu)]
  names(psgiventm) <- string
  psgiventm  
}

likelihood_ratio <- function(string, alphabet = 9) {
  llk_deterministic <- likelihood_d(string, alphabet = alphabet)
  l <- nchar(string)
  lu <- unique(l)
  llk_random <- (1/alphabet^lu)[match(l, lu)]
  llk_random/llk_deterministic
}

# old version:
prob_random <- function(string, alphabet = 9, prior= 0.5){
  l <- nchar(string)
  lu <- unique(l)
  psgivenr <- (1/alphabet^lu)[match(l, lu)]
  psgiventm <- likelihood_d(string, alphabet = alphabet)
  psgivenr*prior/(psgivenr*prior+psgiventm*(1-prior))
#   tmp <- psgivenr*prior/(psgivenr*prior+psgiventm*(1-prior))
#   names(tmp) <- string
#   tmp
}

local_complexity <- function(string, alphabet = 9, span = 5) {
  check_string(string)
  if (length(alphabet) > 1) stop("'alphabet' needs to be of length 1.")
  alphabet <- as.numeric(alphabet)
  if (!(alphabet %in% c(2, 4, 5, 6, 9))) stop("alphabet must be in c(2, 4, 5, 6, 9)")
  if (span < 2 | span > 12) stop("span needs to be between 2 and 12 (inclusive).")
  #browser()
  #l <- nchar(string)
  splitted <- strsplit(string,"")  
  new.string <- lapply(splitted, function(x) rollapply(x, width = span, FUN = paste0, collapse = ""))  
  tmp <- lapply(new.string, function(x) acss(x, alphabet = alphabet)[,paste0("K.", alphabet)])
  tmp <- mapply(function(x,y) {
    names(x) <- y
    return(x)}, tmp, new.string, SIMPLIFY = FALSE)
  #names(tmp) <- make.unique(string)
  names(tmp) <- string
  tmp
}
