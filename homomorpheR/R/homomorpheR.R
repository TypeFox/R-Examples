#' homomorpheR: Homomorphic computations in R
#'
#' \code{homomorpheR} is a start at a rudimentary package for
#' homomorphic computations in R. The goal is to collect homomorphic
#' encryption schemes in this package for privacy-preserving
#' distributed computations; for example, applications of the sort
#' immplemented in package \code{distcomp}.
#'
#' At the moment, only one scheme is implemented, the Paillier
#' scheme. The current implementation makes no pretense at efficiency
#' and also uses direct translations of other implementations,
#' particularly the one in Javascript.
#'
#' For a quick overview of the features, read the
#' \code{\link{homomorpheR}} vignette by running
#' \code{vignette("homomorpheR")}.
#'
#' @references \url{https://en.wikipedia.org/wiki/Homomorphic_encryption}
#' @references \url{https://mhe.github.io/jspaillier/}
#'
#' @importFrom gmp as.bigz
#'
#' @examples
#' keys <- PaillierKeyPair$new(1024) # Generate new key pair
#' encryptAndDecrypt <- function(x) keys$getPrivateKey()$decrypt(keys$pubkey$encrypt(x))
#' a <- gmp::as.bigz(1273849)
#' identical(a + 10L, encryptAndDecrypt(a+10L))
#' x <- lapply(1:100, function(x) random.bigz(nBits = 512))
#' edx <- lapply(x, encryptAndDecrypt)
#' identical(x, edx)
#' @docType package
#' @name homomorpheR
#'
NULL
ONE <- as.bigz(1L)
ZERO <- as.bigz(0L)

#' Return a random big number using the cryptographically secure random number generator
#' from in the \code{sodium} package.
#'
#' @param nBits, the number of bits, which must be a multiple of 8, is not checked for efficiency.
#' @importFrom gmp as.bigz
#' @importFrom sodium random
#' @export
random.bigz <- function(nBits) {
    nBytes <- nBits %/% 8L
    as.bigz(paste0(c("0x", random(n = nBytes)), collapse = ""))
}
