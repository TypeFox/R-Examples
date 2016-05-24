#' Construct a Paillier public key with the given modulus.
#' @docType class
#' @seealso \code{PaillierPrivateKey} which goes hand-in-hand with this object
#' @importFrom R6 R6Class
#' @importFrom gmp add.bigz
#' @importFrom gmp mod.bigz
#' @importFrom gmp mul.bigz
#' @importFrom gmp powm
#' @importFrom gmp as.bigz
#' @field bits the number of bits in the modulus
#' @field n the modulus
#' @field nSquared the square of the modulus
#' @field nPlusOne one more than the modulus
#' @section Methods:
#' \describe{
#'   \item{\code{PaillierPublicKey$new(bits, n)}}{Create a new public key with given \code{bits} and modulus
#'         \code{n}. It also precomputes a few values for more efficient computations}
#'   \item{\code{PaillierPublicKey$encrypt(m)}}{Encrypt a message. The value \code{m} should be less than
#'         the modulus, not checked}
#'   \item{\code{PaillierPublicKey$add(a, b)}}{Return the sum of two encrypted messages \code{a} and \code{b}}
#'   \item{\code{PaillierPublicKey$mult(a, b)}}{Return the product of two encrypted messages \code{a} and
#'              \code{b}}
#' }
#'
#' @export
#' @format An \code{\link{R6Class}} generator object

PaillierPublicKey <- R6Class("PaillierPublicKey",
                            private = list(
                                randomize = function(a) {
                                    repeat {
                                        r <- random.bigz(nBits = self$bits)
                                        ## make sure r <= n
                                        if (r < self$n)
                                            break
                                    }
                                    rn <- powm(r, self$n, self$nSquared)
                                    mod.bigz(mul.bigz(a, rn), self$nSquared)
                                }),
                            public = list (
                                ## fields
                                bits = NULL,
                                n = NULL,
                                nSquared = NULL,
                                nPlusOne = NULL,
                                ## methods
                                initialize = function(bits, n) {
                                    ## bits
                                    self$bits <- bits
                                    ## n
                                    self$n <- n
                                    ## n squared
                                    self$nSquared <- mul.bigz(n, n)
                                    ## n plus 1
                                    self$nPlusOne <- add.bigz(n, ONE)
                                },
                                encrypt = function(m) {
                                    private$randomize(
                                        mod.bigz(
                                            add.bigz(
                                                mul.bigz(self$n, m),
                                                ONE),
                                            self$nSquared))
                                },
                                add = function(a, b) {
                                    mod.bigz(mul.bigz(a, b), self$nSquared)
                                },
                                mult = function(a, b) {
                                    powm(a, b, self$nSquared)
                                })
                            )

#' Construct a Paillier private key with the given secret and a public key
#' @docType class
#' @seealso \code{PaillierPublicKey} which goes hand-in-hand with this object
#' @importFrom R6 R6Class
#' @importFrom gmp add.bigz
#' @importFrom gmp mod.bigz
#' @importFrom gmp mul.bigz
#' @importFrom gmp div.bigz
#' @importFrom gmp sub.bigz
#' @importFrom gmp powm
#' @importFrom gmp inv.bigz
#' @importFrom gmp as.bigz
#' @field pubkey the Paillier public key
#'
#' @section Methods:
#' \describe{
#'   \item{\code{PaillierPrivateKey$new(lambda, pubkey)}}{Create a new private key with given secret
#'         \code{lambda} and the public key}
#'   \item{\code{PaillierPrivateKey$getLambda()}}{Return the secret \code{lambda}}
#'   \item{\code{PaillierPrivateKey$decrypt(c)}}{Decrypt a message. The value \code{c} should be an
#'         encrypted value}
#' }
#'
#' @export
#' @format An \code{\link{R6Class}} generator object
PaillierPrivateKey <- R6Class("PaillierPrivateKey",
                              private = list(
                                  lambda = NA,
                                  ## cached value for decryption
                                  x = NA
                              ),
                              public = list(
                                  ## fields
                                  pubkey = NULL,
                                  ## methods
                                  initialize = function(lambda, pubkey) {
                                      ## lambda
                                      private$lambda <- lambda
                                      self$pubkey <- pubkey
                                      ## x (cached) for decryption
                                      private$x <- inv.bigz(
                                          div.bigz(
                                              sub.bigz(
                                                  powm(pubkey$nPlusOne, private$lambda, pubkey$nSquared),
                                                  ONE),
                                              pubkey$n),
                                          pubkey$n)
                                  },
                                  getLambda = function() {
                                      private$lambda
                                  },
                                  decrypt = function(c) {
                                      mod.bigz(
                                          mul.bigz(
                                              div.bigz(
                                                  sub.bigz(
                                                      powm(c, private$lambda, self$pubkey$nSquared),
                                                      ONE),
                                                  self$pubkey$n),
                                              private$x),
                                          self$pubkey$n)
                                  })
                              )

#' Construct a Paillier public and private key pair given a fixed number of bits
#' @docType class
#' @seealso \code{PaillierPublicKey} and \code{PaillierPrivateKey}
#' @importFrom R6 R6Class
#' @importFrom gmp add.bigz
#' @importFrom gmp sub.bigz
#' @importFrom gmp mod.bigz
#' @importFrom gmp mul.bigz
#' @importFrom gmp powm
#' @importFrom gmp isprime
#' @importFrom gmp sizeinbase
#' @importFrom gmp lcm.bigz
#' @importFrom gmp as.bigz
#' @field pubkey the Paillier public key
#'
#' @section Methods:
#' \describe{
#'   \item{\code{PaillierKeyPair$new(modulusBits)}}{Create a new private key with specified
#'         number of modulus bits}
#'   \item{\code{PaillierKeyPair$getPrivateKey()}}{Return the private key}
#' }
#'
#' @examples
#' keys <- PaillierKeyPair$new(1024)
#' keys$pubkey
#' keys$getPrivateKey()
#'
#' @export
#' @format An \code{\link{R6Class}} generator object
PaillierKeyPair <- R6Class("PaillierKeyPair",
                          private = list(
                              privkey = NULL
                          ),
                          public = list(
                              ## fields
                              pubkey = NULL,
                              ## methods
                              initialize = function(modulusBits) {
                                  repeat {
                                      repeat {
                                          p <- random.bigz(nBits = modulusBits %/% 2)
                                          if (isprime(n = p, reps = 10))
                                              break
                                      }
                                      repeat {
                                          q <- random.bigz(nBits = modulusBits %/% 2)
                                          if (isprime(n = q, reps = 10))
                                              break
                                      }
                                      n <- mul.bigz(p, q)
                                      if ( ( p != q ) && ( sizeinbase(a = n, b = 2) == modulusBits) )
                                          break
                                  }

                                  pubkey <- PaillierPublicKey$new(modulusBits, n)
                                  self$pubkey <- pubkey
                                  lambda <- lcm.bigz(sub.bigz(p, ONE), sub.bigz(q, ONE))
                                  private$privkey <- PaillierPrivateKey$new(lambda, pubkey)
                              },
                              getPrivateKey = function() {
                                  private$privkey
                              })
                          )



