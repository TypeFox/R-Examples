#' Bcrypt password hashing
#'
#' Bcrypt is used for secure password hashing. The main difference with
#' regular digest algorithms such as MD5 or SHA256 is that the bcrypt
#' algorithm is specifically designed to be CPU intensive in order to
#' protect against brute force attacks. The exact complexity of the
#' algorithm is configurable via the \code{log_rounds} parameter. The
#' interface is fully compatible with the Python one.
#'
#' The \code{hashpw} function calculates a hash from a password using
#' a random salt. Validating the hash is done by reshashing the password
#' using the hash as a salt. The \code{checkpw} function is a simple
#' wrapper that does exactly this.
#'
#' \code{gensalt} generates a random text salt for use with \code{hashpw}.
#' The first few characters in the salt string hold the bcrypt version number
#' and value for \code{log_rounds}. The remainder stores 16 bytes of base64
#' encoded randomness for seeding the hashing algorithm.
#'
#' @importFrom openssl rand_bytes
#' @export
#' @param log_rounds integer between 4 and 31 that defines the complexity of
#' the hashing, increasing the cost as \code{2^log_rounds}.
#' @rdname bcrypt
#' @name bcrypt
#' @examples # Secret message as a string
#' passwd <- "supersecret"
#'
#' # Create the hash
#' hash <- hashpw(passwd)
#' hash
#'
#' # To validate the hash
#' identical(hash, hashpw(passwd, hash))
#'
#' # Or use the wrapper
#' checkpw(passwd, hash)
#'
#' # Use varying complexity:
#' hash11 <- hashpw(passwd, gensalt(11))
#' hash12 <- hashpw(passwd, gensalt(12))
#' hash13 <- hashpw(passwd, gensalt(13))
#'
#' # Takes longer to verify (or crack)
#' system.time(checkpw(passwd, hash11))
#' system.time(checkpw(passwd, hash12))
#' system.time(checkpw(passwd, hash13))
gensalt <- function(log_rounds = 12){
  stopifnot(is.numeric(log_rounds))
  csalt <- openssl::rand_bytes(16)
  encode_salt(csalt, as.integer(log_rounds))
}

#' @useDynLib bcrypt R_encode_salt
encode_salt <- function(csalt, log_rounds){
  .Call(R_encode_salt, csalt, log_rounds)
}

#' @useDynLib bcrypt R_hashpw
#' @export
#' @param salt a salt generated with \code{gensalt}.
#' @rdname bcrypt
hashpw <- function(password, salt = gensalt()){
  password <- enc2utf8(password)
  salt <- enc2utf8(salt)
  .Call(R_hashpw, password, salt)
}

#' @export
#' @param password the message (password) to encrypt
#' @param hash the previously generated bcrypt hash to verify
#' @rdname bcrypt
checkpw <- function(password, hash){
  hash <- enc2utf8(hash)
  password <- enc2utf8(password)
  identical(hash, hashpw(password, hash))
}
