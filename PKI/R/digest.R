PKI.digest <- function(what, hash=c("SHA1", "SHA256", "MD5")) {
  hash <- pmatch(hash, c("SHA1", "SHA256", "MD5"))[1]
  if (is.na(hash)) stop("invalid hash specification")
  .Call(PKI_digest, what, hash)
}
