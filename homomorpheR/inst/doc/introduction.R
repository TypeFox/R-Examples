## ----echo=F--------------------------------------------------------------
### get knitr just the way we like it

knitr::opts_chunk$set(
  message = FALSE,
  warning = FALSE,
  error = FALSE,
  tidy = FALSE,
  cache = FALSE
)


## ------------------------------------------------------------------------
keyPair <- PaillierKeyPair$new(modulusBits = 1024)

## ------------------------------------------------------------------------
keyPair

## ------------------------------------------------------------------------
encryptAndDecrypt <- function(x) keyPair$getPrivateKey()$decrypt(keyPair$pubkey$encrypt(x))

## ------------------------------------------------------------------------
a <- gmp::as.bigz(1273849)
identical(a + 10, encryptAndDecrypt(a + 10))

## ------------------------------------------------------------------------
m <- lapply(1:100, function(x) random.bigz(nBits = 512))
md <- lapply(m, encryptAndDecrypt)
identical(m, md)

