PKI.encrypt <- function(what, key, cipher = NULL) .Call(PKI_encrypt, what, key, cipher)

PKI.decrypt <- function(what, key, cipher = NULL) .Call(PKI_decrypt, what, key, cipher)
