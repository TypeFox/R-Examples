signRequest  <- function(url, params, consumerKey, consumerSecret,
                         oauthKey = "", oauthSecret = "", httpMethod = "GET",
                         signMethod = "HMAC", nonce = genNonce(),
                         timestamp = Sys.time(),
                         escapeFun = encodeURI,
                         handshakeComplete=TRUE) {
  ## Sign an request made up of the URL, the parameters as a named character
  ## vector the consumer key and secret and the token and token secret.
  httpMethod <- toupper(httpMethod)
  signMethod <- toupper(signMethod)

  params["oauth_nonce"] <- nonce
  params["oauth_timestamp"] <- as.integer(timestamp)

  token <- oauthKey
  if(!is.null(token) && !is.na(token) && token != "")
     params["oauth_token"] <- token

  params["oauth_consumer_key"] <- consumerKey
  params["oauth_signature_method"] <- switch(signMethod,
                                             HMAC = 'HMAC-SHA1',
                                             RSA  = 'RSA-SHA1',
                                             text = 'PLAINTEXT',
                                             stop("Unsupported signature method: ", signMethod)
                                             )
  params["oauth_version"] <- '1.0'

  args <- escapeFun(normalizeParams(params, escapeFun), post.amp = TRUE)

  if(is.null(oauthSecret))
     oauthSecret <- ""

  okey <- paste(sapply(c(consumerSecret, oauthSecret), escapeFun),
                collapse = "&")
  ## note that we don't escape the args string again.
  odat <- paste(c(sapply(c(httpMethod, url), escapeFun), args),
                collapse = "&")

  sig <- signString(odat, okey, signMethod)

  params["oauth_signature"] <- sig
  ##
  return(params)
}

signString <- function(str, key, method) {
  ## Perform the actual computation to get the signature of the data
  sigFunc <- switch(toupper(method),
                    HMAC = signWithHMAC,
                    RSA = signWithRSA,
                    text = signWithPlaintext,
                    stop("No signature method for ", method)
                    )
  sigFunc(key, str)
}

genNonce <- function(len = 15L + sample(1:16, 1L)) {
  ## Get a random sequence of characters.
  ## Nonce - number used only once.
  els <- c(letters, LETTERS, 0:9, "_")
  paste(sample(els, len, replace = TRUE), collapse = "")
}

signWithHMAC <- function(key, data) {
  ## From Ozaki Toru's code at https://gist.github.com/586468
  blockSize <- 64
  hashlength <- 20
  innerpad   <- rawToBits(as.raw(rep(0x36, blockSize)))
  outerpad   <- rawToBits(as.raw(rep(0x5C, blockSize)))
  zero       <- rep(0 ,64)

  HexdigestToDigest <- function(digest) {
    as.raw(strtoi(substring(digest, (1:hashlength)*2-1,
                            (1:hashlength)*2), base=16))
  }

  mac <- function(pad, text) {
    HexdigestToDigest(digest(append(packBits(xor(key, pad)), text),
                             algo='sha1', serialize=FALSE))
  }
  
  if(nchar(key) >= 64) {
    keyDigested <- digest(key, algo="sha1", serialize=FALSE)
    key <- intToUtf8(strtoi(HexdigestToDigest(keyDigested), base=16))
  }
  key <- rawToBits(as.raw(append(utf8ToInt(key), zero)[1:blockSize]))

  base64(mac(outerpad, mac(innerpad, charToRaw(data))))[1]
}

signWithRSA <- function(key, data) {
  stop("RSA signature not implemented")
}

signWithPlaintext <- function(key, data) {
  key
}

## this function is derived from utils::URLencode
## Characters not in the unreserved character set ([RFC3986] section 2.3) MUST be encoded
##   unreserved = ALPHA, DIGIT, '-', '.', '_', '~'
## cf. http://oauth.net/core/1.0/#encoding_parameters
encodeURI <- function(URI, ...) {
  if (!is.character(URI)) {
    URI
  } else {
    OK <- "[^-A-Za-z0-9_.~]"
    x <- strsplit(URI, "")[[1L]]
    z <- grep(OK, x)
    if (length(z)) {
      y <- sapply(x[z], function(x) paste("%", toupper(as.character(charToRaw(x))),
                                          sep = "", collapse = ""))
      x[z] <- y
    }
      paste(x, collapse = "")
  }
}

## cf. http://tools.ietf.org/html/rfc5849#section-3.4.1.3.2
normalizeParams <- function(params, escapeFun) {
  ## we escape the values of the parameters in a special way that escapes
  ## the resulting % prefix in the escaped characters, e.g. %20 becomes
  ## %2520 as %25 is the escape for %
  names(params) <- sapply(names(params), escapeFun, post.amp = TRUE)
  params <- sapply(params, escapeFun, post.amp = TRUE)
  ## If two or more parameters share the same name, they are sorted by their value.
  params <- params[order(names(params), params)]
  return(paste(names(params), params, sep = "=", collapse = "&"))
}
