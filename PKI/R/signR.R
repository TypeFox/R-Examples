# tar one regular file from raw payload
tar1 <- function(name, what, mode=0x180) {
    size <- length(what)
    header <- raw(512L)
    fn <- charToRaw(name)
    header[seq_along(fn)] <- fn
    header[101:107] <- charToRaw(sprintf("%07o", mode))
    header[137:147] <- charToRaw(sprintf("%011o", as.integer(Sys.time())))
    header[157L] <- charToRaw("0") # regular file
    header[125:135] <- charToRaw(sprintf("%011o", as.integer(size)))
    header[149:156] <- charToRaw(" ")
    checksum <- sum(as.integer(header))%%2^24
    header[149:154] <- charToRaw(sprintf("%06o", as.integer(checksum)))
    header[155L] <- as.raw(0L)
    bsize <- ceiling(size / 512L) * 512L
    padding <- raw(bsize - size)
    c(header, what, padding)
}

chunk <- 4194304L ## 4Mb .. as good as any value ...

PKI.sign.tar <- function(tarfile, key, certificate, output=tarfile) {
    io <- file
    file <- file(tarfile, "rb")
    on.exit(if (!is.null(file)) close(file))
    magic <- readBin(file, raw(), n = 3)
    if (all(magic[1:2] == c(31, 139)) || all(magic[1:2] == c(31, 157)))
        io <- gzfile
    else if (rawToChar(magic[1:3]) == "BZh") 
        io <- bzfile
    else if (rawToChar(magic[1:5]) == "\xfd7zXZ") 
        io <- xzfile
    close(file)
    file <- NULL
    file <- io(tarfile, "rb")
    payload <- raw(0)
    while (length(r <- readBin(file, raw(), chunk))) payload <- c(payload, r)
    close(file)
#   FIXME: if we want the .signature to be visible, we need to strip padding to inject new file !
    file <- NULL
    sign <- PKI.sign(payload, key, "SHA1")
    ## SEQ(BIT STREAM sig, subjectPubKeyInfo[if not cert], cert[optional, if present])
    a <- if (missing(certificate)) 
        ASN1.encode(list(ASN1.item(sign, 3L), ASN1.decode(PKI.save.key(key, "DER", FALSE))))
    else
        ASN1.encode(list(ASN1.item(sign, 3L), ASN1.item(raw(0), 0L), ASN1.decode(attr(certificate, "crt.DER"))))
    payload <- c(payload, tar1(".signature", a), as.raw(rep(0L, 1024L)))
    if (inherits(output, "connection")) {
        writeBin(payload, output)
        return(invisible(output))
    }
    if (is.raw(output)) return(payload)
    file <- io(as.character(output), "wb")
    writeBin(payload, file)
    close(file)
    file <- NULL
    invisible(output)
}

PKI.verify.tar <- function(tarfile, key, silent = FALSE, enforce.cert = FALSE) {
    if (is.raw(tarfile))
        payload <- tarfile
    else {
        io <- file
        file <- file(tarfile, "rb")
        on.exit(if (!is.null(file)) close(file))
        magic <- readBin(file, raw(), n = 3)
        if (all(magic[1:2] == c(31, 139)) || all(magic[1:2] == c(31, 157)))
          io <- gzfile
        else if (rawToChar(magic[1:3]) == "BZh") 
          io <- bzfile
        else if (rawToChar(magic[1:5]) == "\xfd7zXZ") 
          io <- xzfile
        close(file)
        file <- NULL
        file <- io(tarfile, "rb")
        payload <- raw(0)
        while (length(r <- readBin(file, raw(), chunk))) payload <- c(payload, r)
        if (length(payload) < 1024L) stop("invalid tar format")
        close(file)
        file <- NULL              
    }
    fn <- c(charToRaw(".signature"), raw(1))
    i <- length(payload) - 511L
    n <- length(fn) - 1L
    while (i > 0L) {
        if (identical(payload[seq.int(i, i + n)], fn)) break
        i <- i - 512L
    }
    if (i < 1L) {
        if (!silent) warning("no signature found")
        return(FALSE)
    }
    asn <- try(ASN1.decode(payload[seq.int(i + 512L, length(payload))]), silent=TRUE)
    if (!is.list(asn) || length(asn) < 2L || ASN1.type(asn[[1]]) != 3L) {
        if (!silent) warning("bad signature format")
        return(FALSE)
    }
    sig.key <- NULL
    sig.cert <- NULL
    if (is.list(asn[[2]]) && length(asn[[2]]) == 2L) { ## no certificate, jsut a key
        der <- ASN1.encode(asn[[2]]) ## encode back to DER
        sig.key <- try(PKI.load.key(der, "DER", FALSE), silent=TRUE)
        if (!inherits(sig.key, "public.key"))
            sig.key <- NULL
    } else if (length(asn) > 2L && identical(ASN1.type(asn[[3]]), 3L)) { ## certificate
        sig.cert <- try(PKI.load.cert(asn[[3]], "DER"), silent=TRUE)
        if (!inherits(sig.cert, "X509cert"))
            sig.cert <- NULL
        else
            sig.key <- PKI.pubkey(sig.cert)
    }
    if (missing(key)) {
        if (is.null(sig.key)) {
            if (!silent) warning("no key supplied and no valid key or certificate in the signature")
            return(FALSE)
        } else
            key <- sig.key
    }
    if (!identical(enforce.cert, FALSE) && is.null(sig.cert)) {
        if (!silent) warning("certificate required but no certificate found")
        return(FALSE)
    }
    res <- PKI.verify(payload[seq.int(1L, i - 1L)], asn[[1]], key)
    if (enforce.cert) {
        if (inherits(enforce.cert, "X509cert")) {
            key1 <- PKI.save.key(PKI.pubkey(enforce.cert), "DER", FALSE)
            key2 <- PKI.save.key(key, "DER", FALSE)
            if (!identical(key1, key2)) {
                if (!silent) warning("signed by a different certificate")
                return (FALSE)
            }
        }
        if (!isTRUE(res)) FALSE else sig.cert
    } else res
}
