#' SFTP location
#'
#' This class represents a directory on a remote server that is accessed
#' via the SFTP protocol. The only SSH authentication currently supported
#' the public RSA key method.
#'
#' The \code{show} and \code{as.character} methods for this class have been
#' adapted to reveal the URL of the directory it represents.
#'
#' The \code{meta} method returns a \code{data.frame} with filenames and
#' an \code{is.dir} attribute. It currently returns wrong values when the
#' remote directory contains files with spaces in their file name.
#' 
#' @examples
#' getSlots("SftpLocation")
#'
#' @seealso \code{\link{sftpdir}}
#'
#' @name SftpLocation-class
#' @rdname SftpLocation-class
#' @exportClass SftpLocation
setClass(
  Class="SftpLocation", 
  representation=representation(pubkey="character", privatekey="character", uri="character", keypasswd="character"), 
  contains="Location",
  validity=function(object) {
    if(!file.exists(object@pubkey)) stop("could not find public key file '", object@pubkey, "'.")
    if(!file.exists(object@privatekey)) stop("could not find private key file '", object@privatekey, "'.")
    if(strhead(object@uri, 7)!="sftp://") object@uri <- paste(object@uri, "sftp://", object@uri, sep="")
    if(strtail(object@uri, 1)!="/") object@uri <- paste(object@uri, "/", sep="")
  }
)

#' The sftpdir function creates a SftpLocation object.
#' 
#' @param uri          character, the address of the remote host in the form
#'                     sftp://username@@the.host.com:port/remote/dir/
#' @param pubkey       character, pointing to the public key file. Required.
#' @param privatekey   character, pointing to the private key file. Required.
#' @param keypasswd    character, optional key password
#' @param clss         character, optional class name. Default is "SftpLocation".
#'
#' @rdname SftpLocation-class
#' @export
sftpdir <- function(uri, pubkey, privatekey, keypasswd="", clss="SftpLocation") 
   new(clss, uri=uri, pubkey=pubkey, privatekey=privatekey, keypasswd=keypasswd)

#' @rdname show-methods
#' @name show
#' @export
#' @docType methods
#' @aliases show show,SftpLocation-method
setMethod(
  f="show",
  signature="SftpLocation",
  definition=function(object) cat(sprintf("<Sftp Directory %s>\n", object@uri))
)

#' @rdname as.character-methods
#' @name as.character
#' @export
#' @docType methods
#' @aliases as.character,SftpLocation-method
setMethod(
  f="as.character",
  signature="SftpLocation",
  definition=function(x) x@uri
)

#' @rdname meta-methods
#' @name meta
#' @export
#' @docType methods
#' @aliases meta meta,SftpLocation-method
setMethod(
    f="meta",
    signature="SftpLocation",
    definition=function(self) {
        h <- basicTextGatherer()
        opts <- list(
            ssh.public.keyfile = self@pubkey, 
            ssh.private.keyfile = self@privatekey, 
            writefunction = h$update
        ) 
        if(self@keypasswd!="") opts$keypasswd <- self@keypasswd
        RCurl::curlPerform(
            url=self@uri, 
            .opts = opts, curl = RCurl::getCurlHandle()
        )
        ret <- h$value()
        ret <- strsplit(ret, "\n")[[1]]
        ret <- strsplit(ret, " ")
        ret <- lapply(ret, function(s) data.frame(isdir=s[[1]], name=s[[length(s)]]))
        ret <- Reduce(rbind, ret, init=data.frame())
        ret <- transform(ret, isdir=ifelse(substring(isdir,1,1)=="d", TRUE, FALSE))
        ret <- subset(ret, name!="." & name!="..")
        rownames(ret) <- NULL
        return(ret)
    }
)

#' @rdname put-methods
#' @name put
#' @export
#' @docType methods
#' @aliases put put,FileTarget,SftpLocation-method
setMethod(
    f="put",
    signature=c(target="FileTarget", where="SftpLocation"),
    definition=function(target, where, ...) {
        opts <- list(ssh.public.keyfile = where@pubkey, ssh.private.keyfile = where@privatekey, ftp.create.missing.dirs=TRUE)
        if(where@keypasswd!="") opts$keypasswd <- where@keypasswd
        RCurl::ftpUpload(
            what = target@filename,
            to = paste(where@uri, target@name, sep=""),
            .opts = opts
        )
    }
)

#' @rdname put-methods
#' @name put
#' @export
#' @docType methods
#' @aliases put put,character,SftpLocation-method
setMethod(
    f="put",
    signature=c(target="character", where="SftpLocation"),
    definition=function(target, where, ...) put(filetarget(name=basename(target), filename=normalizePath(target)), where, ...)
)

#' @rdname query-methods
#' @name query
#' @export
#' @docType methods
#' @aliases query query,SftpLocation,character-method
setMethod(
    f="query",
    signature=c(self="SftpLocation", resource="character"),
    definition=function(self, resource, verbose=getOption("verbose"), ...) {
        uri <- file.path(self@uri, resource)
        opts <- list(
            ssh.public.keyfile = self@pubkey, 
            ssh.private.keyfile = self@privatekey
        ) 
        if(self@keypasswd!="") opts$keypasswd <- self@keypasswd
        RCurl::getURL(
            url=uri, 
            .opts = opts, curl = RCurl::getCurlHandle()
        )
    }
)
