##
## XML-RPC method: getIntermediateResultsNonBlocking()
##
NgetIntermediateResultsNonBlocking <- function (obj, offset = NULL, convert = TRUE){
    if(requireNamespace("XMLRPC", quietly = TRUE)) {
        if (!(class(obj) == "NeosJob")) {
            stop("\nObject 'obj' is not of class 'NeosJob'.\n")
        }
        call <- match.call()
        jobnumber <- obj@jobnumber
        password <- obj@password
        nc <- obj@nc
        if(is.null(offset)){
            offset <- as.integer(0)
        } else {
            offset <- as.integer(offset)
        }
        ans <- XMLRPC::xml.rpc(url = nc@url, method = "getIntermediateResultsNonBlocking",
                               .args = list(jobnumber = jobnumber, password = password, offset = offset),
                               .convert = FALSE, .opts = nc@curlopts, .curl = nc@curlhandle)
        tmp <- xmlToList(xmlRoot(xmlTreeParse(ans)))
        offset <- as.integer(tmp[2, ])
        if (convert) {
            tmp1 <- tmp[1, ]
            if(!is.null(tmp1$params)){
                tmp1 <- gsub("\\n", "", tmp1)
                class(tmp1) <- "base64"
                ans <- base64(tmp1)
            } else {
                ans <- "\nNothing left to return from NEOS.\n"
            }
        }
        res <- new("NeosOff", ans = ans, offset = offset, jobnumber = jobnumber, password = password,
                   method = "getIntermediateResultsNonBlocking", call = call, nc = nc)
        return(res)
    } else {
        stop("Package 'XMLRPC' not available, please install first from Omegahat.")
    }
}
