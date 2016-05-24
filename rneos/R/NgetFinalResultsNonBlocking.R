##
## XML-RPC method: getFinalResultsNonBlocking()
##
NgetFinalResultsNonBlocking <- function(obj, convert = TRUE){
    if(requireNamespace("XMLRPC", quietly = TRUE)) {
        if (!(class(obj) == "NeosJob")) {
            stop("\nObject 'obj' is not of class 'NeosJob'.\n")
        }
        call <- match.call()
        jobnumber <- obj@jobnumber
        password <- obj@password
        nc <- obj@nc
        ans <- XMLRPC::xml.rpc(url = nc@url, method = "getFinalResultsNonBlocking",
                               .args = list(jobnumber = jobnumber, password = password),
                               .convert = FALSE, .opts = nc@curlopts, .curl = nc@curlhandle)
        if(convert){
            tmp <- xmlValue(xmlRoot(xmlParse(ans)))
            tmp <- gsub("\\n", "", tmp)
            class(tmp) <- "base64"
            ans <- base64(tmp)
        }
        res <- new("NeosAns", ans = ans, method = "getFinalResultsNonBlocking", call = call, nc = nc)
        return(res)
    } else {
        stop("Package 'XMLRPC' not available, please install first from Omegahat.")
    }
}
