##
## XML-RPC method: welcome()
##
Nwelcome <- function(convert = TRUE, nc = CreateNeosComm()){
    if(requireNamespace("XMLRPC", quietly = TRUE)) {
        if(!(class(nc) == "NeosComm")){
            stop("\nObject provided for 'nc' must be of class 'NeosComm'.\n")
        }
        call <- match.call()
        ans <- XMLRPC::xml.rpc(url = nc@url, method = "welcome", .convert = convert,
                               .opts = nc@curlopts, .curl = nc@curlhandle)
        res <- new("NeosAns", ans = ans, method = "welcome", call = call, nc = nc)
        return(res)
    } else {
        stop("Package 'XMLRPC' not available, please install first from Omegahat.")
    }
}
