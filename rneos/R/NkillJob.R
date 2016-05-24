##
## XML-RPC method: killJob()
##
NkillJob <- function(obj, killmsg = NULL, convert = TRUE){
    if(requireNamespace("XMLRPC", quietly = TRUE)) {
        if(!(class(obj) == "NeosJob")){
            stop("\nObject 'obj' is not of class 'NeosJob'.\n")
        }
        call <- match.call()
        jobnumber <- obj@jobnumber
        password <- obj@password
        nc <- obj@nc
        if(is.null(killmsg)){
            killmsg <- ""
        } else {
            killmsg <- as.character(killmsg)
        }
        ans <- XMLRPC::xml.rpc(url = nc@url, method = "killJob",
                               .args = list(jobnumber = jobnumber, password = password, killmsg = killmsg),
                               .convert = convert, .opts = nc@curlopts, .curl = nc@curlhandle)
        res <- new("NeosAns", ans = ans, method = "killJob", call = call, nc = nc)
        return(res)
    } else {
        stop("Package 'XMLRPC' not available, please install first from Omegahat.")
    }
}
