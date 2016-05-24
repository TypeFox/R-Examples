`median.fts` <- function (x, method = c("hossjercroux", "coordinate", "FM", "mode", "RP", "RPD", "radius"), alpha, beta, weight, ...) 
{
    if (class(x)[1] == "fts"|class(x)[1] == "fds"|class(x)[1] == "sfts"){
        method = match.arg(method)
        if (method == "hossjercroux"){
            loc <- L1median2(t(x$y), method = "hossjercroux")
        }
        if (method == "coordinate"){
            loc <- apply(x$y,1,median.default)
        }
        if (method == "FM"){
            loc <- depth.FM(x)$median
        }
        if (method == "mode"){
            loc <- depth.mode(x)$median
        }
        if (method == "RP"){
            loc <- depth.RP(x)$median
        }
        if (method == "RPD"){    
            loc <- depth.RPD(x)$median
        }
        if (method == "radius"){
        	loc <- depth.radius(x, alpha, beta, weight)$median
        }
        if (class(x)[1] == "fds"){
            warning("Object is not a functional time series.")
        }        
        return(list(x = x$x, y = loc))
    }
    else {
         stop("Not a functional object.")
    }
}
