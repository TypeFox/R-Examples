`mean.fts` <- function (x, method = c("coordinate", "FM", "mode", "RP", "RPD", "radius"), 
                         na.rm = TRUE, alpha, beta, weight, ...) 
{
   if (class(x)[1] == "fts"|class(x)[1] == "fds"|class(x)[1] == "sfts"){
       method = match.arg(method)
       if (method == "coordinate"){
          loc <- rowMeans(x$y, na.rm = na.rm)
       }
       if (method == "FM"){
          loc <- depth.FM(x)$mtrim
       }
       if (method == "mode"){
          loc <- depth.mode(x)$mtrim
       }
       if (method == "RP"){
          loc <- depth.RP(x)$mtrim
       }
       if (method == "RPD"){
          loc <- depth.RPD(x)$mtrim    
       }
       if (method == "radius"){
       	  loc <- depth.radius(x, alpha, beta, weight)$mtrim
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
