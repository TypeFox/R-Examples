`var.fts` <- function(x, method = c("coordinate", "FM", "mode", "RP", "RPD", "radius"), trim = 0.25, alpha, weight, ...)
{
   if (class(x)[1] == "fts"|class(x)[1] == "fds"|class(x)[1] == "sfts"){
       functions = t(x$y) 
       n = dim(functions)[1]
       p = dim(functions)[2]
       method = match.arg(method)
       if (method == "coordinate"){
           loc = (n - 1) * apply(functions, 2, var) / n
       }
       if (method == "FM"){
           lista = depth.FM(x, trim = trim)$ltrim
           loc = func.var(functions[lista,])
       }
       if (method == "mode"){
           lista = depth.mode(x, trim = trim)$ltrim
           loc = func.var(functions[lista,])
       }
       if (method == "RP"){
           lista = depth.RP(x, trim = trim)$ltrim
           loc = func.var(functions[lista,])
       }
       if (method == "RPD"){
           lista = depth.RPD(x, trim = trim)$ltrim
           loc = func.var(functions[lista,])
       }   
       if (method == "radius"){
       	   lista = which(depth.radius(x, alpha, trim, weight)$weight==1)
       	   loc = func.var(functions[lista,])
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
