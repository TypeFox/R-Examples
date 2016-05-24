`depth.mode` <- function(data, trim = 0.25, h = 0.15, mdist = NULL, scale = TRUE, x = NULL, ...){
  functions = t(data$y)
  nr <- nrow(functions)
  nc <- ncol(functions)
  if(is.null(x)) 
     x = 1:nc
     if(is.null(nr) && is.null(nc)) 
        stop("A matrix in first parameter is expected")
     if(is.null(mdist)||dim(mdist)!=c(nr,nr)){
        mdist <- matrix(0.0, nrow = nr, ncol = nr)
        for(i in 1:(nr-1)){
            for(j in (i+1):nr){
                mdist[i,j] <- metri.p(functions[i,], functions[j,],...)
                mdist[j,i] <- mdist[i,j]
            }   
         }
     }
     h = quantile(mdist + diag(Inf,nrow(mdist)), probs = h, na.rm = TRUE)
     ans <- apply(mdist/h, 1, skernel.norm)
     if(scale)
        ans = as.vector(scale(ans, center = min(ans), scale = max(ans) - min(ans)))
     k = which.max(ans)
     med = functions[k,]
     lista = which(ans >= quantile(ans, probs = trim, na.rm = TRUE))
     mtrim = apply(functions[lista,], 2, mean)
     return(list("median" = med, "lmed" = k, "mtrim" = mtrim, 
            "ltrim" = lista, "prof" = ans))
}

