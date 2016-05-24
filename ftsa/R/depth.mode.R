`depth.mode` <- function(data, trim = 0.25, h = 0.15, mdist = NULL, scale = TRUE, x = NULL, graph = FALSE, ...){
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
     mtrim = apply(functions[lista,],2,mean)
     if(graph){
        cgray = 1-(ans-min(ans))/(max(ans)-min(ans))
        if(nc == 2){
           plot(range(functions[,1]), range(functions[,2]), type = "n")
           text(functions[,1], functions[,2], round(ans,3), cex = 0.75)
           points(rbind(mtrim), pch = 19, col = gray(2 * trim), cex = 2)
           points(rbind(med), col = 3, pch = 20, cex = 2)
        } 
        else{
           plot(range(x), range(functions), type = "n", xlab = "t",ylab = "X(t)", main = "Mode Depth")
           for(i in 1:nr){
               lines(x, functions[i,], col = gray(cgray[i]))
           }
           lines(x, mtrim,lwd = 2, col = "yellow")
           lines(x,med, col = "red", lwd = 2)                   
           legend("topleft", legend = c(paste("Trim", trim * 100, "\u0025", sep = ""), "Median"), lwd = 2,
        		  col = c("yellow", "red"))  
        }
     }
 return(list("median" = med, "lmed" = k, "mtrim" = mtrim, "ltrim" = lista, "prof" = ans, "dist" = mdist))
}
