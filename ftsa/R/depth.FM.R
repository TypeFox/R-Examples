`depth.FM` <- function(data, trim = 0.25, xeps = 0.00000001, x = NULL, graph = FALSE){
    functions = t(data$y)
    nrow <- dim(functions)[1]
    ncol <- dim(functions)[2]
    if(is.null(nrow) && is.null(ncol))
       stop("I do not have a matrix")
    if(is.null(x)) 
       x = 1:ncol
    d <- matrix(NA, nrow = nrow, ncol = ncol)
    for(i in 1:ncol){
        Fn = ecdf(functions[,i])
        d[,i] = 1 - abs(0.5 - Fn(functions[,i]))
    }
    ans <- apply(d, 1, mean)
    k = which.max(ans)
    med = functions[k,]
    lista = which(ans >= stats::quantile(ans, probs = trim, na.rm = TRUE))
    mtrim = apply(functions[lista,], 2, mean)
    if(graph){
       cgray = 1 - (ans - min(ans)) / (max(ans) - min(ans))
       if(ncol == 2){
          plot(range(functions[,1]), range(functions[,2]), type = "n")
          text(functions[,1], functions[,2], round(ans, 3), cex = 0.75)
          points(rbind(mtrim), pch = 19, col = gray(2 * trim), cex = 2)
          points(rbind(med), col = 3, pch = 20, cex = 2)
       } 
       else{
          plot(range(x), range(functions), type = "n", xlab = "t", ylab = "X(t)", main = "FM Depth")
          for(i in 1:nrow){
              lines(x, functions[i,], col = gray(cgray[i]))
          }
          lines(x, mtrim, lwd = 2, col = "yellow")
          lines(x, med, col = "red", lwd = 2) 
       }    
	       legend("topleft", legend = c(paste("Trim", trim * 100, "\u0025", sep = ""), "Median"), lwd = 2,
             	  col = c("yellow", "red"))  
  }
 return(list("median" = med, "lmed" = k, "mtrim" = mtrim, "ltrim" = lista, "prof" = ans))
}
