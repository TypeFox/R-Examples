lab.mat <-
function(plotlab){
     plotlab <- unique(plotlab)
     nx <- max(as.numeric(substring(plotlab,
           regexpr("X", plotlab)+1,regexpr("Y", plotlab)-1)))
     ny <- max(as.numeric(substring(plotlab,
           regexpr("Y", plotlab)+1,nchar(plotlab))))
     xxx <- expand.grid(1:nx, 1:ny)
     xlab <- paste("X",xxx[,1], sep = "")
     ylab <- paste("Y",xxx[,2], sep = "")
     res00 <- paste(xlab, ylab, sep = "")
     res <- res00
     dim(res) <- c(nx, ny)
     return(res)
}

