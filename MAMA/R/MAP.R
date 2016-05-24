plotpattern <- function (intx, method) 
{
    if (method == 1) {
        matplot(t(intx)[5:8, ], type = "l", ylab = "pvalue", 
            xlab = "pattern ID", lwd = 3, lty = c(1, 2, 3, 5))
        legend(x = "topleft", colnames(intx)[5:8], lty = c(1, 
            2, 3, 5), col = 1:4, lwd = 2)
        abline(h = 0.05, col = "red")
    }
    if (method == 2) {
      plot(1:4, type = "n", yaxt = "n", xlab = "p-value", ylab = "", xlim = c(0, 1), frame.plot = F)
        for (i in 5:8) {
          lines(x = c(0,1), y = c(i - 4, i - 4), col = "black")
          points(x = intx[, i], y = rep(i - 4 ,nrow(intx)), pch = 19, cex = 2,  col = c(1:nrow(intx)))  
        }
      abline(v=0.05, col = "red") 
      legend(x = "topright", legend = c(as.character(intx[, 1])), 
            col = 1:nrow(intx), pch = 19,  bty = "o", inset = 0.1)    
      axis(2, at = c(1:4), labels = names(intx)[5:8])
   
    }
}
test.group.shuffle <- function (x, varname, minSampleNum = 3, method = "t", B = 1000) 
{
out<-list()
cat("Permuting class labels in dataset: \n")
for (j in 1:length(datanames(x)) )
{
    statX <- NULL
    cat(datanames(x)[[j]], "\n")
    for (i in 1:B) {
        
        ALLtype <- clinical(x)[[j]][, grep(varname, names(clinical(x)[[j]]))]
        stat <- entitybuild2(expr.mat = GEDM(x)[[j]], 
            ALLtype = ALLtype, type = levels(ALLtype), minSampleNum = minSampleNum, 
            method = method, random = TRUE)
        statX <- cbind(statX, stat)
    }
    entity <- entitybuild2(expr.mat = GEDM(x)[[j]], 
        	ALLtype = ALLtype, type = levels(ALLtype), minSampleNum = minSampleNum, 
        	method = method, random = FALSE)
    statX <- statX * (-1)
    entity$stat <- entity$stat * (-1)
    temp <- entity$entity[1]
    entity$entity[1] <- entity$entity[2]
    entity$entity[2] <- temp
    out[[j]]<-list(entity = entity, statX = statX)
}
return(out)
}


 
