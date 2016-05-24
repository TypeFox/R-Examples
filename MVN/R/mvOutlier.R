mvOutlier <-
function (data, qqplot = TRUE, alpha = 0.5, tol = 1e-25, method = c("quan", "adj.quan"), label = TRUE, position = NULL, offset = 0.5)
{
    
        if (!is.data.frame(data) && !is.matrix(data)) stop('Input must be one of classes \"data frame\" or \"matrix\"')
  
        if (dim(data)[2] < 2 || is.null(dim(data))) {stop("number of variables must be equal or greater than 2")}
        
  
        dataframe=as.data.frame(data)
        dname <- deparse(substitute(data))
        method <- match.arg(method)
        
        n <- dim(data)[1]
        p <- dim(data)[2]
        
    
        
        covr <- covMcd(data, alpha = alpha)
        mah <- mahalanobis(data, center = covr$center, cov = covr$cov, tol = tol)
        d <- mah
        sortMah <- data.frame(sort(mah, decreasing = TRUE)) # sorted Mahalanobis' distances in increasing order
        
        out <-  cbind(rownames(sortMah),round(sortMah,3), NA)
        colnames(out) <- c("Observation","Mahalanobis Distance", "Outlier")

        
        if (method=="adj.quan"){
            crt <- arw(x=data, m0=covr$center, c0=covr$cov, alpha = 0.025)$cn
            
            for(i in 1:n){
                {
                if (sortMah[i,] > crt){
                    out[i,3] <- "TRUE"
                    
                } else
                {
                    out[i,3] <- "FALSE"
                }}
            }
            
            
            
            if (qqplot){
                d <- mah
                r <- rank(d)
                chi2q <- qchisq((r-0.5)/n,p)
                
                colors = NULL
                for (i in 1:n) {
                    if (d[i] > crt) colors[i] = "red" else colors[i] = "black"
                }
                
                plot(d, chi2q , pch = 16, main = "Adjusted Chi-Square Q-Q Plot",
                xlab = "Robust Squared Mahalanobis Distance",ylab="Chi-Square Quantile", col=colors)
                abline(v=crt, lwd = 2, col = "blue")
                tbl = table(out[,3])
                
                legend("topleft",legend=c(paste("Outliers (n=",if(is.na(tbl[2])) 0 else tbl[2],")",sep=""),paste("Non-outliers (n=",if(is.na(tbl[1])) 0 else tbl[1],")",sep="")),
                col=c("red","black"), pch=16, bty="n",)
                
                if(label && is.element("TRUE", out[, 3])) {
                    labelOutlier <- rownames(out)[out[,3] == TRUE]
                    xCoord <- out[out[,3] == TRUE,2]
                    yCoord <- sort(chi2q,decreasing = T)[1:length(xCoord)]
                    text(xCoord, yCoord, labelOutlier, pos = position, offset = offset)
                 }

		if (max(d) >= crt) {text(crt-0.2,2,paste("Quantile: ", round(crt,3)), srt=90, pos=3, col="blue")}
            }
                      
            newData <- out[out$Outlier %in% "FALSE",]
            ind <- sort(row.names(newData))
            newData <- data[ind,]
            
            result <- list(out, newData)
            names(result) <- c("outlier", "newData")
            
        }
        
        if (method=="quan"){
            
            chiSq <- qchisq(0.975, p)

            for(i in 1:n){
                {
                    if (sortMah[i,] > chiSq){
                        out[i,3] <- "TRUE"
                        
                    } else
                    {
                        out[i,3] <- "FALSE"
                    }}
            }
            
            if (qqplot){
                d <- mah
                r <- rank(d)
                chi2q <- qchisq((r-0.5)/n,p)
                
                
                colors = NULL
                for (i in 1:n) {
                    if (d[i] > chiSq) colors[i] = "red" else colors[i] = "black"
                }
                
                plot(d, chi2q , pch = 16, col=colors, main = "Chi-Square Q-Q Plot",
                xlab = "Robust Squared Mahalanobis Distance",ylab="Chi-Square Quantile")
                abline(v=chiSq, lwd = 2, col = "red")
                
                tbl = table(out[,3])
                
                legend("topleft",legend=c(paste("Outliers (n=",if(is.na(tbl[2])) 0 else tbl[2],")",sep=""),paste("Non-outliers (n=",if(is.na(tbl[1])) 0 else tbl[1],")",sep="")),
                col=c("red","black"), pch=16, bty="n",)
                
                if(label && is.element("TRUE", out[, 3])) {
                     labelOutlier <- rownames(out)[out[,3] == TRUE]
                     xCoord <- out[out[,3] == TRUE,2]
                     yCoord <- sort(chi2q,decreasing = T)[1:length(xCoord)]
                     text(xCoord, yCoord, labelOutlier, pos = position, offset = offset)
                }
                
		if (max(d) >= chiSq) {text(chiSq-0.2,2,paste("Quantile: ", round(chiSq,3)),srt = 90, pos = 3, col="red")}
            }

            
            newData <- out[out$Outlier %in% "FALSE",]
            ind <- sort(row.names(newData))
            newData <- data[ind,]
            
            result <- list(out, newData)
            names(result) <- c("outlier", "newData")
            
        }
        
    return(result)

}

