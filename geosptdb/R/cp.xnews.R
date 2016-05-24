assign("cp.xnews",
       function(newdata,eigenvalues,data,trend, ...){                                                
         data[length(eigenvalues)+1, ] <- newdata[1,]
         d <- gower.dist(data[1:length(eigenvalues),],data[length(eigenvalues)+1,],KR.corr=F, ...)
         b <- diag(tcrossprod(trend))                                                                        
         x.new <- (1/2) * diag(eigenvalues[1:ncol(trend)]^(-1)) %*% t(trend) %*% (b - d)
         x.new
       }
)