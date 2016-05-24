plot.gnsc <-
function(x, which.lambda=NULL, ...){
    if(x$nlambda == 1){  
       cat("There is only one thresholding value provided.")
       gnsc.heatmap(x$Thresh.mat[,,1], Colv=NA, Rowv=NA,col = greenred(256), cexRow = .7, cexCol = .7, 
       trace = "none",key=FALSE, main="Classes Significance Level", xlab="sample class",ylab="variable class",
       margins=c(5,7),keysize=1)       
    }
    if(x$nlambda > 1){
       par(mfrow=c(1,2))
       plot(c(1:x$nlambda),c(x$errors),type="b",xlab="Threshold Values", ylab="Predicted Error", main="Threshold Level vs. Predicted Error")       
       if(is.null(which.lambda)) which.lambda = floor(x$nlambda/2)
       if(which.lambda>x$nlambda) which.lambda = floor(x$nlambda/2)
       gnsc.heatmap(x$Thresh.mat[,,which.lambda], Colv=NA, Rowv=NA,col = greenred(256), cexRow = .7, cexCol = .7, 
       trace = "none",key=FALSE, main="Classes Significance Level", xlab="sample class",ylab="variable class",
       margins=c(5,7),keysize=1, dendrogram="none")              
    }
}
