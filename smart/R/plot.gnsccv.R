plot.gnsccv <-
function(x, ...){
    if(x$nlambda == 1){  
       cat("There is only one thresholding value provided.")
       gnsc.heatmap(x$Thresh.mat[,,1], Colv=NA, Rowv=NA,col = greenred(256), cexRow = .7, cexCol = .7, 
       trace = "none",key=FALSE, main="Classes Significance Level", xlab="sample class",ylab="variable class",
       margins=c(5,7),keysize=1, dendrogram="none")       
    }
    if(x$nlambda > 1){
       par(mfrow=c(1,2))
       plot(c(1:x$nlambda),c(x$errors),type="b",xlab="Threshold Values", ylab="Predicted Error", main="Threshold Level vs. Predicted Error")       
       gnsc.heatmap(x$Thresh.mat[,,x$lambda.min], Colv=NA, Rowv=NA,col = greenred(256), cexRow = .7, cexCol = .7, 
       trace = "none",key=FALSE, main="Classes Significance Level", xlab="sample class",ylab="variable class",
       margins=c(5,7),keysize=1, dendrogram="none")              
    }
}
