predict.RWBP <-
function(object,top_k=3, type="raw",...){
topOutliers_k <- top_k
OutScores <- object$OutScore

# returns classification
if (type=="raw"){ 
OutScores <- cbind(OutScores ,c(rep(0,object$n))) #initialize 0 scores for all
OutScores[OutScores[,1] %in% head(OutScores,n = topOutliers_k )[,1],3]<-1 #change top records score to 1
OutScores <- OutScores[order(OutScores[,1],decreasing=F), ] #sort list by original index
results<- data.frame(object$data,OutScores[,3]) #Mark outliers in the original data
colnames(results)[ncol(results)] <- "class" # added "class" column: 1 for outlier, 0 for normal
}
# returns probabilities
if (type=="prob"){ 
probs <- cbind(OutScores ,1-(OutScores[,2]-min(OutScores[,2]))/(max(OutScores[,2])-min(OutScores[,2]))) #calculate normalized score value (probability)
probs <- probs [order(probs[,1],decreasing=F), ] #sort list by original index
results <- data.frame(object$data,probs[,3]) #add probability to the original data
colnames(results)[ncol(results)] <- "prob" # added "prob" column which holds the records probability to be an outlier
}
results
}
