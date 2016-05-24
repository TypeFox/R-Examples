game.pred <-
function(res,home,away,neutral.site=FALSE){
if(neutral.site==FALSE&res$home.field==TRUE){
cat("Normal Distribution for Scores:\n")
if(is.null(res$n.ratings.offense)){
cat("N/A for this object.\n\n")}else{
cat(paste("Predicted score for ",home,": ",round(res$n.mean[match("LocationHome",names(res$n.mean))]+res$n.ratings.offense[match(home,names(res$n.ratings.offense))]-res$n.ratings.defense[match(away,names(res$n.ratings.defense))],2),"\n",sep=""))
cat(paste("Predicted score for ",away,": ",round(res$n.mean[match("LocationAway",names(res$n.mean))]+res$n.ratings.offense[match(away,names(res$n.ratings.offense))]-res$n.ratings.defense[match(home,names(res$n.ratings.defense))],2),"\n",sep=""))
cat("\n")
}
cat("Poisson Distribution for Scores:\n")
if(is.null(res$p.ratings.offense)){
cat("N/A for this object.\n\n")}else{
cat(paste("Predicted score for ",home,": ",round(exp(res$p.mean[match("LocationHome",names(res$p.mean))]+res$p.ratings.offense[match(home,names(res$p.ratings.offense))]-res$p.ratings.defense[match(away,names(res$p.ratings.defense))]),0),"\n",sep=""))
cat(paste("Predicted score for ",away,": ",round(exp(res$p.mean[match("LocationAway",names(res$p.mean))]+res$p.ratings.offense[match(away,names(res$p.ratings.offense))]-res$p.ratings.defense[match(home,names(res$p.ratings.defense))]),0),"\n",sep=""))
cat("\n")
}
cat("Binary Distribution for Outcomes:\n")
if(is.null(res$b.ratings)){
cat("N/A for this object.\n\n")}else{
cat(paste("Probability of ",home," defeating ",away,": ",round(pnorm(res$b.mean+res$b.ratings[match(home,names(res$b.ratings))]-res$b.ratings[match(away,names(res$b.ratings))]),3),"\n",sep=""))
cat("\n")
}
if(res$method=="NB.mov"){
cat("Normal Distribution for Margin of Victory:\n")
if(is.null(res$n.ratings.mov)){
cat("N/A for this object.\n\n")}else{
cat(paste("Predicted margin of victory for ",home," over ",away,": ",round(res$n.mean[match("LocationHome",names(res$n.mean))]+res$n.ratings.mov[match(home,names(res$n.ratings.mov))]-res$n.ratings.mov[match(away,names(res$n.ratings.mov))],2),"\n",sep=""))
cat("\n")
}
}
}else if(neutral.site==TRUE&res$home.field==TRUE){
cat("Normal Distribution for Scores:\n")
if(is.null(res$n.ratings.offense)){
cat("N/A for this object.\n\n")}else{
cat(paste("Predicted score for ",home,": ",round(res$n.mean[match("LocationNeutral Site",names(res$n.mean))]+res$n.ratings.offense[match(home,names(res$n.ratings.offense))]-res$n.ratings.defense[match(away,names(res$n.ratings.defense))],2),"\n",sep=""))
cat(paste("Predicted score for ",away,": ",round(res$n.mean[match("LocationNeutral Site",names(res$n.mean))]+res$n.ratings.offense[match(away,names(res$n.ratings.offense))]-res$n.ratings.defense[match(home,names(res$n.ratings.defense))],2),"\n",sep=""))
cat("\n")
}
cat("Poisson Distribution for Scores:\n")
if(is.null(res$p.ratings.offense)){
cat("N/A for this object.\n\n")}else{
cat(paste("Predicted score for ",home,": ",round(exp(res$p.mean[match("LocationNeutral Site",names(res$p.mean))]+res$p.ratings.offense[match(home,names(res$p.ratings.offense))]-res$p.ratings.defense[match(away,names(res$p.ratings.defense))]),0),"\n",sep=""))
cat(paste("Predicted score for ",away,": ",round(exp(res$p.mean[match("LocationNeutral Site",names(res$p.mean))]+res$p.ratings.offense[match(away,names(res$p.ratings.offense))]-res$p.ratings.defense[match(home,names(res$p.ratings.defense))]),0),"\n",sep=""))
cat("\n")
}
cat("Binary Distribution for Outcomes:\n")
if(is.null(res$b.ratings)){
cat("N/A for this object.\n\n")}else{
cat(paste("Probability of ",home," defeating ",away,": ",round(pnorm(res$b.ratings[match(home,names(res$b.ratings))]-res$b.ratings[match(away,names(res$b.ratings))]),3),"\n",sep=""))
cat("\n")
}
cat("Normal Distribution for Margin of Victory:\n")
if(is.null(res$n.ratings.mov)){
cat("N/A for this object.\n\n")}else{
cat(paste("Predicted margin of victory for ",home," over ",away,": ",round(res$n.ratings.mov[match(home,names(res$n.ratings.mov))]-res$n.ratings.mov[match(away,names(res$n.ratings.mov))],2),"\n",sep=""))
cat("\n")
}
}else{
 cat("Normal Distribution for Scores:\n")
if(is.null(res$n.ratings.offense)){
cat("N/A for this object.\n\n")}else{
cat(paste("Predicted score for ",home,": ",round(res$n.mean[1]+res$n.ratings.offense[match(home,names(res$n.ratings.offense))]-res$n.ratings.defense[match(away,names(res$n.ratings.defense))],2),"\n",sep=""))
cat(paste("Predicted score for ",away,": ",round(res$n.mean[1]+res$n.ratings.offense[match(away,names(res$n.ratings.offense))]-res$n.ratings.defense[match(home,names(res$n.ratings.defense))],2),"\n",sep=""))
cat("\n")
}
cat("Poisson Distribution for Scores:\n")
if(is.null(res$p.ratings.offense)){
cat("N/A for this object.\n\n")}else{
cat(paste("Predicted score for ",home,": ",round(exp(res$p.mean[1]+res$p.ratings.offense[match(home,names(res$p.ratings.offense))]-res$p.ratings.defense[match(away,names(res$p.ratings.defense))]),0),"\n",sep=""))
cat(paste("Predicted score for ",away,": ",round(exp(res$p.mean[1]+res$p.ratings.offense[match(away,names(res$p.ratings.offense))]-res$p.ratings.defense[match(home,names(res$p.ratings.defense))]),0),"\n",sep=""))
cat("\n")
}
cat("Binary Distribution for Outcomes:\n")
if(is.null(res$b.ratings)){
cat("N/A for this object.\n\n")}else{
cat(paste("Probability of ",home," defeating ",away,": ",round(pnorm(res$b.ratings[match(home,names(res$b.ratings))]-res$b.ratings[match(away,names(res$b.ratings))]),3),"\n",sep=""))
cat("\n")
}
cat("Normal Distribution for Margin of Victory:\n")
if(is.null(res$n.ratings.mov)){
cat("N/A for this object.\n\n")}else{
cat(paste("Predicted margin of victory for ",home," over ",away,": ",round(res$n.ratings.mov[match(home,names(res$n.ratings.mov))]-res$n.ratings.mov[match(away,names(res$n.ratings.mov))],2),"\n",sep=""))
cat("\n")
}

}
}
