summarybest.model <- function(x,N=NULL){
if(is.null(N)) N <- x$top
temp <-  sapply(x$BestModels,FUN=function(x){res <- x[1:N]})
name <- c("Rank","nVisits","FirstVisit","nEvalBefore1st","ModeSize","logCondPost","postProb","jeffries")
res <- temp[,name]
res1 <- matrix(as.numeric(res),nrow=N,ncol=length(name),byrow=FALSE)
result <- data.frame(round(res1[,1:5],0),round(res1[,6],2),signif(res1[,7],digits=3),round(res1[,8],2),temp[,"modelName"])
names(result) <- c(name,"modelName") 
names(result)[which(names(result)=="jeffries")] <- "jeffrey"
return(list(result=result,N=N))
}
