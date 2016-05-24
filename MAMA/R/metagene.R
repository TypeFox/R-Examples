metagene<-function(x, results)
{
metagene<-list()
for (i in 1:length(results))
{
if (x %in% rownames(results[[i]])) metagene[[i]]<-results[[i]][x,]
if ("METRADISC" %in% class(results[[i]])) {
  metagene[[i]] <- c(results[[i]]$RQ[x,], results[[i]]$MCtest[x,]) 
}
if ("ES.GeneMeta" %in% class(results[[i]])  )
  metagene[[i]] <- c(results[[i]]$theScores[x,], results[[i]]$ScoresFDR$two.sided[x,])
class(metagene[[i]])<- c(class(results[[i]]),paste(class(results[[i]]),".gene", sep=""))
}
names(metagene)<-names(results)
class(metagene)<-"metagene"
return(metagene)
}