sim.genetree <- function(n,numbsim){
statistics<-vector()
for (i in 1:numbsim){
	statisticstemp<-c(1:4)*0
	lineages<-c(1:n)*0+1
	while (length(lineages)>1){
			merge<-sort(sample(1:length(lineages),2))
			statisticstemp<-stats(lineages,merge[1],merge[2],statisticstemp)
			lineages[merge[1]]<-lineages[merge[1]]+lineages[merge[2]]
			lineages<-lineages[-merge[2]]	
	}
	statistics<-rbind(statistics, statisticstemp)
}
rownames(statistics)<-c(1:length(statistics[,1]))
colnames(statistics)<-c("Colless","s","Sackin","cherries")
statistics
}