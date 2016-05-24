drug.trial.table<-c(120,60,10,10)
names(drug.trial.table)<-c("Marked","Moderate", "Slight","Stationary")
drug.trial.table
barplot(drug.trial.table,main="Responses to new drug",xlab="response",ylab="frequency")
freq1way(drug.trial.table,c(0.5,0.25,0.1,0.15))

