Heatmap <- function(Dataset, Id, Outcome, Time, ...){

Id <- Dataset[,paste(substitute(Id))]
Outcome <- Dataset[,paste(substitute(Outcome))]
Time <- Dataset[,paste(substitute(Time))]
  
Data <- cbind(Id, Outcome, Time)

Dataset.Wide <-
  data.frame(cbind(Data[,1], Data[,2], Data[,3]))
colnames(Dataset.Wide) <- c("Id", "Outcome", "Time")
Data_Wide <- 
  reshape(data = Dataset.Wide, timevar = "Time", idvar="Id", direction="wide")

Data <- 
  data.frame(Data)
min_2_vals <- max(as.numeric(names((table((Data$Time))[table((Data$Time))>=2]))))
Ylab <- Xlab <- 
  unique(Data$Time)[unique(Data$Time)[1]:min_2_vals]

correl <- cor(Data_Wide, use = "pairwise.complete.obs")
correl <- 
  correl[2: min_2_vals, 2: min_2_vals]
row <- dim(correl)[1]
correl <- (matrix(correl, nrow = row))

psych::cor.plot(correl, labels=1:row, show.legend = TRUE, ...)
}
