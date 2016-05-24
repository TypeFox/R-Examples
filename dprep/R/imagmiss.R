imagmiss <-
function(data,name="")
{

#Function to create a graph of the observations of the dataset
#leaving white lines where data is missing.
#The main idea is to use the original dataset to create 
#a temporary dataset containing 1 if a value is found or
#0 is the value is missing. The temporary data set is graphed by column
#changing color for each feature and leaving a blank horizontal line if
#value is missing.

#Uses the R function image from the base library.

#data:  the dataset
#name:the name of the dataset as desired in the graph title

  ncol=dim(data)[2]
  nrow=dim(data)[1]

  xaxis=colnames(data)
    xaxis = xaxis[-ncol]
    ticks = 1:(dim(data)[2] - 1)
    data = as.matrix(data)
    data = data[, -ncol]
    cat("Report on missing values for ",name,":\n")
    cat("\nNumber of missing values overall: ")
    cat(sum(is.na(data)))
    cat("\nPercent of missing values overall: ")
   cat((sum(is.na(data))/(dim(data)[1]*dim(data)[2]))*100)
   cat("\nFeatures with missing values (percent): \n")
   print(colSums(is.na(data))[colSums(is.na(data))!=0]/dim(data)[1]*100)
   #cat("\n",length(which(colSums(is.na(data))!=0)))
   cat("\nPercent of features with missing values: ")
   cat(length(which(colSums(is.na(data))!=0))/dim(data)[2]*100)
   cat("\nNumber of instances with missing values: ")
   cat(length(which(rowSums(is.na(data))!=0, arr.ind=TRUE)))
   cat("\nPercent of instances with missing values: ")
   cat(length(which(rowSums(is.na(data))!=0))/dim(data)[1]*100)
   cat("\n")
  data[which(data!="NA")]=1
  data[-which(data!="NA")]=0

  ncol1=ncol-1

  for(i in 1:ncol1)
   {
     data[data[,i]!=0,i]=i
   }

  x=1:ncol1
  y=1:nrow
data=apply(data,1,as.numeric)
  graph.title=paste("Distribution of missing values by variable for - ",name)
 
  image(x,y,data,col=c(0,topo.colors(100)),xlab="features",ylab="instances",axes=FALSE,main=(graph.title),cex.main=.7)
  axis(1,labels=xaxis,at=ticks)
}
