getbestchunksize <-
function(filename,MemoryAllowed=0.5,TestedRows=1000,AdjFactor=0.095,silent=TRUE){
#Function that tests data size and adjusts memory for best chunking of large dataset
#This is done by reading in a number of rows(1000 by default)and then measuring the size of the memory
#used.  Memory allwed is specified in Gb.  The adjfactor is a factor used to adjust memory for overhead
#in the biglm fitting functions.
  
#get column names
columnnames<-names(read.csv(filename, nrows=2,header=TRUE))
#read in rows and test size
datapreview<-read.csv(filename, nrows=TestedRows,header=TRUE)
datamemsize<-object.size(datapreview)
optimalchunksize=floor(((MemoryAllowed*1000000000)/datamemsize[1])*TestedRows *AdjFactor)
if (silent!=TRUE){
print("Total memory usage for 1000 lines:")
print(datamemsize)
print("Chunksize for dataframe after adjustment factor:")
print(optimalchunksize)
}
return(optimalchunksize)
}

