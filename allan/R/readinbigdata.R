readinbigdata <-
function(filename, chunksize,...){
#This function sets a connection to the large dataset with the specified chunksize.
#The return value is either the next chunk of data or NULL if there is no additional data left.
#Additionally if a reset=TRUE flag is passed, then the data stream goes back to the beginning.
#This was originally done to accommodate the bigglm data function option
#Taken mostly from help from biglm package

#initialize connection
conn<-NULL
function(reset){
if(reset){
if(!is.null(conn)) close(conn)
conn<<-file (description=filename, open="r")
#print("new connection open")
} else{
#make sure header isn't read for other cases and assign next block of data
rval<-read.csv(conn, nrows=chunksize,header=FALSE,skip=1,...)

#print("rows processed")
#print(dim(rval))

if (nrow(rval)==0) {
close(conn)
conn<<-NULL
rval<-NULL
#print("end of file reached")
}
#print(reset)
#print(dim(rval))
return(rval)
}
}
}

