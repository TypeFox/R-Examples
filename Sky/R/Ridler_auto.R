Ridler_auto <-
function(path1,path2=TRUE,write=TRUE,pixel)
{

R<-list()
setwd(path1) #set up the directory where R needs to look for pictures
ll<-list.files(path = path1)  #list of files in the light folder (those could correpond to month, year.. so on)
for( j in 1:length(ll))   #Loop through the pictures within each of these folders
{
    if(path2==TRUE)
    {
    file<-ll[j]
    l<-list.files(path = file)
    eval(parse(text=paste("setwd('",path1,"/",file,"')",sep=""))) # enter the folder
        }
    if(path2==FALSE)
    {
        l<-ll
    }
    Result<-matrix(data=NA,nrow=length(l),ncol=2) #set up the matrix to store the results
    for(k in 1:length(l))   #Loop through the pictures within the folder
    {   r<-Ridler(readImage(l[k]),pixel,p=FALSE)
    	Result[k,1]<-l[k]
    	Result[k,2]<-r
    }
    setwd(path1)
    Result<-data.frame(Result)
    Result[,1]<-as.character(Result[,1])
    Result[,2]<-as.numeric(as.character(Result[,2]))
    colnames(Result)<-c("Name","Sky")
    R[[j]]<-Result
  if(write==TRUE & path2==TRUE)
    {
    eval(parse(text=paste("write.csv(Result,'Result_",file,".csv')",sep="")))
    names(R)[j]<-file
    }
     if(write==FALSE & path2==TRUE)
    {
    names(R)[j]<-file
    }
     if(write==TRUE & path2==FALSE)
    {
    write.csv(R[[1]],"Result.csv")
    }
 }
 if(path2==TRUE)
 {
    return(R)
 }
  if(path2==FALSE)
 {
    return(R[[1]])
 }
}
