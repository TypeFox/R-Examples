### Name: read.results
### Title: Reads in entire set of WASIM-simulation results that is specified in data.types
### Aliases: read.results
### Keywords: ~load ~read ~import ~WASIM data 

### ** Examples

##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--    or do  help(data=index)  for the standard data sets.

## The function is currently defined as
`read.results` <-
function (record, path="wasim/output",  ts.length,
subcatchments,
data.types.prefix = data.types$prefix,
data.types.has_stat = data.types$has_stat,
ending,
endings=rep(ending,NROW(data.types.prefix)),
na.values=c("999", "999.00"),
read.files.nr=1:nrow(data.types),
remove.missing=FALSE
)
{
t.data<-array(dim=c(subcatchments+1,NROW(data.types),ts.length))
    
toberemoved=NULL

for (i in read.files.nr ){
 name <- paste(sep="", path, record,"/",data.types$prefix[i],endings[i])
 cat("i: ", i, " ", name)
 if (!file.exists(name))
 {
   cat("...not found\n")
   toberemoved=c(toberemoved,i)   #mark this entry for later removal
 } else {
   cat("\n")

   #t.data[,i,] <- myRead.table(name, subcatchments, has_stat=data.types.has_stat[i], ts.length, na.values) #load file
   try(t.data[1:(subcatchments+data.types.has_stat[i]),i,] <- myRead.table(name, subcatchments, has_stat=data.types.has_stat[i], ts.length, na.values)) #load file
   
 }
}

if (remove.missing & length(toberemoved)>0)
{
  assign("data.types",data.types[-toberemoved,],inherits=TRUE)     #remove entry, if file was not found
  t.data=t.data[,-toberemoved,]
  cat(paste("\nremoved",length(toberemoved),"entries from data list where data files are missing.\n"))
}
 
return(t.data)
}

