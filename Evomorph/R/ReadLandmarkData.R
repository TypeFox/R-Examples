ReadLandmarkData <- function(datafile, m_byscale=T){

if (!(file.exists(datafile))){stop("Could not find data file. Check name or path")}

dv<-readLines(datafile)

if (!(grepl("LM",dv[1]))){stop("Incorrect file format")}

if (!(grepl("[0-9]+",dv[1]))){stop("Incorrect file format")}

checkt<-grep("^LM=[0-9]+",dv,value=T)

checkt<-as.numeric(str_extract(checkt,"[0-9]+"))

if ((abs(max(checkt)-min(checkt)))){stop("Not all especimes have the same number of landmarks")}

num_landmark<-as.numeric(str_extract(dv[1],"[0-9]+"))
num_esp<-length(grep("^LM",dv))
x<-na.omit(as.numeric(str_extract(dv,"(^-?[:digit:]+[:punct:]*[:digit:]+)")))
x<-as.vector(x)
y<-na.omit(as.numeric(str_extract(dv," (-?[:digit:]+[:punct:]*[:digit:]+)")))
y<-as.vector(y)
if (length(x)!=num_esp*num_landmark || length(y)!=num_esp*num_landmark) {stop("Missing landmarks coordinates")}
scalev<-na.omit(str_extract(dv,"SCALE=([:digit:]+[:punct:]*[:digit:]+)"))
scalev<-as.numeric(str_extract(scalev,"([:digit:]+[:punct:]*[:digit:]+)"))
scalev<-as.vector(scalev)
if (length(scalev)!=num_esp){print ("Missing scale values")}
factorsv<-rep(1:num_esp,each=num_landmark)
sx<-split(x,ceiling(seq_along(x)/num_landmark))
sy<-split(y,ceiling(seq_along(y)/num_landmark))
if (m_byscale){for(i in 1:num_esp){sx[[i]]=sx[[i]]*scalev[i];sy[[i]]=sy[[i]]*scalev[i];}}
arrayc<-array(dim = c(num_landmark,2,num_esp))
for(i in 1:num_esp){
arrayc[,1,i]=sx[[i]]
arrayc[,2,i]=sy[[i]]
}
dfc<-data.frame(x = unlist(sx) ,y = unlist(sy), esp= factorsv)

lm<-list(df=dfc,array=arrayc,scale=scalev)

return(lm)

}