write.fstat<-function(dat,fname="genotypes.dat"){
# write a genotype data frame to a text file in fstat format
# dat is a dataframe with pop of origin in first column, and genotypes after
# fname is a string containing the name of the file to which the data set will be saved
nl<-dim(dat)[2]-1
if (!is.numeric(dat[,1])) dat[,1]<-as.integer(as.factor(dat[,1]))
np<-max(dat[,1])
dum<-getal.b(dat[,-1])
nu<-max(dum[,,1],dum[,,2],na.rm=TRUE)
nd<-3
if (nu < 10) nd<-1 else {if (nu<100 ) nd<-2}
write(c(np,nl,nu,nd),fname)
write(names(dat)[-1],fname,append=TRUE,ncolumns=1)
dat[is.na(dat)]<-0
write.table(dat,fname,append=TRUE,row.names=FALSE,col.names=FALSE,quote=FALSE,sep="\t")
}
#################################################################################
subsampind <- function(dat,sampsize=10){
  # subsample sampsize individuals from each population
  # dat is a data frame with pop of origin in first column and genotypes in the following ones
  # sampsize is the sample size, per pop of origin
  
  x<-table(dat[,1])
  if (min(x)<sampsize) warning(paste("At least one population contains [",min(x),"] individuals, less than the requested number [",sampsize,"].",
     "\n Populations with fewer than [",sampsize,"] individuals are returned unchanged",sep=" "))
  csx<-cumsum(x)
  ll<-c(1,csx[-length(csx)]+1) #
  rx<-cbind(ll,csx)
  dum<-which(rx[,2]-rx[,1]<sampsize)
  ldum<-length(dum)
  retainedb<-NULL
  if (ldum==0){
    retained<-apply(rx,1,fun<-function(y) sample(y[1]:y[2],replace=FALSE,size=sampsize))
  } else {
    retained<-apply(rx[-dum,],1,fun<-function(y) sample(y[1]:y[2],replace=FALSE,size=sampsize))
    if (ldum==1) {
      retainedb<-rx[dum,1]:rx[dum,2]
      } 
    else {
      retainedb<-apply(rx[dum,],1,fun<-function(y) (y[1]:y[2]))
      }
  }
  retain<-sort(c(unlist(retained),unlist(retainedb)))
  ddat<-dat[order(dat[,1]),]
  return(data.frame(ddat[retain,]))
}
################################
write.struct<-function(dat,ilab=NULL,pop=NULL,fname="dat.str"){
#write a file suitable to be read as an input for the program structure
dum<-getal.b(dat[,-1])
dum[is.na(dum)]<- -9
nind<-dim(dum)[1]
nloc<-dim(dum)[2]
popid<-dat[,1]
locnames<-paste("L",names(dat)[-1],sep="",collapse=" ")
write(locnames,fname)
for (i in 1:nind){
write(paste(popid[i],paste(dum[i,,1],collapse=" "),collapse=" "),fname,append=T)
write(paste(popid[i],paste(dum[i,,2],collapse=" "),collapse=" "),fname,append=T)
}
}

