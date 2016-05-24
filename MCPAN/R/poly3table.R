"poly3table"<-function(time, status, f, tumour=NULL, symbol="*")
{

N<-length(time)

if(N!=length(status))
{stop("time and status must be vectors of equal length")}

if(length(symbol)!=1)
  {stop("symbol must be a single value")}
 
statussign<-character(length=N)

if(is.null(tumour))
 {
  if(is.logical(status))
   {
   wtumour<-which(status==TRUE)
   statussign[wtumour]<-symbol
   statussign[-wtumour]<-""
   }
  else{

  if(!is.logical(status))
    {stop("status must be a logical vector;
      if not, please specify the value specifying presence of tumour
      in argument tumour")}
    }
 }
else
 {
 
  if(length(tumour)!=1)
   {stop("tumour must be a single value")}

  wtumor<-which(status==tumour)
  if(length(wtumor)==0)
   {warning("tumor was not found in status")}

  statussign[wtumor]<-symbol
  statussign[-wtumor]<-""
 }

outvec<-paste(time, statussign, sep="")

datlist<-split(outvec, f=as.factor(f))

shortendat<-function(x)
 {
  tab<-table(x)
  nam<-names(tab)
  ntab<-as.numeric(tab)
  wg1<-which(ntab>1)
  nchar<-as.character(ntab)
  nchar2<-paste("(", nchar,")", sep="")
  nchar2[-wg1]<-""
  nchar3<-paste(paste(nam, nchar2, sep=""), collapse=", ")
  return(nchar3)
 }

datlist2<-lapply(X=datlist, FUN=function(x){shortendat(x)})

return(datlist2)

}
