fdata=function(mdata,argvals=NULL,rangeval=NULL,names=NULL,fdata2d=FALSE){
call<-match.call()
if (length(call[[2]])>1) nam<-"fdataobj"
else nam<-call[[2]]
if (is.array(mdata))  {
    dm<-dim(mdata)
    len.dm<-length(dm)
    if (len.dm<3) fdata2d=FALSE
    else {
      fdata2d=TRUE
      if (len.dm>3) stop("Not implemented yet for dim(mdata)>",3)
    }
}
if (is.list(argvals)) fdata2d=TRUE
if (is.list(rangeval)) fdata2d=TRUE

if (fdata2d){
out=list("data"=NULL)
if (length(class(mdata))>1) class(mdata)<-class(mdata)[1]
out<-switch(class(mdata),
matrix={
        out[["data"]]<-array(NA,dim=c(1,dm[1],dm[2]))
        out[["data"]][1,,]<-mdata
                out},
data.frame={
        out[["data"]]<-array(NA,dim=c(1,dm[1],dm[2]))
        out[["data"]][1,,]<-mdata
        out},
fdata=stop("The data is a fdata class object"),
array={
        out[["data"]]=mdata
        out}
)
dm<-dim(out[["data"]] )
len.dm<-length(dm)
if (is.null(argvals)) {
 argvals<-list()
 if (len.dm>2){
  len.argvals<-len.dm-1
  for (i in 2:len.dm) {
    if (is.null(dimnames(out[["data"]][i])))  argvals[[i-1]]<-1:dm[i]
  }
 }
else {
  len.argvals<-len.dm
  for (i in 1:len.dm) {
  if (is.null(dimnames(out[["data"]][i])))  argvals[[i]]<-1:dm[i]
# nam<-as.character(call[[2]])
# print(length(call[[2]]))
 if (is.null(names(argvals[[i]]))) names(argvals[[i]])<-paste(nam,1:dm[i],sep="")
  }
 }
 out[["argvals"]]<-argvals
}
else {
  #verificar dimensiones
  if (is.list(argvals)) {
     len.argvals<-length(argvals)
#     print(len.argvals)
     if (len.dm>2){
       for (i in 1:len.argvals) {
         if (length(argvals[[i]])!=dm[i+1]) stop("Incorrect dimension in between mdata and argvals arguments")
    }   }
    else {
       for (i in 1:len.argvals) {
        if (length(argvals[[i]])!=dm[i]) stop("Incorrect dimension in between mdata and argvals arguments")
    }  }
  }
else stop("The argument argvals must be a list")
#print("peta name")
  if (is.null(names(argvals))) names(argvals)<-paste(drop(nam),1:len.argvals,sep="")
  out[["argvals"]]<-argvals
}
##########################
if (is.null(rangeval))  {
 rangeval<-list()
 for (i in 1:len.argvals) {
    rangeval[[i]]<-range(argvals[i])
  }
 if (is.null(names(rangeval))) names(rangeval)<-names(argvals)
 out[["rangeval"]]<-rangeval
 len.rangeval<-length(out[["rangeval"]])
}
else {
  if (is.list(rangeval)) {
     len.rangeval<-length(rangeval)
     if (len.rangeval!=len.argvals) stop("Incorrect dimension in rangeval argument")
  }
else stop("The argument reangeval must be a list")
  if (is.null(names(rangeval))) names(rangeval)<-names(argvals)
  out[["rangeval"]]<-rangeval
}
if (is.null(names))  {
 names<-list()
 names[["xlab"]]<-paste("argvals ",1:len.argvals,sep="")
 names[["ylab"]]<-paste("values of ",nam,len.argvals,sep="")
 names[["main"]]<-paste(nam,len.argvals,sep="")
 out[["names"]]<-names
}
else {
  if (is.list(names)) {
     len.names<-length(names[["xlab"]])
    if (len.rangeval!=len.names) stop("Incorrect dimension in names argument")
  }
else stop("The argument names must be a list")
 out[["names"]]<-names
}
if (is.null(dimnames(out$data))) dimnames(out$data)<-list("data"=1:dm[1],"argvals.x"=argvals[[1]],"argvals.y"=argvals[[2]])
class(out)=c("fdata","fdata2d")
return(out)
}   #### fdata1d
else{
out=list("data"=NULL)
if (length(class(mdata))>1) class(mdata)<-class(mdata)[1]
out<-switch(class(mdata),
matrix={
        out[["data"]]=mdata
        out},
data.frame={
            out[["data"]]=as.matrix(mdata)
            out},
fdata=stop("The data could not be converted into fdata class"),
numeric={
         out[["data"]]=matrix(mdata,nrow=1)
         out},
integer={
         out[["data"]]=matrix(mdata,nrow=1)
         out},
fd={
   r = mdata$basis[[3]]
   if (is.null(argvals))
     argvals= seq(r[1], r[2], len =mdata$basis$nbasis)
    nb<- length(argvals)
    tt = argvals
#   tt = seq(r[1], r[2], len = length(mdata$fdnames$time))

   out[["data"]] = t(eval.fd(tt, mdata))
   if (!is.null(mdata$fdnames$reps)) rownames(out[["data"]]) = mdata$fdnames$reps
   else      rownames(out[["data"]]) =1:nrow( out[["data"]])
   if (!is.null(mdata$fdnames$time)) {
      colnames(out[["data"]]) = 1:ncol( out[["data"]])
      }
   else      {   colnames(out[["data"]]) =1:ncol( out[["data"]])   }
   out
   },
fds={
 out[["data"]] =t(mdata$y)
 if (is.null(mdata$x))       out[["argvals"]] =1:ncol(out[["data"]])
 else { 
   if (is.numeric(mdata$x)) out[["argvals"]]=mdata$x 
   else out[["argvals"]]=seq(mdata$x[1], mdata$x[length(mdata$x)],
  len = length(mdata$x))
  }
  out[["rangeval"]]<-range(out[["argvals"]])  
  out[["names"]]<-list("main"=deparse(nam),"xlab"=mdata$xname,"ylab"=mdata$yname)
  class(out)<-"fdata"
  return(out)
  },
fts={
 out[["data"]] =t(mdata$y)
 if (is.null(mdata$x))       out[["argvals"]] =1:ncol(out[["data"]])
 else { 
   if (is.numeric(mdata$x)) out[["argvals"]]=mdata$x 
   else out[["argvals"]]=seq(mdata$x[1], mdata$x[length(mdata$x)],
  len = length(mdata$x))
  }
  out[["rangeval"]]<-range(out[["argvals"]])  
  out[["names"]]<-list("main"=deparse(nam),"xlab"=mdata$xname,"ylab"=mdata$yname)
  class(out)<-"fdata"
  return(out)
  },
sfts={
 out[["data"]] =t(mdata$y)
 if (is.null(mdata$x))       out[["argvals"]] =1:ncol(out[["data"]])
 else { 
   if (is.numeric(mdata$x)) out[["argvals"]]=mdata$x 
   else out[["argvals"]]=seq(mdata$x[1], mdata$x[length(mdata$x)],
  len = length(mdata$x))
  }
  out[["rangeval"]]<-range(out[["argvals"]])  
  out[["names"]]<-list("main"=deparse(nam),"xlab"=mdata$xname,"ylab"=mdata$yname)
  class(out)<-"fdata"
  return(out)
  }
)
nc<-nrow(out[["data"]])
np<-ncol(out[["data"]])
if (is.null(argvals)) {
 if (is.null(colnames(out[["data"]]))) {out[["argvals"]]=1:ncol(out[["data"]])}
   else {
#   if (!any(is.na(as.numeric(colnames(out[["data"]]))))) {
#    out[["argvals"]]=as.numeric(colnames(out[["data"]]))   }
#    else    out[["argvals"]]=1:ncol(out[["data"]])}   }
out[["argvals"]]=1:ncol(out[["data"]])}   }
else     out[["argvals"]]=argvals
lentt=length(out[["argvals"]])
if (is.null(rangeval)) rangeval=range(out[["argvals"]])
out[["rangeval"]]<-rangeval
if ((np!=lentt) && (nc==lentt)) {
         out[["data"]]=matrix(out[["data"]],ncol=nc)
         nc<-1
         print("Warning: The transposed data is returned")      }
else    out[["data"]]=out[["data"]]
if (is.null(dimnames(mdata))) {
#rownames(out[["data"]])<-1:nc
colnames(out[["data"]])=round(out[["argvals"]],4)
}
out[["names"]]<-list("main"="fdataobj","xlab"="t","ylab"="X(t)")
if (!is.null(names$main)) out$names$main<-names$main
if (!is.null(names$xlab)) out$names$xlab<-names$xlab
if (!is.null(names$ylab)) out$names$ylab<-names$ylab
class(out)="fdata"
return(out)
}
}