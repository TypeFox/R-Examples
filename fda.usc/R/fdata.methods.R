######################################################
######################################################
Math.fdata<-function (x,...) {
   if (!inherits(x, "fdata")) stop("Objects must be of class fdata")
   x$data<-callGeneric(x$data@.Data,...)
   x
}
######################################################
######################################################
Ops.fdata<-function (e1, e2 = NULL) {
inhe.fdata1<- inherits(e1, "fdata")
inhe.fdata2<- inherits(e2, "fdata")
if (!inhe.fdata1 && !inhe.fdata2)
  stop("Neither argument are fdata object .")
if (inhe.fdata1 && inhe.fdata2) {
   fdataobj<-e1
   if (!all(e1[["argvals"]] == e2[["argvals"]]))  stop("Error in argvals")
   if (all(dim(e1[["data"]])==dim(e2[["data"]]))) {
     fdataobj$data<-callGeneric(e1$data@.Data, e2$data@.Data)   }
   else {
     if (nrow(e2[["data"]])==1) {
           fdataobj$data<-sweep(e1$data, 2, e2$data, FUN = .Generic)
     }
     if (nrow(e1[["data"]])==1) {
           fdataobj$data<-sweep(e2$data, 2, e1$data, FUN = .Generic)
     } else    stop("Error in data dimenstion")
   }
   return(fdataobj)
   }
else {
if (!inhe.fdata1 && inhe.fdata2) {
 if (all(dim(e1)==dim(e2[["data"]]))) {
  fdataobj<-e2
  fdataobj$data<-callGeneric(e1@.Data, e2$data@.Data)
  }
 else {if (nrow(e2[["data"]])==1) {
           dat<-sweep(e1$data, 2, e2$data, FUN = .Generic)
     }}
 return(fdataobj)
}
if (inhe.fdata1 && !inhe.fdata2) {
 if (all(dim(e1)==dim(e2[["data"]]))) {
   fdataobj<-e1
   fdataobj$data<-callGeneric(e1$data@.Data,e2@.Data)
  }
 else {if (nrow(e2[["data"]])==1) {
           dat<-sweep(e1$data, 2, e2$data, FUN = .Generic)
     }}
 return(fdataobj)
}}
}

# sweep(e1$data, 2, e2$data, FUN = "-")
######################################################
######################################################
Summary.fdata<-function(...,na.rm = FALSE) {
   aa<-match.call()
   fdataobj<-aa[[2]]
   if (!inherits(fdataobj, "fdata")) stop("Objects must be of class fdata")
   do.call(.Generic,list(fdataobj$data,na.rm=na.rm))
}

######################################################
######################################################



