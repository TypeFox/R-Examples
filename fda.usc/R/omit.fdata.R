#################################################################
#################################################################
omit.fdata<-function(fdataobj,y=NULL){
nas<-apply(fdataobj$data,1,count.na)
if (!is.null(y)) {
nas.g<-is.na(y)
if (is.null(names(y))) names(y)<-1:length(y)
if (any(nas) & !any(nas.g)) {
   bb<-!nas
   cat("Warning: ",sum(nas)," curves with NA are omited\n")
   fdataobj$data<-fdataobj$data[bb,]
   y<-y[bb]
   }
else {
if (!any(nas) & any(nas.g)) {
   cat("Warning: ",sum(nas.g)," values of group with NA are omited \n")
   bb<-!nas.g
   fdataobj$data<-fdataobj$data[bb,]
   y<-y[bb]
   }
else {
if (any(nas) & any(nas.g))  {
   bb<-!nas & !nas.g
   cat("Warning: ",sum(!bb)," curves  and values of group with NA are omited \n")
   fdataobj$data<-fdataobj$data[bb,]
   y<-y[bb]
   }
}}}
else {
if (any(nas)) {
   bb<-!nas
   cat("Warning: ",sum(nas)," curves with NA are omited\n")
   fdataobj$data<-fdataobj$data[bb,]
   }}
return(list(fdataobj,y))
}
#################################################################
#################################################################

omit2.fdata<-function(fdataobj,index.na=FALSE){
if (is.fdata(fdataobj)) {
   nas<-apply(fdataobj$data,1,count.na)
   ind.na<-which(!nas)
   fdataobj<-fdataobj[ind.na]
   if (index.na) return(list("fdataobj"=fdataobj,"index.na"=ind.na))
   else   return(fdataobj)
   }
else {
#names.list<-names(is.list(fdataobj))
if (is.list(fdataobj)) {
 ind.na2<-list()
 list.na<-list()
 ind.na3<-NULL
 n2<-NULL
 nobj<-length(fdataobj)
 if (is.null(names(fdataobj))) namesobj<-1:nobj
 else namesobj<-names(fdataobj)
 for (i in 1:nobj) {
 list.na[[i]]<-NULL
 x<-fdataobj[[i]]
 ind.na2[[i]]<-switch(class(x),
     "fdata"={     n<-nrow(x); nas<-apply(x$data,1,count.na)},
     "matrix"={    n<-nrow(x); nas<-apply(x,1,count.na)},
   "data.frame"={  n<-nrow(x); nas<-apply(x,1,count.na)},
   "vector"={    n<-length(x); nas<-is.na(x)},
   "factor"={    n<-length(x); nas<-is.na(x)},
   "integer"={   n<-length(x); nas<-is.na(x) },
   "numeric"={   n<-length(x); nas<-is.na(x) }  )
  if (sum(ind.na2[[i]])>0) {
     list.na[[i]]<-which(ind.na2[[i]])
     ind.na3<-union(ind.na3,list.na[[i]])
     }
  else list.na[[i]]<-NULL
  if (is.null(n2)) n2<-n
  else if (n2!=n) stop("The  elements  of the list have incorrect dimensions")
   }
  for (i in 1:length(fdataobj)) {
   fdataobj[[i]]<-switch(class(fdataobj[[i]]),
   "fdata"={   fdataobj[[i]][-ind.na3]},
   "matrix"={       fdataobj[[i]][-ind.na3,]},
   "data.frame"={       fdataobj[[i]][-ind.na3,]},
   "vector"={   fdataobj[[i]][-ind.na3]},
   "factor"={     fdataobj[[i]][-ind.na3]},
   "integer"={     fdataobj[[i]][-ind.na3]},
   "numeric"={     fdataobj[[i]][-ind.na3]}                )
  }
  names(fdataobj)<-namesobj
  if (index.na) return(list("fdataobj"=fdataobj,"index.na"=list.na,"ind.na"=ind.na3))
  else   return(fdataobj)
}
else {stop("No fdata or list class object introduced")}
}
}
#################################################################
#################################################################
missing.fdata<-function(fdataobj,basis=NULL){
   tt<-fdataobj$argvals
   rtt<-fdataobj$rangeval
   n<-nrow(fdataobj) 
   np<-length(tt)

   nas <- apply(fdataobj$data, 1, count.na)
   if (any(nas))  cat("Warning: ", sum(nas), " curves with NA are omited\n")
   nas <- which(nas)
   xall<-fdataobj
#   if (is.null(basis)) basis<-create.bspline.basis(rangeval = rtt, nbasis = max(5,floor(tt/4)))
   if (is.null(basis)) {
#      basis<-create.bspline.basis(rangeval = rtt, nbasis = length(tt))
      basis<-create.bspline.basis(rangeval = rtt, nbasis  = max(5,floor(tt/4)))
                       }
   for (i in nas) {
#  	cat(i)
	# is.na(xna$data[i,])
	curve<- which(!is.na(fdataobj$data[i,]))
	xall$data[i,-curve]<-eval.fd(tt,Data2fd(tt[curve],fdataobj$data[i,curve],basis))[-curve,1]
	}
xall
}

