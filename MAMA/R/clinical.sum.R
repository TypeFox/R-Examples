#datname<-ColonData@datanames

clinical.sum<-function(x){
CL<-clinical(x)
datname<-datanames(x)
N=length(CL)
param<-sapply(CL,colnames)
param<-sort(unique(unlist(param)))

  S<-function(x, par){
  res<-NA
  c<-grep(par,names(x))
  if (length(c)!=0) res<-summary(x[,c]) else res<-summary(as.factor(rep(NA, nrow(x))))
  return(res)
  }
res<-list()
for (j in param) res[[j]]<-lapply(CL, S, j)

  C<-function(x, par){
  res<-NA
  c<-grep(par,names(x))
  if (length(c)!=0 & is.numeric(x[,c])) res<-x[,c] 
  return(res)
  }
con<-list()
for (j in param) con[[j]]<-lapply(CL, C, j)
con<-lapply(con, unlist)
CON<-unlist(lapply(con, function(x) !all(is.na(x))))
con.var<-names(con)[CON]

fill<-function(x,  datname){
  n<-length(x)
  categ<-unique(names(unlist(x)))
  M<-matrix(NA, ncol=length(categ), nrow=n+1)
  colnames(M)<-categ
  rownames(M)<-c(datname, "Total")
  for (i in 1:n) {
    w<-match( names(x[[i]]), categ )
    w<-w[!is.na(w)]
    M[i,w]<-x[[i]]
  }
  M[n+1,]<-apply(M[1:n,],2, function(x) sum(x, na.rm=T))
  M[is.na(M)]<-0
  #na.pos<-match("NA's",colnames(M)) 
  #if (length(na.pos)!=0){M[,ncol(M)]<-M[,na.pos]; M<-M[,-na.pos]} else M<-M[,-ncol(M)]
  return(M)
}

tab<-lapply(res, fill,datname)
if (length(con.var)>0 ) {dum<-nrow(tab[[con.var]]); tab[[con.var]][dum,]<-summary(con[[con.var]])}

tab.per<-lapply(tab, function(x) prop.table(as.table(x),1))
if (length(names(con)[CON])>0 ) tab.per[[names(con)[CON]]]<-NULL
T<-list(absolute=tab, relative=tab.per)
return(T)
}

