object.extract <-
function(func,variable=TRUE,subset.q=TRUE,exception=NULL){
alph<-c(letters,toupper(letters))
num<-c(0:9,".")
symb<-c("(",")","|","+","-","*","/","~","^",">","<","=",",",":","%")
sqb<-c("[","]")
if(subset.q==TRUE){sqbX<-sqb}else{sqbX<-NULL}
#split string and indexing
td1<-c("|",strsplit(func,"")[[1]],"|")
names(td1)<-seq(1,length(td1),1)
i<-1;j1<-1;V.List<-F.List<-Variable<-Function<-NULL
repeat{
   if(td1[i]%in%alph){
       j1<-i
       repeat{
       if(subset.q==TRUE & td1[j1+1]=="["){ j2<-j1; repeat{if(td1[j2+1]=="]") break; j2<-j2+1};  V.List<-c(V.List,paste(td1[i:(j2+1)],collapse="",sep="")); j1<-j2+1}
       if(td1[j1+1]%in%c(symb,"[","]")) break
       j1<-j1+1}
       if(td1[j1+1]=="("){Function<-paste(td1[i:j1],collapse="",sep="")}else{if(subset.q==TRUE & td1[j1+1]=="]"){Variable<-paste(td1[i:(j1+1)],collapse="",sep="")}else{Variable<-paste(td1[i:j1],collapse="",sep="")}    }
       V.List<-c(V.List,Variable);F.List<-c(F.List,Function)
    i<-j1+1
    }else{i<-i+1}
if(i==length(td1)) break
}
V.List<-unique(V.List);V.List<-setdiff(V.List,exception)
F.List<-unique(F.List);F.List<-setdiff(F.List,exception)
if(variable==TRUE){return(V.List)}else{return(F.List)}}
