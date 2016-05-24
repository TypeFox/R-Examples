linerulecoalitional <-
function(n=NA,a=NA,av=NA,d=NA,h=NA,m=NA){
STI<-STIcoo(n,a,av,d,h,m)[[1]]
STI<-STI[,n+2]
coa<-cbind(as.matrix(coalitions(n)[[1]]),STI)
for (i in 2:nrow(coa)){
aux<-which(coa[i,1:n]==1)
aux1<-linerule(length(aux),a,av[aux],d[aux],h[aux],m=NA)
for (j in 1:length(aux)){
coa[i,aux[j]]<-aux1[j]
}
}
colnames(coa)<-c(1:n,"Coalitional cost")
rownames(coa)<-rep(" ",2^n)
return(coa)}
