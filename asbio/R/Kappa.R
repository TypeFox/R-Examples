Kappa<-function(class1,reference){
t<-table(class1,reference)

##################################
##assymetric matrix?

##if necessary, create more columns
tc<-match(colnames(t),rownames(t))
mtc<-matrix(ncol=ncol(t),nrow=length(tc[tc=="NA"]),0)
nrn<-colnames(t)[is.na(tc)==TRUE]
rownames(mtc)<-nrn
t1<-rbind(t,mtc)

##if necessary, create more rows
tr<-match(rownames(t1),colnames(t1))
mtr<-matrix(nrow=nrow(t1),ncol=length(tr[tr=="NA"]),0)
ncn<-rownames(t1)[is.na(tr)==TRUE]
colnames(mtr)<-ncn
t2<-cbind(t1,mtr)

##line up classes
sr<-sort(rownames(t2))
mr<-match(sr,rownames(t2))
t3<-t(t2[mr,])
sc<-sort(rownames(t3))
mc<-match(sc,rownames(t3))
t4<-t(t3[mc,])

###################################
agree<-diag(t4)
prod1<-apply(t4,1,sum)
prod2<-agree/prod1
user1<-apply(t4,2,sum)
user2<-agree/user1
N<-sum(t4)
k1<-sum(agree)
k2<-sum(prod1*user1)
khat<-((N*k1)-k2)/(N^2-k2)
result<-list()
result$ttl_agreement<-(sum(diag(t4))/sum(t4))*100
result$user_accuracy<-round(user2*100,1)
result$producer_accuracy<-round(prod2*100,1)
result$khat<-round(khat*100,1)
result$table<-t
result
}
