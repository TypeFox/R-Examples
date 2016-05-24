perm.fact.test<-function(Y,X1,X2,X3=NA, perm=100, method = "a"){

if(all(is.na(X3))){init.model<-anova(lm(Y~X1*X2))}
if(all(!is.na(X3))){init.model<-anova(lm(Y~X1*X2*X3))}


r<-length((init.model)$"F value")-1
F.init<-init.model$"F value"[1:r]
MS <- init.model$"Mean Sq"


# (a)
if(method == "a"){
F.perm<-matrix(nrow=r,ncol=perm)
if(all(is.na(X3))){
for(i in 1:perm){
Y.new<-sample(Y,replace=FALSE)
F.perm[,i]<-anova(lm(Y.new~X1*X2))$"F value"[1:r]}}

if(all(!is.na(X3))){
for(i in 1:perm){
Y.new<-sample(Y,replace=FALSE)
F.perm[,i]<-anova(lm(Y.new~X1*X2*X3))$"F value"[1:r]}}
}


# (b)
if(method == "b"){
MS.perm<-matrix(nrow=r + 1,ncol=perm)
if(all(is.na(X3))){
for(i in 1:perm){
Y.new<-sample(Y,replace=FALSE)
MS.perm[,i]<-anova(lm(Y.new~X1*X2))$"Mean Sq"}}

if(all(!is.na(X3))){
for(i in 1:perm){
Y.new<-sample(Y,replace=FALSE)
MS.perm[,i]<-anova(lm(Y.new~X1*X2*X3))$"Mean Sq"}}

F.perm <- matrix(nrow = r, ncol = perm)
for(i in 1:perm){
F.perm[,i] <- MS.perm[,i][1:r]/MS.perm[,i][r+1]
}
}

p1<-matrix(nrow=r,ncol=perm)
for(i in 1:r){
p1[i,]<-F.perm[i,]>=F.init[i]}

p2<-apply(p1,1,function(x){length(x[x==TRUE])})

p.val<-(p2+1)/perm
p.val<-ifelse(p.val>1,1,p.val)

if(all(is.na(X3))){
Table<-data.frame(Initial.F=init.model$"F value",Df = init.model$"Df",row.names=c("X1","X2","X1:X2","Residual"),pval=c(p.val,NA))}
if(all(!is.na(X3))){
Table<-data.frame(Initial.F=init.model$"F value",Df = init.model$"Df",row.names=c("X1","X2","X3","X1:X2","X1:X3","X2:X3","X1:X2:X3","Residual"),pval=c(p.val,NA))}
res<-list()
res$Table<-Table
res
}