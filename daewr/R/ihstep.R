ihstep <- function(y,des) {
#get names of linear terms
lin<-colnames(des)
#make replacement table for use in getting quadratic names
values<-c("A","B","C","D","E","F","G","H","I","J","K","L","M","N","O","P","Q","R","S","T","U","V","W","X","Y","Z")
repl<-c("I(A^2)","I(B^2)","I(C^2)","I(D^2)","I(E^2)","I(F^2)","I(G^2)","I(H^2)","I(I^2)","I(J^2)","I(K^2)","I(L^2)","I(M^2)","I(N^2)","I(O^2)","I(P^2)","I(Q^2)","I(R^2)","I(S^2)","I(T^2)","I(U^2)","I(V^2)","I(W^2)","I(X^2)","I(Y^2)","I(Z^2)")
repl.tab <- cbind(values, repl)
#get quadratic variable names
indx <- match(lin, repl.tab[, 1], nomatch = 0)
quad<-lin
quad[indx != 0] <- repl.tab[indx, 2]
quad<-paste(quad,collapse='+')
# create the data frame for analysis
dat<-data.frame(y=y,des)
#gets the model matrix
lm1<-lm(y~(.)^2,data=dat)
mm<-model.matrix(lm1)
fact<-colnames(mm)
fact<-fact[-1]
fact<-paste(fact,collapse='+')
mod<-paste(c(fact,quad),collapse='+')
lmtest<-lm(y~A+B+A:B+I(A^2)+I(B^2),data=dat)
lm2<-lm(reformulate(termlabels=mod, response='y'),data=dat) # This works
mm<-model.matrix(lm2)
#deletes the constant column from the model matrix
mm<-mm[,2:ncol(mm)]
#computes absolute correlation of each column in model matrix with response 
cm<-abs(cor(mm,(dat$y)-mean(dat$y)))
# gets the index of the maximum correlation
idm<-which.max(cm)
#get vector of model terms
varname<-unlist(dimnames(cm))
term<-varname[idm]
#checks to see if term is an interaction or quadratic term
iquad=FALSE
#gets first letter in term
t1<-substr(term,1,1)
if(t1=="I") {t2=term; iquad=TRUE}
#gets second letter in an interaction
if(!iquad) {t2<-substr(term,3,3) } else {t1<-substr(term,3,3)}
cmp<-FALSE
if(t2!= ""){cmp=TRUE}
if(cmp){t3<-substr(term,1,3)}
#gets vector of terms in the model
if(cmp) {if(iquad){mt<-c(t1,t2)} else {mt<-c(t1,t2,t3)}} else {mt<-t1}
# fits model model in factors with maximum correlation
d1<-data.frame(y=y,mm[,mt])
# if only one term in the model get the correct name
if(t2=="") {
nms<-names(d1)
nms[2]<-mt
names(d1)<-nms
            }
m1<-lm(y~(.),data=d1)
result<-summary(m1)
print(result)
return(mt)
}
