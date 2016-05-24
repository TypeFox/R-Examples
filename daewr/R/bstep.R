bstep <- function(y,des,prvm) {
dat<-cbind(y,des)
#get model from last step
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
#creates the data frame
dat<-data.frame(y=y,des)
#gets the model matrix of the full model
#gets the model matrix
lm1<-lm(y~(.)^2,data=dat)
mm<-model.matrix(lm1)
fact<-colnames(mm)
fact<-fact[-1]
fact<-paste(fact,collapse='+')
mod<-paste(c(fact,quad),collapse='+')
lm2<-lm(reformulate(termlabels=mod, response='y'),data=dat) # This works
mm<-model.matrix(lm2)
#deletes the constant column from the model matrix
mm<-mm[,2:ncol(mm)]
#creates data frame with terms from previous model
d1<-data.frame(y=y,mm[,prvm])
#fits the previous model
m3<-lm(y~(.),data=d1)
#mm<-model.matrix(m1)
d3<-d1
#mm<-model.matrix(y~(.)^2,data=dat)
d3<-data.frame(y=y,mm[,prvm]) ###This gives me an error on second pass : . problem???
m3<-lm(y~(.),data=d3)
#get coefficients
c3<-coefficients(m3)
#get standard error of coefficients
se3<-sqrt(sum((m3$residuals)^2)/(nrow(d3)-ncol(d3)))*sqrt(diag(solve(t(model.matrix(m3))%*%model.matrix(m3))))
#get t-values
t3<-c3/se3
#get F-values
f3<-t3^2
#get p-values
p3<-1-pf(f3,1,(nrow(d3)-ncol(d3))) 
#get index of term with higest p-value
p3<-p3[-1]
mxp<-max(p3)
idm<-which.max(p3)
#get name of term with highest p-value
#get vector of model terms
varname<-unlist(names(p3))
term<-varname[idm]
#remove term with higest p-value from data frame dat
trm<-varname[-idm]
#replace . with : in interaction terms
trm<-gsub("\\.",":",trm)
# fix quadratic terms
nm<-length(trm)
lquad=FALSE
i=1
 while (i<=nm) {
     if (substr(trm[i],1,1)=="I") {lquad=TRUE; nquad=i}
     if(lquad) {trm[nquad]<-gsub(":2:","^2)",trm[nquad])
     trm[nquad]<-gsub(":","(",trm[nquad])}
     i<-i+1
                   }
#fit the model with term with highest p-value removed
d3f<-d3[,-(idm+1)]
m3f<-lm(y~(.),data=d3f)
result<-summary(m3f)
print(result)
return(trm)
                  }