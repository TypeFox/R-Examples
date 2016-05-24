BinToOrd <-
function(prop.vec.bin, ordPmat, Mlocation, bin.data){
myord = validation.ordPmat(ordPmat)
J = myord$J
K = myord$K
ord.data=matrix(NA,nrow(bin.data), ncol(bin.data) )


for (j in 1:J) {
p=numeric(0)
for (k in 1:K[j]) {
if (k<Mlocation[j]) { p[k] = ordPmat[k,j]/(1-prop.vec.bin[j]) }
else { p[k] = ordPmat[k,j]/(1-prop.vec.bin[j]) }
}
w1 = bin.data[,j]==0
if ((Mlocation[j]-1)==1){ord.data[w1,j] = 1}
else{
ord.data[w1,j] = sample(1:(Mlocation[j]-1), sum(w1), prob=p[1:(Mlocation[j]-1)], replace=TRUE)
}
if (Mlocation[j]==K[j]){ord.data[!w1,j]=K[j]}
else{
ord.data[!w1,j] = sample(Mlocation[j]:K[j], sum(!w1), prob=p[Mlocation[j]:K[j]], replace=TRUE)
}
}
return(list(y=ord.data, Corr=cor(ord.data)))
}
