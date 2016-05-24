genPBONdata<-
function(n, no_pois,no_bin,no_ord,no_norm,inter.mat=NULL,lamvec=NULL,prop_vec_bin=NULL,prop_vec_ord=NULL,nor.mean=NULL,nor.var=NULL){

n1=no_pois; n2=no_bin; n3=no_ord; n4=no_norm

d=n1+n2+n3+n4

xx1=rmvnorm(n, rep(0, d), inter.mat)

if (n1!=0) {
PP=matrix(0, n, n1)
for (i in 1: length(lamvec)) {
PP[,i]=qpois(pnorm(xx1[,i]), lamvec[i])}
} else PP=NULL


if (n2!=0){
BB=matrix(0,n, n2) 
for (j in (n1+1):(n1+n2)){
for (i in 1:n){
if (1*xx1[i,j]>qnorm(1-prop_vec_bin[j-n1])) BB[i,j-n1]=1 else BB[i,j-n1]=0}
}} else BB=NULL


if (n3!=0){
OO=matrix(0, n, n3)
for (i in 1:length(prop_vec_ord)){
y=xx1[,(n1+n2+i)]
pvec=prop_vec_ord[[i]]
yord = numeric(n)
for (r in 1:length(pvec)){
if (r !=length(pvec)) {
t1 = qnorm(pvec[r])
t2 = qnorm(pvec[r+1] )
yord[(t1<y)&(y<=t2)]= r
} else {
yord[y>qnorm(pvec[r])]= r
}
}
OO[,i]=yord}
} else OO=NULL

if (n4!=0){
NN=t(t(xx1[, (n1+n2+n3+1):d])*sqrt(nor.var)+nor.mean)
} else NN=NULL

data=cbind(PP, BB, OO, NN)

final.corr=cor(data)

result<-list( 
n.rows=n, prob.bin=prop_vec_bin, prob.ord=prop_vec_ord, nor.mean=nor.mean, nor.var=nor.var, lamvec=lamvec,
no.pois=n1, no.bin=n2, no.ord=n3, no.norm=n4, data=data)

return(result)
}

