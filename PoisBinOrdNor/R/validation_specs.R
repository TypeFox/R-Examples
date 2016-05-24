validation_specs<-function(no.pois,no.bin,no.ord,no.norm,corr.mat=NULL,prop.vec.bin=NULL,prop.vec.ord=NULL,lamvec=NULL,nor.mean=NULL,nor.var=NULL){

nPois = length(lamvec)
n1=no.pois; n2=no.bin; n3=no.ord; n4=no.norm
d=n1+n2+n3+n4

if (no.pois!=nPois) {stop("Dimension of lamvec does not match the number of Poisson variables!\n") }
if (sum(lamvec<0)>0){stop("Lambda values cannnot be negative!\n") }

if (no.pois<0) { stop("Number of Poisson variables cannnot be negative!\n") }
if (no.bin<0|(floor(no.bin) != no.bin)) { stop("Number of binary variables must be a nonnegative integer!\n") }
if (no.ord<0) { stop("Number of ordinal variables cannnot be negative!\n") }
if (no.norm<0) { stop("Number of normal variables cannnot be negative!\n") }

if(!is.null(prop.vec.bin)) {
if (no.bin == 0) {stop("Proportion vector is specified while no.bin=0")
} else if((min(prop.vec.bin)<=0)|(max(prop.vec.bin)>=1)) {stop("Proportions for binary variables must be between 0 and 1!\n")
} else if(length(prop.vec.bin) != no.bin) {stop("Proportion vector is misspecified, dimension is wrong!\n")}
} else if (is.null(prop.vec.bin)) {if (no.bin > 0) {stop("Proportion vector is not specified while no.bin > 0")}}


if(no.ord !=0){
if(sum(is.na(prop.vec.ord))>0){stop("NAs are not permitted in prop.vec.ord")}
if(length(prop.vec.ord)!=no.ord) {warning("The value of no.ord and the number of probability vectors in the list do not match."); warn=TRUE} 

for (c in 1:length(prop.vec.ord)) {
pord=prop.vec.ord[[c]]
if ( max(pord)>=1) {warning("A probability value cannot be greater than 1. The probability vector ", i, " is not valid."); warn=TRUE }
if ( !(min(pord) >0) ) {warning("A probability value cannot be less than 0. The probability vector ", i, " is not valid."); warn=TRUE }
if ( sum( abs( pord[order(pord)]-pord) ) > 0 ) {warning("The probability vector ", i, " is not in the form of cumulative probabilities") ; warn=TRUE }
}}

is.wholenumber <-function(x, tol = .Machine$double.eps^0.5)  abs(x-round(x))<tol
if(!is.wholenumber(no.norm)){ stop("Number of normal variables cannnot be fractional number!\n") }
if((no.norm<0)|(floor(no.norm) != no.norm)) {stop("Number of normal variables \nmust be an integer whose value or 0!\n")}
if(!is.null(nor.mean) & no.norm==0) {stop("Mean vector for the normal part is specified while no.norm=0!\n")}
if(!is.null(nor.var) & no.norm==0) {stop("Vector of variances for the normal part is specified while no.norm=0!\n")}
if(is.null(nor.mean) & no.norm>0) {stop("Mean vector for the normal part is not specified while no.norm>0!\n")}
if(is.null(nor.var) & no.norm>0) {stop("Vector of variances for the normal part is not specified while no.norm>0!\n")}
    
if (!is.null(nor.mean)& !is.null(nor.var) & no.norm>0) {
if (length(nor.mean) != no.norm) {stop("Mean vector for the normal part is misspecified, \ndimension is wrong!\n")}
if (length(nor.var) != no.norm) {stop("Vector of variances for the normal part is misspecified, \ndimension is wrong!\n")}
if (min(nor.var) <= 0) {stop("Variances must be positive!\n")}
}


if(ncol(corr.mat)!=d){stop("Dimension of correlation matrix does not match the number of variables!\n")} 
    
if(is.positive.definite(corr.mat)==FALSE) {stop("Specified correlation matrix is not positive definite! \n")}

if(isSymmetric(corr.mat)==FALSE) {stop("Specified correlation matrix is not symmetric! \n")}
  
if(sum(corr.mat>1)>0) {stop("Correlation values cannot be greater than 1! \n") }
  
if(sum(corr.mat<(-1) )>0) {stop("Correlation values cannot be less than -1! \n")}
  
if(sum(diag(corr.mat)!=1)>0) {stop("All diagonal elements of the correlation matrix must be 1! \n") }


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


p=prop.vec.bin; q=1-p
sigma=corr.mat
L_sigma=diag(d); U_sigma=diag(d)

## Find lower and upper bounds for poisson-poisson combinations 
u = runif(100000, 0, 1)

if(no.pois>0){
for (i in 1:n1){
for (j in 1:n1){
if (i!=j) L_sigma[i,j]=L_sigma[j,i]=cor(qpois(u, lamvec[i]), qpois(1-u, lamvec[j]))
if (i!=j) U_sigma[i,j]=U_sigma[j,i]=cor(qpois(u, lamvec[i]), qpois(u, lamvec[j]))}}
}

## Find lower and upper bounds for poisson-binary combinations 


if(no.pois>0 & no.bin>0){
x=rbinom(100000, no.bin, prop.vec.bin)
z=rpois(100000, lamvec)
for(i in (n1+1):(n1+n2)){
for(j in 1:n1){
if (i!=j) L_sigma[i,j]=L_sigma[j,i]=cor(x[order(x, decreasing=TRUE)], z[order(z)])
if (i!=j) U_sigma[i,j]=U_sigma[j,i]=cor(x[order(x)], z[order(z)])}}
}

## Find lower and upper bounds for poisson-ordinal combinations 

if(no.pois>0 & no.ord>0){
for(i in (n1+n2+1):(n1+n2+n3)){
y=rnorm(100000, 0,1)
pvec=prop.vec.ord[[i-n1-n2]]
yord = numeric(length(y))
for (r in 1:length(pvec)){
if (r !=length(pvec)) {
t1 = qnorm(pvec[r])
t2 = qnorm(pvec[r+1] )
yord[(t1<y)&(y<=t2)]= r
} else {
yord[y>qnorm(pvec[r])]= r
}}
yord=yord+1
z=rpois(100000, lamvec)
for(j in 1:n1){
if (i!=j) L_sigma[i,j]=L_sigma[j,i]=cor(yord[order(yord, decreasing=TRUE)], z[order(z)])
if (i!=j) U_sigma[i,j]=U_sigma[j,i]=cor(yord[order(yord)], z[order(z)])}}
}


## Find lower and upper bounds for poisson-normal combinations 

if(no.pois>0 & no.norm>0){  
for(i in (n1+n2+n3+1):d){
for(j in 1:n1){
y=rnorm(100000, 0,1)
z=rpois(100000, lamvec)
if (i!=j) L_sigma[i,j]=L_sigma[j,i]=cor(z[order(z,decreasing=TRUE)],y[order(y)])
if (i!=j) U_sigma[i,j]=U_sigma[j,i]=cor(z[order(z)],y[order(y)])}}
}

## Find lower and upper bounds for binary-binary combinations 

if(no.bin>0){
for (i in (n1+1):(n1+n2)){
for (j in (n1+1):(n1+n2)){
if (i!=j) L_sigma[i,j]=L_sigma[j,i]=max(-sqrt((p[i-n1]*p[j-n1])/(q[i-n1]*q[j-n1])), -sqrt((q[i-n1]*q[j-n1])/(p[i-n1]*p[j-n1])))
if (i!=j) U_sigma[i,j]=U_sigma[j,i]=min(sqrt((p[i-n1]*q[j-n1])/(q[i-n1]*p[j-n1])), sqrt((q[i-n1]*p[j-n1])/(p[i-n1]*q[j-n1])))}}
}

## Find lower and upper bounds for binary-ordinal combinations 

if(no.bin>0&no.ord>0){
for(i in (n1+n2+1):(n1+n2+n3)){
for(j in (n1+1):(n1+n2)){
x=rbinom(100000, no.bin, prop.vec.bin)
if (i!=j) L_sigma[i,j]=L_sigma[j,i]=cor(yord[order(yord, decreasing=TRUE)], x[order(x)])
if (i!=j) U_sigma[i,j]=U_sigma[j,i]=cor(yord[order(yord)], x[order(x)])}}
}

## Find lower and upper bounds for binary-normal combinations 

if(no.bin>0&no.norm>0){
for(i in (n1+n2+n3+1):d){
for(j in (n1+1):(n1+n2)){
if (i!=j) L_sigma[i,j]=L_sigma[j,i]= -dnorm(qnorm(p[j-n1]))/sqrt(p[j-n1]*q[j-n1])
if (i!=j) U_sigma[i,j]=U_sigma[j,i]= dnorm(qnorm(p[j-n1]))/sqrt(p[j-n1]*q[j-n1])}}
}

## Find lower and upper bounds for ordinal-ordinal combinations

if(no.ord>0){
for(i in (n1+n2+1):(n1+n2+n3)){
for(j in (n1+n2+1):(n1+n2+n3)){

w=rnorm(100000, 0,1)
pvec=prop.vec.ord[[j-n1-n2]]
word = numeric(length(w))
for (r in 1:length(pvec)){
if (r !=length(pvec)) {
t1 = qnorm(pvec[r])
t2 = qnorm(pvec[r+1] )
word[(t1<w)&(w<=t2)]= r
} else {
word[w>qnorm(pvec[r])]= r
}}
word=word+1
if (i!=j) L_sigma[i,j]=L_sigma[j,i]=cor(yord[order(yord,decreasing=TRUE)],word[order(word)]) 
if (i!=j) U_sigma[i,j]=U_sigma[j,i]=cor(yord[order(yord)],word[order(word)])}}
}

## Find lower and upper bounds for ordinal-normal combinations 

if(no.ord>0 & no.norm>0){
for(j in (n1+n2+1):(n1+n2+n3)){ #ordinal
for(i in (n1+n2+n3+1):d){ #normal
z=rnorm(100000, 0,1)
if (i!=j) L_sigma[i,j]=L_sigma[j,i]=cor(word[order(word,decreasing=TRUE)],z[order(z)]) 
if (i!=j) U_sigma[i,j]=U_sigma[j,i]=cor(word[order(word)],z[order(z)])}}
}


# Find lower & upper bounds for normal-normal combinations

if(no.norm>0){
for (i in (n1+n2+n3+1):d){
for (j in (n1+n2+n3+1):d){
if (i!=j) L_sigma[i,j]=L_sigma[j,i]=-1
if (i!=j) U_sigma[i,j]=U_sigma[j,i]= 1}}
}

valid.state=TRUE
for (i in 1:d){for (j in 1:d){ if(j >=i){
 if(sigma[i,j] < L_sigma[i,j] | sigma[i,j] > U_sigma[i,j]){ 
cat("Range violation! Corr[",i,",",j,"] must be between", round(L_sigma[i,j],3),"and",round(U_sigma[i,j],3),"\n")
valid.state=FALSE }
}}}

if(valid.state==TRUE)  cat("All correlations are in feasible range! \n")
if(valid.state==FALSE) stop("All correlations must be in feasible range!")
return(TRUE)

}
