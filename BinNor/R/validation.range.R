validation.range <-
function(no.bin, no.nor, prop.vec.bin=NULL, corr.mat){

d = no.bin+no.nor
sigma=corr.mat
p=prop.vec.bin ; q=1-p

#Create lower and upper bound matrices for correlations

L_sigma=diag(d) ; U_sigma=diag(d)

## Find lower and upper bounds for binary-binary combinations 
#(STEP 2 in the algorithm)
if(no.bin>0){
for (i in 1:no.bin){
for (j in 1:no.bin){
if (i!=j) L_sigma[i,j]=L_sigma[j,i]=max(-sqrt((p[i]*p[j])/(q[i]*q[j])),
-sqrt((q[i]*q[j])/(p[i]*p[j])))
if (i!=j) U_sigma[i,j]=U_sigma[j,i]=min(sqrt((p[i]*q[j])/(q[i]*p[j])),
sqrt((q[i]*p[j])/(p[i]*q[j])))}}

## Find lower and upper bounds for binary-normal combinations 
#(STEP 3 in the algorithm)
}

if(no.bin>0&no.nor>0){
for(i in (no.bin+1):d){
for(j in 1:no.bin){
L_sigma[i,j]=L_sigma[j,i]= -dnorm(qnorm(p[j]))/sqrt(p[j]*q[j])
U_sigma[i,j]=U_sigma[j,i]= dnorm(qnorm(p[j]))/sqrt(p[j]*q[j])}}
}
#Lower & upper bounds for normal-normal combinations

if(no.nor>0){
for (i in (no.bin+1):d){
     for (j in (no.bin+1):d){
     L_sigma[i,j]=L_sigma[j,i]=-1
     U_sigma[i,j]=U_sigma[j,i]= 1}}
}
## Check if the correlations are in the feasible range 
#(STEP 4 in the algorithm)

valid.state=TRUE
for (i in 1:d){for (j in 1:d){ if(j >=i){
 if(sigma[i,j] < L_sigma[i,j] | sigma[i,j] > U_sigma[i,j]){ 
cat("Range violation! Corr[",i,",",j,"] must be between", round(L_sigma[i,j],3),
"and",round(U_sigma[i,j],3),"\n")
valid.state=FALSE}}}}

if(valid.state==TRUE)  cat("All correlations are in feasible range! \n")
if(valid.state==FALSE) stop("All correlations must be in feasible range!")
}

