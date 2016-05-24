validation.CorrMat <-
function(prop.vec.bin, CorrMat) {


    if (is.positive.definite(CorrMat) == FALSE) {
           stop("Specified correlation matrix is not positive definite! \n")
        }
        if (isSymmetric(CorrMat) == FALSE) {
            stop("Specified correlation matrix is not symmetric! \n")
        }

if (length(prop.vec.bin )!=ncol(CorrMat)){
stop("Dimensions of binary probability vector and correlation matrix do not match")
}



d = length(prop.vec.bin )
sigma=CorrMat
p=prop.vec.bin ; q=1-p
no.bin=d

#Create lower and upper bound matrices for correlations

L_sigma=diag(d) ; U_sigma=diag(d)

#Find lower and upper bounds

for (i in 1:no.bin){
for (j in 1:no.bin){
if (i!=j) L_sigma[i,j]=L_sigma[j,i]=max(-sqrt((p[i]*p[j])/(q[i]*q[j])),
-sqrt((q[i]*q[j])/(p[i]*p[j])))
if (i!=j) U_sigma[i,j]=U_sigma[j,i]=min(sqrt((p[i]*q[j])/(q[i]*p[j])),
sqrt((q[i]*p[j])/(p[i]*q[j])))
}
}


valid.state=TRUE
for (i in 1:d){for (j in 1:d){ if(j >=i){
if(sigma[i,j] < L_sigma[i,j] | sigma[i,j] > U_sigma[i,j]){ 
stop("Range violation occurred in the binary data generation routine! \n 
  Corr[",i,",",j,"] must be between ", round(L_sigma[i,j],3),
 " and ",round(U_sigma[i,j],3),
 "\n This algorithm cannot generate ordinal data with specified correlations")
valid.state=FALSE}}}}

#if(valid.state==TRUE)  cat("All correlations are in feasible range! \n")
if(valid.state==FALSE) stop("All correlations must be in feasible range!")

}
