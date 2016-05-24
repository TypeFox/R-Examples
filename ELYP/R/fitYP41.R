fitYP41 <- function(Y, d, Z, beta1=1, beta2= -1, maxiter=60){
#### Do two things: (1) for the data, and given beta1, beta2;
#### compute the baseline that max the EL. i.e. NPMLE.
#### (2) Given the baseline and the 2 betas, compute the (log) EL value.
####  there is no alpha, So Z is a vector or nx1 matrix.

temp1 <- YP41(y=Y, d=d, Z=Z, b1=beta1, b2=beta2, k=maxiter)

ELval <- ELcomp(Haz=temp1$Hazw, Sur=temp1$Survival, gam=temp1$gam)

list(EmpLik=ELval, BaselineH=temp1$Hazw)
}

