`path.analysis` <-
function(corr.x,corr.y) {
if (ncol(corr.y)>1)corr.y<-t(corr.y)
Direct<-solve(corr.x,corr.y)
n<-ncol(corr.x)
Coeff <- corr.x
for ( i in 1:n) {
for ( j in 1:n) {
Coeff[i,j]<-Direct[j]*corr.x[i,j]
}
}
Residual<-1-t(Direct)%*%corr.y
cat("Direct(Diagonal) and indirect effect path coefficients",
    "\n======================================================\n")
print(Coeff)
cat("\nResidual Effect^2 = ",Residual,"\n")
output<-list(Coeff=Coeff,Residual=as.numeric(Residual))
invisible(output)
}

