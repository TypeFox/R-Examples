dependCox.reg.CV <-
function(t.vec,d.vec,X.mat,K=5,G=20){

### Grid search on CV ####
tau_grid=seq(0.0001,0.9,length=G)
alpha_grid=2*tau_grid/(1-tau_grid)
C_grid=NULL
for(i in 1:length(alpha_grid)){
  C_grid[i]=cindex.CV(t.vec,d.vec,X.mat,alpha_grid[i])
}
alpha_hat=alpha_grid[C_grid==max(C_grid)][1]

######### univariate Cox with dependent censoring #########
p=ncol(X.mat)
P_d=numeric(p)
beta_d=numeric(p)
SD_d=numeric(p)
Z_d=numeric(p)

for(j in 1:p){
  res=dependCox.reg(t.vec,d.vec,X.mat[,j],alpha=alpha_hat,var=TRUE)
  beta_d[j]=res[1]
  SD_d[j]=res[2]
  Z_d[j]=res[1]/res[2]
  P_d[j]=1-pchisq(Z_d[j]^2,df=1)
}

plot(tau_grid,C_grid,xlab="Kendall's tau ( alpha )",ylab="CV( alpha )",type="b",lwd=3)
points(tau_grid[C_grid==max(C_grid)][1],max(C_grid),col="red",pch=17,cex=2)

list(beta_hat=beta_d,SD=SD_d,Z=Z_d,P=P_d,alpha=alpha_hat)

}
