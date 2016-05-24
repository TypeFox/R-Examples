
Selvar <- function(Y,Z,X,thresh_vect,nb_boot=80,nb_repli=50,CI_level=0.95,nb_cores=1)
{
#### Z has to be scaled before using this function
  
N=dim(Z)[2]
n=length(Y)

###Parameter
alpha=1
level=1-CI_level

###Initialisation
corre=matrix(0,N)
etachap_fin=c()
CI_low=c()
CI_up=c()
selec_ind=list()

### Projection on the orthogonal of the image of X
tXX=t(X)%*%X
inv_g=ginv(tXX)
Pg=diag(n)-X%*%inv_g%*%t(X)
eige=eigen(Pg)
vp=eige$values
dim_imX_orth=round(sum(vp))
Ag=eige$vectors
Ag_t=Ag[,1:dim_imX_orth]
cat("Dimension of the orthogonal of the image of X =",dim_imX_orth,"\n")
cat("It has to be equal to n=",length(Y),"-the number of null eigenvalues of the projection\n")
#### To check
tAgAg=t(Ag_t)%*%Ag_t
cat("This sum :",sum(tAgAg-diag(dim_imX_orth)),"has to be small\n")

#### Projection
Z_proj=t(Ag_t)%*%Z
Y_proj=t(Ag_t)%*%Y

n=length(Y_proj)
##correlation
for (j in 1:N)
{
  corre[j]=abs(sum(Y_proj*Z_proj[,j]))
}
corre_sort=sort(corre,decreasing=TRUE,index.return=TRUE)

Top_indices=(corre_sort$ix)[1:n]

##### Estimation de eta Top indices #####
Z_Top=Z_proj[,Top_indices]

stabsel.glmnet <- function(i) 
{
  cat("+")
  b_sort <- sort(sample(1:n,floor(n/2)))
  resultat_glmnet=glmnet(Z_Top[b_sort,],Y_proj[b_sort],family="gaussian",alpha=alpha)
  nb_colonnes=dim(resultat_glmnet$beta)[2]
  ind_glmnet=which(resultat_glmnet$beta[,nb_colonnes]!=0)
  return(tabulate(ind_glmnet,n))
}
res.cum <- Reduce("+", mclapply(1:nb_repli, stabsel.glmnet, mc.cores=nb_cores))
cat('\n')

for (j in 1:length(thresh_vect))
{
  thresh=thresh_vect[j]
  ind_glmnet=which(res.cum/nb_repli>thresh)
  res_HiLMM=estim_herit(Y_proj,Z_Top[,ind_glmnet])
  etachap_fin[j]=res_HiLMM$heritability
  sigma2_hat=res_HiLMM$sig2
  cat("threshold=",thresh,", heritability=",etachap_fin[j],"\n",sep='')
  
  ###IC bootstrap ###

  res_boot_corr=bootstrap_corr(Y_proj,Z_Top[,ind_glmnet],nb_boot,etachap_fin[j],sigma2_hat,level,nb_cores)
  CI_low[j]=res_boot_corr$CI_low
  CI_up[j]=res_boot_corr$CI_up
  selec_ind[[j]]=Top_indices[ind_glmnet]
  cat("CI_low=",CI_low[j],", CI_up=",CI_up[j],"\n",sep='')
}

list(heritability=etachap_fin,CI_low=CI_low,CI_up=CI_up,selec_ind=selec_ind)
}
