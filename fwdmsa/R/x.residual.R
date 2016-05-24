X.residual <- function(X,samp,member){
#computes residuals
#each observation remains in same restscore group for all items which is indicated by 'member'
N <- dim(X)[1]
J <- dim(X)[2]
resid.Xtot <- matrix(,N,J)
for(j in 1:J){
for(g in 1:max(member)){
 samp.member <- samp[which(member[samp]==g)]
 groupmember <- which(member==g)
 resid.Xtot[groupmember,j] <- X[groupmember,j]-mean(X[samp.member,j])
}
}
return(resid.Xtot)
}
