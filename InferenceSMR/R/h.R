h <-
function(t_i,delta_i,beta,time,z, S0_i, S1_i)
{
valeur  = matrix(0,ncol=1,nrow=length(z));
x = seq(1:length(t_i))[t_i<=time & delta_i==1];
#valeur = as.matrix(z)*sum(1/S0_i[x]) - apply(as.matrix(S1_i[,x]),1,sum);
valeur = as.matrix(z)*sum(1/S0_i[x]) - rowSums(as.matrix(S1_i[,x]));
valeur = c(exp(as.matrix(t(beta))%*%as.matrix(z)))*valeur/length(t_i);
}
