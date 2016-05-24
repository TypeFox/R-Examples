covariance <-
function(t_i,delta_i,beta,time1,z1,time2,z2, S0_i, S1_i, S02_i, Omega_1)
{
z1 = t(z1);
z2 = t(z2);
n=length(t_i);
cov = 0;
x = seq(1:n)[t_i<=min(time1,time2) & delta_i!=0];
cov = sum(1/S02_i[x]);
cov = 1/n*exp(t(beta)%*%(z1+z2))*cov + as.matrix(t(h(t_i,delta_i,beta,time1,z1, S0_i, S1_i)))%*%Omega_1%*%as.matrix(h(t_i,delta_i,beta,time2,z2, S0_i, S1_i));
return(cov/n);
}
