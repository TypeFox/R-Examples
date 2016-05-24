function1 <-
function(x)
{
if(class(x)!="coxph"){return("This function requires an object of class coxph")} #Makes sure the object is of the good class

omega_1 = x$n*x$var;
beta = as.matrix(x$coefficients);

#Creating objects that will be used in the variance calcul

S0<-function(t_i,covar_i,beta,t)
{
1/length(t_i)*exp(as.matrix(t(beta))%*%as.matrix(covar_i))%*%as.matrix(t_i>=t)
}

S1<-function(t_i,covar_i,beta,t)
{
1/length(t_i)*covar_i%*%(as.matrix(t_i>=t)*t(exp(t(beta)%*%as.matrix(covar_i))))
}

n<-x$n;
ti<-coxph.detail(x)$y[,2];
S0_i<-numeric(n);
covar = t(coxph.detail(x)$x);
deceased=coxph.detail(x)$y[,3];
S1_i<-matrix(0,ncol=n,nrow=length(covar[,1]));
for(i in 1:n)
{
S0_i[i]<-S0(ti,covar,beta,ti[i]);
S1_i[,i]<-S1(ti,covar,beta,ti[i])/S0_i[i]**2;
}
S02_i=S0_i**2;
#S0 and S1 are not really the objects found in Lin and Flemming (1994)
#The objects were modified for computationnal efficiency

return(list(ti,deceased,covar,omega_1, beta, S0_i, S1_i, S02_i));
}
