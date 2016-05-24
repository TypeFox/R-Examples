var.expDeath <-
function(contribution,incid,cox,fuzz = 0.01, Poisson = FALSE, covnames)
{
f1 = function1(cox); 
times_np = f1[[1]]; deceased = f1[[2]];
covar = f1[[3]]; Omega_1 = f1[[4]];beta = f1[[5]];
S0_i = f1[[6]]; S1_i = f1[[7]]; S02_i = f1[[8]]; 
increment = as.numeric(contribution[[4]]);
p_incid = as.numeric(contribution[[2]]);
p_surv = as.numeric(contribution[[3]]);
ti = as.matrix(contribution[[1]][1]);
cov_incid = as.matrix(contribution[[1]][2:(1+p_incid)]);
cov_surv = as.data.frame(contribution[[1]][(2+p_incid):(1+p_incid+p_surv)]);
colnames(cov_surv) = covnames;
freq = as.matrix(contribution[[1]][,2+p_incid+p_surv]);


equal_row<-function(x,y,fuzz = 0)
{
(sum(abs(x - y)) < fuzz);
}

surv_p<-function(x,cov,time)
{
y = survfit(x,newdata=cov, se.fit = F);
return(approx(x=y$time,y=y$surv,xout=time,rule=2)$y);
}

variance1 = 0;
variance2 = 0;

if(Poisson == FALSE)
{
for(i in seq(1,length(freq))[ti!=0])
{
g = seq(1,length(incid[,1]),1)[apply(incid[,2:(1+p_incid)],1,function(x){equal_row(x,cov_incid[i,],fuzz)})];
#g = seq(1,length(incid[,1]),1)[apply(incid[,2:(1+p_incid)], 1, isTRUE(all.equal), cov_incid[i,])];
if(incid[g,1]==0) next;

for(j in seq(i,length(freq))[ti[seq(i,length(freq))]!=0])
{
g_prime = seq(1,length(incid[,1]),1)[apply(incid[,2:(1+p_incid)],1,function(x){equal_row(x,cov_incid[j,],fuzz)})];
#g_prime = seq(1,length(incid[,1]),1)[apply(incid[,2:(1+p_incid)], 1, isTRUE(all.equal), cov_incid[j,])];
if(incid[g_prime,1]==0) next                                        
facteur = 2*freq[i]*freq[j]; #freq individus sont dans cette catégories
if(i==j) facteur = 1/2*facteur;
if(g==g_prime)
{
variance1 = variance1 + facteur*((1 - surv_p(cox,cov_surv[i,],ti[i]))*
(1 - surv_p(cox,cov_surv[j,],ti[j]))*
incid[g,1]*(1-incid[g,1])/incid[g,p_incid+2] +
(incid[g,1]*(1-incid[g,1])/incid[g,p_incid+2] + incid[g,1]**2)*
surv_p(cox,cov_surv[i,, drop = F],ti[i])*
surv_p(cox,cov_surv[j,, drop = F],ti[j])*
covariance(times_np, deceased, beta, ti[i], cov_surv[i,], ti[j], cov_surv[j,], S0_i, S1_i, S02_i, Omega_1));
}
else #cas g != g'
{
variance2 = variance2 + facteur*incid[g,1]*incid[g_prime,1]*
surv_p(cox,cov_surv[i,, drop = F],ti[i])*
surv_p(cox,cov_surv[j,, drop = F],ti[j])*
covariance(times_np, deceased, beta, ti[i],cov_surv[i,],ti[j],cov_surv[j,], S0_i, S1_i, S02_i, Omega_1);
}
}
}
}

if(Poisson == TRUE)
{
for(i in seq(1,length(freq))[ti!=0])
{
g = seq(1,length(incid[,1]),1)[apply(incid[,2:(1+p_incid)],1,function(x){equal_row(x,cov_incid[i,],fuzz)})];
#g = seq(1,length(incid[,1]),1)[apply(incid[,2:(1+p_incid)], 1, isTRUE(all.equal), cov_incid[i,])];

if(incid[g,1]==0) next;

for(j in seq(i,length(freq))[ti[seq(i,length(freq))]!=0])
{
g_prime = seq(1,length(incid[,1]),1)[apply(incid[,2:(1+p_incid)],1,function(x){equal_row(x,cov_incid[j,],fuzz)})];
#g_prime = seq(1,length(incid[,1]),1)[apply(incid[,2:(1+p_incid)], 1, isTRUE(all.equal), cov_incid[j,])];
if(incid[g_prime,1]==0) next                                        
facteur = 2*freq[i]*freq[j]; #freq individus sont dans cette catégories
if(i==j) facteur = 1/2*facteur;
if(g==g_prime)
{
variance1 = variance1 + facteur*((1 - surv_p(cox,cov_surv[i,],ti[i]))*
(1 - surv_p(cox,cov_surv[j,],ti[j]))*
incid[g,1]/incid[g,p_incid+2] +
(incid[g,1]/incid[g,p_incid+2] + incid[g,1]**2)*
surv_p(cox,cov_surv[i,, drop = F],ti[i])*
surv_p(cox,cov_surv[j,, drop = F],ti[j])*
covariance(times_np, deceased, beta, ti[i],cov_surv[i,],ti[j],cov_surv[j,], S0_i, S1_i, S02_i, Omega_1));
}
else #cas g != g'
{
variance2 = variance2 + facteur*incid[g,1]*incid[g_prime,1]*
surv_p(cox,cov_surv[i,, drop = F],ti[i])*
surv_p(cox,cov_surv[j,, drop = F],ti[j])*
covariance(times_np, deceased, beta, ti[i],cov_surv[i,],ti[j],cov_surv[j,], S0_i, S1_i, S02_i, Omega_1);
}
}
}
}
return(increment**2*(variance1+variance2));
}
