est.expDeath <-
function(contribution,incid,cox,fuzz = 0.01,covnames)
{
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
(sum(abs(x - y)) <= fuzz);
}
surv_p<-function(x,cov,time)
{
y = survfit(x, newdata=cov, se.fit = F);
return(approx(x=y$time,y=y$surv,xout=time,rule=2)$y);
}

attendu = 0;
for(i in seq(1,length(freq))[ti!=0])
{
g = seq(1,length(incid[,1]),1)[apply(incid[,2:(1+p_incid)],1,function(x){equal_row(x,cov_incid[i,],fuzz)})];
if(incid[g,1]==0) next;
attendu = attendu + increment*freq[i]*incid[g,1]*(1 - surv_p(cox, cov_surv[i,, drop = F], ti[i]));
}
return(attendu);
}
