contrib <-
function(start_follow, end_follow, incid_cov, surv_cov, follow_up, increment)
{
start_follow2 = increment*ceiling(start_follow/increment); #Rounding up at increment the follow-up times
end_follow2 = increment*floor(end_follow/increment);
end_follow2[end_follow == follow_up] = increment*round(end_follow[end_follow == follow_up]/increment);
follow_up2 = increment*round(follow_up/increment);
x = seq(1:length(start_follow2))[start_follow2 <= end_follow2]; 
start_follow2 = start_follow2[x]; end_follow2 = end_follow2[x];
incid_cov = incid_cov[x,]; surv_cov = surv_cov[x,];
follow_up2 = follow_up2[x];
cov=cbind(incid_cov,surv_cov);

p=length(c(surv_cov[1,],incid_cov[1,]));

ti2<-matrix(unlist(apply(cbind(start_follow2,end_follow2,follow_up2),1,function(x){seq(x[3]-x[1],x[3] - x[2],-increment)})),ncol=1,byrow=F);
#We want the follow up time left

temp<-matrix(0,nrow=length(ti2),ncol=p);
for(i in 1:p)
{
temp[,i]=unlist(apply(cbind(start_follow2,end_follow2,follow_up2,cov[,i]),1,function(x){rep(x[4],length(seq(x[3] - x[1],x[3] - x[2],-increment)))}));
}

w1<-matrix(unlist(apply(cbind(start_follow2,end_follow2,follow_up2),1,function(x){rep(x[3],length(seq(x[3] - x[1],x[3] - x[2],-increment)))})),ncol=1,byrow=F);

cont = cbind(ti2,w1,temp);
cont[,2] = apply(cont,1,function(cont)
{
if(cont[1] < increment/2)#We want to test if start_follow2 == 0, but we don't trust == function
{
cont[2] = 0.5;
}
else if(cont[2] - cont[1] < increment/2)#We want to test if start_follow2 == follow_up2
{
cont[2] = 0.5;
}
else cont[2] = 1;
});

by_cov2 = vector(mode='list',length=p+1) 
by_cov2[[1]] = cont[,1];
for(i in 1:p){by_cov2[[i+1]] = cont[,i+2]}

return(list(aggregate(cont[,2],by_cov2,sum),length(incid_cov[1,]),length(surv_cov[1,]),increment))
}
