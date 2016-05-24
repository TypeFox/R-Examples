adaptive.case.test <-
function(no.exposed.in.matched.set,narrowcase,case.exposed,size.matched.set,Gamma,narrowmultiplier=1,alpha=.05){
if(narrowmultiplier!=round(narrowmultiplier)){
stop("narrowmultiplier must be an integer");
}
# Only work with discordant matched sets
discordant=(no.exposed.in.matched.set>0)*(no.exposed.in.matched.set<size.matched.set);
critmat=jointcrit(mhxct.weighted(no.exposed.in.matched.set[((discordant==1)*(narrowcase==1))==1],size.matched.set[((discordant==1)*(narrowcase==1))==1],gamma=Gamma,m=narrowmultiplier),mhxct.weighted(no.exposed.in.matched.set[(discordant==1)*(narrowcase==0)==1],size.matched.set[(discordant==1)*(narrowcase==0)==1],gamma=Gamma,m=1));


critical.values=adaptive.test.critical.value.func(critmat,alpha=alpha);
t1.critical.value=critical.values$t1;
t1plust2.critical.value=critical.values$t1plust2;
t1=narrowmultiplier*sum((case.exposed[narrowcase==1]==1)*(discordant[narrowcase==1]==1));
t1plust2=t1+sum((case.exposed[narrowcase==0]==1)*discordant[narrowcase==0]==1);
reject=(t1>=t1.critical.value)|(t1plust2>=t1plust2.critical.value);
if(reject==TRUE){
overall.test.result="Reject H0";
}
if(reject==FALSE){
overall.test.result="Accept H0";
}

# Modify t1 so that it is the same as in the paper
modt1=t1/narrowmultiplier;
modt1.critical.value=t1.critical.value/narrowmultiplier;
narrowtestdef="T1";
broadtestdef=paste(narrowmultiplier,"*T1+T2",sep="");
narrowtest.critval=ceiling(modt1.critical.value);
broadtest.critval=t1plust2.critical.value;
narrowtest.obsval=modt1;
broadtest.obsval=t1plust2;
testdef=c(narrowtestdef,broadtestdef);
critvals=c(narrowtest.critval,broadtest.critval);
obsvals=c(narrowtest.obsval,broadtest.obsval);
testinfomat=data.frame(cbind(testdef,critvals,obsvals),row.names=c("Narrow Test Statistic","Narrow + Marginal Test Statistic"));
names(testinfomat)=c("Definition","Critical Value","Observed Value")

list(alpha=alpha,Gamma=Gamma,testinfomat=testinfomat,overall.test.result=overall.test.result)
}

