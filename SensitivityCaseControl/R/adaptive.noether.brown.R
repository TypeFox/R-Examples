adaptive.noether.brown <-
function(y,Gamma,alpha=.05,lambda1=1/3,lambda2=2/3){
nobrown=function(y,fr1=1/3,fr2=2/3){
        yo<-rev(y[order(abs(y))])
        n<-length(y)
        rk<-rank(abs(y))
        top<-rk>=((1-fr1)*n)
        high<-rk>=((1-fr2)*n)
        I1<-sum(top)
        I2<-sum(high-top)
        B1<-sum((y>0)[top])
        T1<-sum((y>0)[high])+B1
        list(I1I2=c(I1,I2),B1T1=c(B1,T1))
}

nobrownstats=nobrown(y,fr1=lambda1,fr2=lambda2)
critmat=jointcrit(mhxct.weighted(rep(1,nobrownstats$I1I2[1]),rep(2,nobrownstats$I1I2[1]),gamma=Gamma,m=2),mhxct.weighted(rep(1,nobrownstats$I1I2[2]),rep(2,nobrownstats$I1I2[2]),gamma=Gamma,m=1));
t1=nobrownstats$B1T1[1]*2; # This is Noether's test statistic multiplied by 2
t1plust2=nobrownstats$B1T1[2]; # This is Brown's test statistic
critical.values=adaptive.test.critical.value.func(critmat,alpha=alpha);
t1.critical.value=critical.values$t1;
t1plust2.critical.value=critical.values$t1plust2;
reject=(t1>=t1.critical.value)|(t1plust2>=t1plust2.critical.value);
if(reject==TRUE){
overall.test.result="Reject H0";
}
if(reject==FALSE){
overall.test.result="Accept H0";
}
critval.noether=ceiling(t1.critical.value/2);
critval.brown=t1plust2.critical.value;
critvals=c(critval.noether,critval.brown);
obsval.noether=nobrownstats$B1T1[1]
obsval.brown=nobrownstats$B1T1[2];
obsvals=c(obsval.noether,obsval.brown);
testinfomat=data.frame(cbind(critvals,obsvals),row.names=c("Noether Test Statistic","Brown Test Statistic"));
names(testinfomat)=c("Critical Value","Observed Value");

list(alpha=alpha,Gamma=Gamma,testinfomat=testinfomat,overall.test.result=overall.test.result);
}

