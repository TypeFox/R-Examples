qqnorm(react)
t.test(react)
wilcox.test(react)
wilcox.test(vital.capacity~group, data=vitcap)
attach(intake) ; opar <- par(mfrow=c(2,2))
plot(post ~ pre) ; abline(0,1)
plot((post+pre)/2, post - pre, 
    ylim=range(0,post-pre)); abline(h=0)
hist(post-pre)
qqnorm(post-pre)
detach(intake)
par(opar)
shapiro.test(react)
shapiro.test(react[-c(1,334)])
qqnorm(react[-c(1,334)])
attach(ashina)
t.test(vas.active, vas.plac, paired=TRUE)
t.test((vas.active-vas.plac)[grp==1], 
       (vas.plac-vas.active)[grp==2])
t.test(rnorm(25))$p.value       #repeat 10x
t.test(rt(25,df=2))$p.value     #repeat 10x
t.test(rexp(25), mu=1)$p.value  #repeat 10x
x <- replicate(5000, t.test(rexp(25), mu=1)$p.value)
qqplot(sort(x),ppoints(5000),type="l",log="xy")
rm(list=ls())
while(search()[2] != "package:ISwR") detach()
