curve(pt(x,25,ncp=3), from=0, to=6)
abline(v=qt(.975,25))
if (.make.epsf) dev.copy2eps(file="noncentral-t.ps")
pt(qt(.975,25),25,ncp=3)
power.t.test(delta=0.5, sd=2, sig.level = 0.01, power=0.9)
power.t.test(n=450, delta=0.5, sd=2, sig.level = 0.01)
power.t.test(delta=0.5, sd=2, sig.level = 0.01, power=0.9,
alt="one.sided")
power.t.test(delta=10, sd=10*sqrt(2), power=0.85, type="paired")
power.prop.test(power=.85,p1=.15,p2=.30)
### no data sets in exercises
rm(list=ls())
while(search()[2] != "package:ISwR") detach()
