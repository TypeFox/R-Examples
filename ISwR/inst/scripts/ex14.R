library(survival)    
attach(graft.vs.host)
plot(survfit(Surv(time,dead)~gvhd))
survdiff(Surv(time,dead)~gvhd)
summary(coxph(Surv(time,dead) ~ gvhd)) # for comparison
summary(coxph(Surv(time,dead) ~ 
              gvhd + log(index) + donage + rcpage + preg))
attach(melanom)
cox1 <- coxph(Surv(days, status==1) ~ 
              log(thick) + sex + strata(ulc))
new <- data.frame(sex=2, thick=c(0.1, 0.2, 0.5))
svfit <-  survfit(cox1,newdata=new)
plot(svfit[2], ylim=c(.985, 1))
summary(coxph(Surv(obsmonths, dead)~age+sex, data=stroke))
summary(coxph(Surv(obsmonths, dead)~sex, data=stroke))
with(stroke, tapply(age,sex,mean))
stroke.trim <- function(t1, t2) 
   subset(transform(stroke,
                    entry=t1, exit=pmin(t2, obsmonths), 
                    dead=dead & obsmonths <= t2), 
          entry < exit)
stroke2 <- do.call(rbind, mapply(stroke.trim, 
       c(0,0.5,2,12), c(0.5,2,12,Inf), SIMPLIFY=F))
summary(coxph(Surv(entry, exit, dead)~age+sex, data=stroke2))
rm(list=ls())
while(search()[2] != "package:ISwR") detach()
