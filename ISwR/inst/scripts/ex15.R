bcmort2 <- within(bcmort,{
  period <- area <- cohort
  levels(period) <- rep(c("1991-2001","1981-1991"), each=2)
  levels(area) <- rep(c("Cph+Frb","Nat"),2) 
})
bcfit <- glm(bc.deaths ~ (age + period + area)^2, poisson,
             offset=log(p.yr), data=bcmort2)
summary(bcfit)
drop1(bcfit, test="Chisq")
confint(bcfit, parm="period1981-1991:areaNat")
stroke.trim <- function(t1, t2) 
   subset(transform(stroke,
                    entry=t1, exit=pmin(t2, obsmonths), 
                    dead=dead & obsmonths <= t2), 
          entry < exit)
stroke2 <- do.call(rbind, mapply(stroke.trim, 
       c(0,0.5,2,12), c(0.5,2,12,Inf), SIMPLIFY=F))
summary(glm(dead~sex+age+factor(entry), poisson, 
       offset=log(exit-entry), data=stroke2))  
rm(list=ls())
while(search()[2] != "package:ISwR") detach()
