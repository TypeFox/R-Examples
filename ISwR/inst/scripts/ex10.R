attach(thuesen)
f <- cut(blood.glucose, c(4, 7, 9, 12, 20))
levels(f) <- c("low", "intermediate", "high", "very high")
bcmort2 <- within(bcmort,{
  period <- area <- cohort
  levels(period) <- rep(c("1991-2001","1981-1991"), each=2)
  levels(area) <- rep(c("Cph+Frb","Nat"),2) 
})
summary(bcmort2)
ashina.long <- reshape(ashina, direction="long", 
                    varying=1:2, timevar="treat")
ashina.long <- within(ashina.long, {
     m <- matrix(c(2,1,1,2),2)
     id <- factor(id)
     treat <- factor(treat)
     grp <- factor(grp)
     period <- factor(m[cbind(grp,treat)])
     rm(m)
})
within(ashina.long,
  period2 <- ifelse(treat != "active", 
             as.numeric(grp), 3 - as.numeric(grp))
)
stroke.trim <- function(t1, t2) 
   subset(transform(stroke,
                    entry=t1, exit=pmin(t2, obsmonths), 
                    dead=dead & obsmonths <= t2), 
          entry < exit)
stroke2 <- do.call(rbind, mapply(stroke.trim, 
       c(0,0.5,2,12), c(0.5,2,12,Inf), SIMPLIFY=F))
table(stroke$dead)        
table(stroke2$dead)        
rm(list=ls())
while(search()[2] != "package:ISwR") detach()
