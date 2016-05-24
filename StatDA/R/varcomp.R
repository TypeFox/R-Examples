varcomp <- function(a1,a2,f1,f2){
# PF, 28.9.2006
# estimate the variance components for ANOVA
# a1,a2 ... analytical duplicates
# f1,f2 ... field duplicates
#
a=na.omit(cbind(a1,a2))
a1=a[,1]
a2=a[,2]
alen <- length(a1)

adiff <- a1 - a2
sams <- sum(adiff * adiff)/(2 * alen)
tms <- var(as.vector(a))
tdf <- alen * 2 - 1
tss <- tms * tdf
wss <- sams * alen
bss <- tss - wss
bdf <- alen - 1
bms <- bss/bdf
bvar <- (bms - sams)/2
tvar <- bvar + sams

f=na.omit(cbind(f1,f2))
f1=f[,1]
f2=f[,2]
flen <- length(f1)
fdiff <- f1 - f2
sams.f <- sum(fdiff * fdiff)/(2 * flen)
tms.f <- var(as.vector(f))
tdf.f <- flen * 2 - 1
tss.f <- tms.f * tdf.f
wss.f <- sams.f * flen
bss.f <- tss.f - wss.f
bdf.f <- flen - 1
bms.f <- bss.f/bdf.f
bvar.f <- (bms.f - sams.f)/2
tvar.f <- bvar.f + sams.f
if ((tvar.f==0) || (tvar==0)){
  warning("Field or analytical variation zero!")
  out=c(NA,NA,NA,NA)
}
else{
# F-Test for euqality of variances
pval=1-pf(bvar/bvar.f,bdf,bdf.f)
if (pval<0.05) {
  warning("Variances not equal!")
#  tvar.f=(tvar.f+tvar)/2
}

pct.regional <- (100 * bvar.f)/tvar.f
pct.site <- (100 * (sams.f-sams))/tvar.f
pct.analytical <- (100 * sams)/tvar.f

out=c(pct.regional,pct.site,pct.analytical,pval)
#list(pct.regional=pct.regional,pct.site=pct.site,pct.analytical=pct.analytical,
#    p.value=pval)
}
out
}

