threewayanova <-
function(Y,n,m,p){

Y=as.matrix(Y)
ss=SUM(Y)$ssq
#	compute grand mean
mm=mean(Y)
ssm=n*m*p*mm^2
#	compute residuals after subtraction of grand mean
Y=Y-mm

#	compute main effects
ma=SUM(Y)$mr			
Y=permnew(Y,n,m,p)
mb=SUM(Y)$mr
Y=permnew(Y,m,p,n)
mc=SUM(Y)$mr
Y=permnew(Y,p,n,m)
ssa=m*p*sum(ma^2)
ssb=n*p*sum(mb^2)
ssc=n*m*sum(mc^2)

#	compute residuals after subtraction of main effects
Y=Y-matrix(ma,n,m*p)		
Y=permnew(Y,n,m,p)-matrix(mb,m,n*p)
Y=permnew(Y,m,p,n)-matrix(mc,p,n*m)
Y=permnew(Y,p,n,m)

#	compute second order interactions
mbc=SUM(Y)$mc
Y=permnew(Y,n,m,p)
mac=SUM(Y)$mc
Y=permnew(Y,m,p,n)
mab=SUM(Y)$mc
Y=permnew(Y,p,n,m)
ssab=p*sum(mab^2)
ssac=m*sum(mac^2)
ssbc=n*sum(mbc^2)

#	compute residuals after subtraction of second order interacions
Y=Y-t(matrix(mbc,m*p,n))
Y=permnew(Y,n,m,p)
Y=Y-t(matrix(mac,n*p,m))
Y=permnew(Y,m,p,n)
Y=Y-t(matrix(mab,n*m,p))
Y=permnew(Y,p,n,m)
ssabc=SUM(Y)$ssq
ss=ss-ssm
cat(paste("Total ssq after subtraction of grand mean  = ", ss), fill=TRUE)
cat(paste("SS_a        =   ",round(ssa,digits=6)       ,"(", round(ssa/ss*100,digits=2) ,")"),fill=TRUE)
cat(paste("SS_b        =   ",round(ssb,digits=6)       ,"(", round(ssb/ss*100,digits=2) ,")"),fill=TRUE)
cat(paste("SS_c        =   ",round(ssc,digits=6)       ,"(", round(ssc/ss*100,digits=2) ,")"),fill=TRUE)
cat(paste("SS_ab       =   ",round(ssab,digits=6)      ,"(", round(ssab/ss*100,digits=2) ,")"),fill=TRUE)
cat(paste("SS_ac       =   ",round(ssac,digits=6)      ,"(", round(ssac/ss*100,digits=2) ,")"),fill=TRUE)
cat(paste("SS_bc       =   ",round(ssbc,digits=6)      ,"(", round(ssbc/ss*100,digits=2) ,")"),fill=TRUE)
cat(paste("SS_abc      =   ",round(ssabc,digits=6)     ,"(", round(ssabc/ss*100,digits=2) ,")"),fill=TRUE)

out=list()
out$SS.a=ssa
out$SS.b=ssb
out$SS.c=ssc
out$SS.ab=ssab
out$SS.ac=ssac
out$SS.bc=ssbc
out$SS.abc=ssabc
return(out)
}
