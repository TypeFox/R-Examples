IH.summary <-
function(dd,L,p=0.95,gam=0.95,bcol=NA){
# from  allss  2 Oct04 Revised from 7 Jul 04
#  
#  INPUT:
#    dd is 2col data frame or matrix  with
#    x  in column 1  and det=1 for  detect 0 for nondetect column 2
#    gam is confidence level (one sided) for confidence intervals
#    p is quantile for UTL-p-gam
#    bcol is column number of by variable (if any)
# REQUIRES lnorm.ml()  efclnp() kmms() efraction.ml()
#
ihsum<-function(dd,L,p,gam){

nt <- length(dd[,1])
ndet <- sum(dd[,2])
if( ndet < 2) stop("At least 2 detects required for most summary statistics")
PCndet <- round( 100*(nt-ndet)/nt,1)
du <- dd[,1:2] 

ys <- lnorm.ml(du)
GM <- exp(ys$mu)
GSD <- exp(ys$sigma)

#  Method 1 is equivalent to modified  Cox direct Method
#  approximate CLs for lognormal mean
#   m = number of non-detects
lex <- ys$logEX  # same as mu + 0.5*sig^2
tv <- qt(gam, ys$m-1)
   EX<-exp(ys$logEX)
   EXL<-exp( lex - tv*ys$se.logEX)
   EXU<-exp( lex + tv*ys$se.logEX)
#  calculate Product limit estimate of CDF for left censored data
#  RH is (progressively left censored version of correlation coef)
#  see Verrill and Johnson JASA 1988 Equation 4.4
#  if no censoring Rsq= RH^2 will be approx equal to Shapiro-Wilk W
pl<-plend(dd)
  # calculate "plotting position"
nd<- dim(pl)[1]
   pp<-1:nd 
   pp[1]<-pl$p[1]
   for (j in 2:nd) pp[j]<- (pl$p[j-1] + pl$p[j])/2.0
Rsq <- cor(qnorm(pp),log(pl$a) )^2
#  calculate UTL and CLs
xcl<- as.numeric(percentile.ple(dd,p,gam,TRUE))[1]  # based on PLE
xp<- percentile.ml(ys,p,gam,FALSE)[1:3]    # baed on ML
Xpo<- list(Xp.obs= xcl)	#  pth percentile of x
#
zL<- (log(L)-ys$mu)/ys$sig; noel<-paste("z_L_",L,sep="")
xmax<-max(du[,1])
m5 <- nptl(nt,p,gam)
if( is.na(m5) )NPTL<-NA    else NPTL<- rev(sort(du[,1]))[m5]
#  calculate Kapla-Meier mean se.mean and CLs
km<-kmms(dd,gam)
par1<-ys[c(1,5,2,6)]
par2<- list("GM"=GM,"GSD"=GSD,"EX"=EX,"EX.LCL"=EXL,"EX.UCL"=EXU)
out<-c(par1,par2,km[1:4],Xpo,xp)
#  calculate excedance fraction for L using lognormal model
#    approximate large sample results
ef<- efraction.ml(ys,gam,L,dat=FALSE)
#  non-parametric excedance fraction and CLs for L
npf<- efclnp(du,gam,L )

par3<-list("NpUTL"=NPTL,"Maximum"=xmax,"NonDet%"=round(PCndet,2),
  "n"=nt,"Rsq"=Rsq,"m"=ndet)
par4<-list("m2logL"=ys$m2,"L"=L,"p"=p,"gamma"=gam)
out<-c(out,par3,ef[1:3],npf[1:3],par4)  
out
}

if(is.na(bcol) ){stats<- cbind(unlist(ihsum(dd[,1:2],L,p,gam ))) 
        cname <- deparse(substitute(dd)); rn<- rownames(stats)}
else{
gval<- sort( unique(dd[,bcol])  ) 
cname<- as.character( gval  )
t1<-dd[dd[,bcol]==gval[1],1:2]
t1<- cbind( unlist(ihsum(t1,L,p,gam)) )

rn<-  rownames(t1)
stats<- cbind(t1)
	jj<-2
    while(jj <= length(gval) ){
    t1<-dd[dd[,bcol]==gval[jj],1:2]
    t1<- unlist(ihsum(t1,L,p,gam))
    jj<- jj + 1
    stats<-cbind(stats,t1)
    }
}
#stats<- data.frame(stats)
dimnames(stats)<-list(rn,cname)
stats
}

