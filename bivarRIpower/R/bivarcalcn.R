bivarcalcn <-
function(power,powerfor='RI',timepts,d1,d2,p,p1,s1,s2,r,r1){

### 1. Calculate correlation between bivariate outcomes under null / alternative hypotheses - concurrent and lagged ###

coryy_con0=(p*d1*d2 + r*s1*s2)/sqrt((d1^2 + s1^2)*(d2^2 + s2^2))
coryy_lag0=(p*d1*d2)/sqrt((d1^2 + s1^2)*(d2^2 + s2^2))

coryy_con1=(p1*d1*d2 + r1*s1*s2)/sqrt((d1^2 + s1^2)*(d2^2 + s2^2))
coryy_lag1=(p1*d1*d2)/sqrt((d1^2 + s1^2)*(d2^2 + s2^2))


### 2. Set up V (covariance matrix) ###

si1=s1^2*diag(timepts) 
si2=r*s1*s2*diag(timepts) 
si3=s2^2*diag(timepts) 
si=rbind(cbind(si1,si2),cbind(si2,si3))
ztop=   cbind(rep(1,timepts),rep(0,timepts))
zbottom=cbind(rep(0,timepts),rep(1,timepts))
z=matrix(rbind(ztop,zbottom),ncol=2)
d=matrix(c(d1^2,p*d1*d2,p*d1*d2,d2^2),ncol=2)

v = si + z%*%d%*%t(z)

vi=solve(v)


### 3. Calculate matrices of 1st derivatives of V with respect to variance parameters(d1,d2,p,s1,s2,r) ###

zeroblock=matrix(rep(0,timepts*timepts),nrow=timepts) 

sub11_d1=matrix(rep(2*d1,timepts*timepts),nrow=timepts)
sub12_d1=matrix(rep(p*d2,timepts*timepts),nrow=timepts)
dv_d1=(rbind(cbind(sub11_d1,sub12_d1),cbind(sub12_d1,zeroblock)))

sub12_d2=matrix(rep(p*d1,timepts*timepts),nrow=timepts)
sub22_d2=matrix(rep(2*d2,timepts*timepts),nrow=timepts)
dv_d2=(rbind(cbind(zeroblock,sub12_d2),cbind(sub12_d2,sub22_d2)))

sub12_p=matrix(rep(d1*d2,timepts*timepts),nrow=timepts)
dv_p=(rbind(cbind(zeroblock,sub12_p),cbind(sub12_p,zeroblock)))

sub11_s1=2*s1*diag(timepts)
sub12_s1=r*s2*diag(timepts)
dv_s1=(rbind(cbind(sub11_s1,sub12_s1),cbind(sub12_s1,zeroblock)))

sub12_s2=r*s1*diag(timepts)
sub22_s2=2*s2*diag(timepts)
dv_s2=(rbind(cbind(zeroblock,sub12_s2),cbind(sub12_s2,sub22_s2)))

sub12_r=s1*s2*diag(timepts)
dv_r=(rbind(cbind(zeroblock,sub12_r),cbind(sub12_r,zeroblock)))


### 4. Calculate information matrix ###

#Note: Leave out n/2 multiplier for now, include in 9.

#Step 5a: Matrices for information matrix.
id1d1 = vi %*% dv_d1 %*% vi %*% dv_d1
id1d2 = vi %*% dv_d1 %*% vi %*% dv_d2 
id1p  = vi %*% dv_d1 %*% vi %*% dv_p 
id1s1 = vi %*% dv_d1 %*% vi %*% dv_s1 
id1s2 = vi %*% dv_d1 %*% vi %*% dv_s2 
id1r  = vi %*% dv_d1 %*% vi %*% dv_r
id2d2 = vi %*% dv_d2 %*% vi %*% dv_d2
id2p  = vi %*% dv_d2 %*% vi %*% dv_p
id2s1 = vi %*% dv_d2 %*% vi %*% dv_s1  
id2s2 = vi %*% dv_d2 %*% vi %*% dv_s2
id2r  = vi %*% dv_d2 %*% vi %*% dv_r 
ipp   = vi %*% dv_p %*% vi %*% dv_p 
ips1  = vi %*% dv_p %*% vi %*% dv_s1
ips2  = vi %*% dv_p %*% vi %*% dv_s2 
ipr   = vi %*% dv_p %*% vi %*% dv_r
is1s1 = vi %*% dv_s1 %*% vi %*% dv_s1 
is1s2 = vi %*% dv_s1 %*% vi %*% dv_s2 
is1r  = vi %*% dv_s1 %*% vi %*% dv_r
is2s2 = vi %*% dv_s2 %*% vi %*% dv_s2 
is2r  = vi %*% dv_s2 %*% vi %*% dv_r
irr   = vi %*% dv_r %*% vi %*% dv_r
 

#Step 5: Traces of matrices from Step 4

tid1d1=0
tid1d2=0
tid1p =0
tid1s1=0
tid1s2=0
tid1r =0
tid2d2=0
tid2p =0
tid2s1=0
tid2s2=0
tid2r =0
tipp  =0
tips1 =0
tips2 =0
tipr  =0
tis1s1=0
tis1s2=0
tis1r =0
tis2s2=0
tis2r =0
tirr  =0


timepts2=2*timepts #note: 1:2*timepts, not from 1 to 2x number time points, need to create variable
for (i in 1:timepts2)
{
tid1d1 = tid1d1 + id1d1[i,i]
tid1d2 = tid1d2 + id1d2[i,i]
tid1p  = tid1p  + id1p[i,i]
tid1s1 = tid1s1 + id1s1[i,i]
tid1s2 = tid1s2 + id1s2[i,i]
tid1r  = tid1r  + id1r[i,i]
tid2d2 = tid2d2 + id2d2[i,i]
tid2p  = tid2p  + id2p[i,i]
tid2s1 = tid2s1 + id2s1[i,i]
tid2s2 = tid2s2 + id2s2[i,i]
tid2r  = tid2r  + id2r[i,i]
tipp   = tipp   + ipp[i,i]
tips1  = tips1  + ips1[i,i]
tips2  = tips2  + ips2[i,i]
tipr   = tipr   + ipr[i,i]
tis1s1 = tis1s1 + is1s1[i,i]
tis1s2 = tis1s2 + is1s2[i,i]
tis1r  = tis1r  + is1r[i,i]
tis2s2 = tis2s2 + is2s2[i,i]
tis2r  = tis2r  + is2r[i,i]
tirr   = tirr   + irr[i,i]
}

imatrix=matrix(rbind(
c(tid1d1,tid1d2,tid1p ,tid1s1,tid1s2,tid1r),
c(tid1d2,tid2d2,tid2p ,tid2s1,tid2s2,tid2r),
c(tid1p ,tid2p ,tipp  ,tips1 ,tips2 ,tipr ),
c(tid1s1,tid2s1,tips1 ,tis1s1,tis1s2,tis1r),
c(tid1s2,tid2s2,tips2 ,tis1s2,tis2s2,tis2r),
c(tid1r ,tid2r ,tipr  ,tis1r ,tis2r ,tirr )),nrow=6)



### 6. Invert information matrix ### 

invi=solve(imatrix)
invip=invi[3,3]   # Inverse info of p (random intercept correlation) under null hypothesis
invir=invi[6,6]   # Inverse info of r (residual correlation)         under null hypothesis



### 7. Multiply information matrix by derivative vectors to obtain variance of correlations on outcome pairs, concurrent and lagged ###

#Concurrent time points
dvector=numeric(6) 
eta=p*d1*d2 + r*s1*s2
omega=sqrt((d1^2 + s1^2)*(d2^2 + s2^2))

dvector[1]=p*d2/omega - d1*eta*(d2^2 + s2^2)/omega^3
dvector[2]=p*d1/omega - d2*eta*(d1^2 + s1^2)/omega^3
dvector[3]=d1*d2/omega
dvector[4]=r*s2/omega - s1*eta*(d2^2 + s2^2)/omega^3
dvector[5]=r*s1/omega - s2*eta*(d1^2 + s1^2)/omega^3
dvector[6]=s1*s2/omega

inviyy_con=t(dvector)%*%invi%*%dvector

#Lagged time points
dvectorl=numeric(6)
dvectorl[1]=p*d2/omega - p*d1^2*d2*(d2^2 + s2^2)/omega^3
dvectorl[2]=p*d1/omega - p*d1*d2^2*(d1^2 + s1^2)/omega^3
dvectorl[3]=d1*d2/omega
dvectorl[4]=-1*p*d1*d2*s1*(d2^2 + s2^2)/omega^3
dvectorl[5]=-1*p*d1*d2*s2*(d1^2 + s1^2)/omega^3
dvectorl[6]=0

inviyy_lag=t(dvectorl)%*%invi%*%dvectorl


### 8. Power calculations using normal approximation and variances from 1. - 5. 

if (powerfor=='RI') 
{
powerp=power
Zbetap=qnorm(1-powerp, mean = 0, sd = 1, lower.tail = TRUE, log.p = FALSE)
n = 2*invip*(1.96-Zbetap)^2 / (p-p1)^2
powerp=round(100*powerp,digits=1)

mur=1.96-sqrt(r1*r1/(invir*2/n))
powerr=round(100*(1-pnorm(mur, mean = 0, sd = 1, lower.tail = TRUE, log.p = FALSE)),digits=1)

muyycon=1.96-sqrt(coryy_con1*coryy_con1/(inviyy_con*2/n))
poweryycon=round(100*(1-pnorm(muyycon, mean = 0, sd = 1, lower.tail = TRUE, log.p = FALSE)),digits=1)

muyylag=1.96-sqrt(coryy_lag1*coryy_lag1/(inviyy_lag*2/n))
poweryylag=round(100*(1-pnorm(muyylag, mean = 0, sd = 1, lower.tail = TRUE, log.p = FALSE)),digits=1)
}

if (powerfor=='RESIDUAL') 
{
powerr=power
Zbetar=qnorm(1-powerr, mean = 0, sd = 1, lower.tail = TRUE, log.p = FALSE)
n = 2*invir*(1.96-Zbetar)^2 / (r-r1)^2
powerr=round(100*powerr,digits=1)

mup=1.96-sqrt(p1*p1/(invip*2/n))
powerp=round(100*(1-pnorm(mup, mean = 0, sd = 1, lower.tail = TRUE, log.p = FALSE)),digits=1)

muyycon=1.96-sqrt(coryy_con1*coryy_con1/(inviyy_con*2/n))
poweryycon=round(100*(1-pnorm(muyycon, mean = 0, sd = 1, lower.tail = TRUE, log.p = FALSE)),digits=1)

muyylag=1.96-sqrt(coryy_lag1*coryy_lag1/(inviyy_lag*2/n))
poweryylag=round(100*(1-pnorm(muyylag, mean = 0, sd = 1, lower.tail = TRUE, log.p = FALSE)),digits=1)
}

if (powerfor=='YYcon') 
{
poweryycon=power
Zbetayy=qnorm(1-poweryycon, mean = 0, sd = 1, lower.tail = TRUE, log.p = FALSE)
n = 2*inviyy_con*(1.96-Zbetayy)^2 / (coryy_con0-coryy_con1)^2

mup=1.96-sqrt(p1*p1/(invip*2/n))
powerp=round(100*(1-pnorm(mup, mean = 0, sd = 1, lower.tail = TRUE, log.p = FALSE)),digits=1)

mur=1.96-sqrt(r1*r1/(invir*2/n))
powerr=round(100*(1-pnorm(mur, mean = 0, sd = 1, lower.tail = TRUE, log.p = FALSE)),digits=1)

muyylag=1.96-sqrt(coryy_lag1*coryy_lag1/(inviyy_lag*2/n))
poweryylag=round(100*(1-pnorm(muyylag, mean = 0, sd = 1, lower.tail = TRUE, log.p = FALSE)),digits=1)
}

if (powerfor=='YYlag') 
{
poweryylag=power
Zbetayy=qnorm(1-poweryylag, mean = 0, sd = 1, lower.tail = TRUE, log.p = FALSE)
n = 2*inviyy_lag*(1.96-Zbetayy)^2 / (coryy_lag1-coryy_lag1)^2

mup=1.96-sqrt(p1*p1/(invip*2/n))
powerp=round(100*(1-pnorm(mup, mean = 0, sd = 1, lower.tail = TRUE, log.p = FALSE)),digits=1)

mur=1.96-sqrt(r1*r1/(invir*2/n))
powerr=round(100*(1-pnorm(mur, mean = 0, sd = 1, lower.tail = TRUE, log.p = FALSE)),digits=1)

muyycon=1.96-sqrt(coryy_con1*coryy_con1/(inviyy_con*2/n))
poweryycon=round(100*(1-pnorm(muyycon, mean = 0, sd = 1, lower.tail = TRUE, log.p = FALSE)),digits=1)
}


### 9. Variances obtained by inverting information matrices

varp=2*invip/n
varr=2*invir/n
varyycon=2*inviyy_con/n
varyylag=2*inviyy_lag/n


### 10. Output results to screen

textvarpar1  ="Variance parameters        "
textvarpar2  ="Clusters                 = "
textvarpar3  ="Repeated measurements    = "
textvarpar4  ="Standard deviations        "
textvarpar5  =" 1st random intercept    = "
textvarpar6  =" 2nd random intercept    = " 
textvarpar7  =" 1st residual term       = "
textvarpar8  =" 2nd residaul term       = "
textvarpar9  ="Correlations               "
textvarpar10 =" RI under H_o            = "
textvarpar11 =" RI under H_a            = "
textvarpar12 =" Residual under H_o      = "
textvarpar13 =" Residual under H_a      = "
textvarpar14 =" Con obs  under H_o      = "
textvarpar15 =" Con obs  under H_a      = "
textvarpar16 =" Lag obs  under H_o      = "
textvarpar17 =" Lag obs  under H_a      = "

textvarpar18 ="Correlation variances under H_o "

textpower    ="Power (%) for correlations "
textdash     ="------------------------------- "
textpowerp   ="Random intercept         = "
textpowerr   ="Residual                 = "
textpowerc   ="Concurrent observations  = "
textpowerl   ="Lagged observations      = "

cat(textvarpar1,"\n",textdash,"\n",textvarpar2,round(n,digits=1),"\n",textvarpar3,timepts,"\n",textvarpar4,"\n",
textvarpar5,d1,"\n",textvarpar6,d2,"\n",textvarpar7,s1,"\n",textvarpar8,s2,"\n",textvarpar9,"\n",
textvarpar10,p,"\n",textvarpar11,p1,"\n",textvarpar12,r,"\n",textvarpar13,r1,"\n",
textvarpar14,coryy_con0,"\n",textvarpar15,coryy_con1,"\n",textvarpar16,coryy_lag0,"\n",textvarpar17,coryy_lag1,"\n","\n",
textvarpar18,"\n",textdash,"\n",textpowerp,varp,"\n",textpowerr,varr,"\n",textpowerc,varyycon,"\n",textpowerl,varyylag,"\n","\n",
textpower,"\n",textdash,"\n",textpowerp,powerp,"%","\n",textpowerr,powerr,"%","\n",textpowerc,poweryycon,"%","\n",
textpowerl,poweryylag,"%","\n",sep="")
}