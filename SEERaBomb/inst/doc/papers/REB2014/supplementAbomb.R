rm(list=ls()) 
# The first step is to show that there are lower AIC EAR models than the CML ERR
# and EAR models in: Hsu, W. L., Preston, ... and Mabuchi, K. The Incidence of
# Leukemia, Lymphoma ... among Atomic Bomb Survivors: 1950-2001. Radiat. Res.
# 179, 361-382 (2013). To use this script please place lsshempy.csv in
# /data/abomb or change this line
d<-read.csv("~/data/abomb/lsshempy.csv",header=T)
head(d,2)
d=d[d$mar_an>=0,] # take only those with known doses 
sum(d$cml) # check to see that we have 75 cases of CML
# viewing d as an environment, use within to remove columns not needed here
d=within(d,{rm(distcat,agxcat, agecat, dcat, time,upyr,subjects,year,
               nhl,hl,mye,all,oll,alltot,cll,hcl,clltot,atl,aml,oml,amol,amltot,
               othleuk,noncll,leuktot,hldtot,mar_ag,mar_an)}) 
head(d,2)
# the next line uses within to transform columns, as in the beginning of the
# epicure log file in Hsu et al's supplement 3
d=within(d, {nic = as.numeric(gdist > 12000) # nic= not in city
             over4gy = 1 - un4gy  
             rm(gdist,un4gy)  # get rid of things no longer needed as we go
             tsx = (age -agex) 
             lt25 = log(tsx/25) 
             a = log(age/70) 
             a55 = log(age/55) ;
             age55=age-55
             agex30=agex-30
             ax30 = log(agex/30) ;
             py10k = pyr/10000 ;
             py = pyr;
             s=sex-1
             c=city-1
             sv=mar_ad10/1000;
             hiro = as.numeric(city == 1) 
             naga = as.numeric(city ==2)
             rm(sex,city,pyr,mar_ad10)
})


head(d)
bk=d[d$sv<0.05,]
library(bbmle)
# ERR model of CML in Table S2 (in Supplement 1) is:
# back=0.22*exp(-0.06*s + (1.38 + 1.75*s)*a)
# Induced=(5.24*d)*exp(-1.50*c - 1.59*lt25 - 1.42*a55)
# where
# c = 0 for Hiroshima and 1 for Nagasaki        
# s = 0 for men and 1 for women                 
# a = ln(age/70) 
# a55 = ln(age/55)              
# lt25 = ln(tsx/25) 

# first compare their polynomial in age background model to our preferred exponential in age model
summary(pow1<-mle2(cml~dpois(lambda=py10k*A*exp(cs*s+(ca+csa*s)*a)),
                     start=list(A=.22,cs=-0.06,ca=1.38,csa=1.75),data=bk) ) 
summary(exp1k<-mle2(cml~dpois(lambda=py10k*A*exp(cs*s+k*age55)),
                     start=list(A=.22,cs=-0.6,k=0.04),data=bk) ) 
summary(exp2k<-mle2(cml~dpois(lambda=py10k*A*exp(cs*s+(k+ks*s)*age55)),
                 start=list(A=.22,cs=-0.6,k=0.03,ks=.02),data=bk) ) 
AIC(pow1,exp1k,exp2k)
ICtab(pow1,exp1k,exp2k)
# AIC favors exp(k*age) over  age^k for background; also, sex diff in k is not justified
BIC(pow1,exp1k,exp2k) # BIC seems to penalize too much for adding a parameter

# now add in their ERR induced term
# Induced=(5.24*d)*exp(-1.50*c - 1.59*lt25 - 1.42*a55)  [according to Table S2 in Supplement 1], but 
# with extra terms according to the epicure log file in Supplement 3. 
# Since the Supplement 3 model exactly reproduces Table S9 in Supplement 2, we go with it:
summary(err<-mle2(cml~dpois(lambda=py10k*(A1*exp(c1s*s+c1a*a+c1sa*a*s+c1nic0*hiro*nic+c1nic1*naga*nic)*
                          (1+sv*A2*exp(c2c*c + c2t*lt25 + c2a*a55+c2o4g*over4gy)   ) ) ),
                 start=list(A1=.22,c1s=-0.06,c1a=1.38,c1sa=1.75,c1nic0=-0.2,c1nic1=-0.86,
                            A2=5.24,c2c=-1.50,c2t=-1.59,c2a=-1.42,c2o4g=-0.3),data=d) ) 
# note that many parameters are highly insignificant by the Wald test, so this looks overparameterized!

AIC(err)
logLik(err)
deviance(err)
prd=predict(err)
-2*sum(d$cml*log(prd))
-2*sum(d$cml*log(prd))+2*length(coef(err)) # 11-parameter model => diff of 22
# [1] 754.454 # exactly as in epicure log file of supplement. In computing the 
# deviance here, I'm OK with ditching the third term below since it only depends
# on the data, and since it vanishes anyway since the cells are all 0 or 1 (in this data)
# and the factorial takes these all to 1 which the log then takes to zero, but was wondering
# how they can leave out the second term, i.e. -prd
-2*sum(d$cml*log(prd)-prd-log(factorial(d$cml)))
-2*sum(d$cml*log(prd)-prd)
sum(prd) # my only guess is that the sum of the expected means = sum of observed cases
# so that it too really only depends on the data, at the level of the whole sum. 

# EAR model of CML. From Table S2 in Supplement 2
# back=0.25*exp(-0.28*s + 1.60*a + 1.38*(s*a))
# Induced=
# (0.68*d)*exp(-1.49*c - 1.63*lt25 - 0.20*a55) men
# (0.57*d)*exp(-1.49*c - 1.63*lt25 + 2.10*a55) women
# but the epicure log shows more to the model, as follows
summary(ear<-mle2(cml~dpois(lambda=py10k*(A1*exp(c1s*s+c1a*a+c1sa*a*s+c1nic0*hiro*nic+c1nic1*naga*nic)
                        +sv*A2*exp(c2c*c + c2s*s+c2t*lt25 + c2a*a55+c2sa*s*a55+c2o4g*over4gy)   ) ) ,
                 #                                        method="L-BFGS-B",
                 start=list(A1=.22,c1s=-0.06,c1a=1.38,c1sa=1.75,c1nic0=-0.2,c1nic1=-0.86,
                            A2=5.24,c2c=-1.50,c2s=-0.18,c2t=-1.59,c2a=-1.42,c2sa=2.3,c2o4g=-0.3),data=d) ) 
# still lots of insignificant parameters.
prd=predict(ear)
-2*sum(d$cml*log(prd))
-2*sum(d$cml*log(prd))+2*length(coef(ear)) # 13-parameter model, so add 26 => worse model 
# 753.7032  # also exactly as in the epicure log file supplement => a better fit than the ERR model

AIC(err,ear)
ICtab(err,ear) # err wins by AIC
anova(err,ear) # err wins by likelihood ratio test (i.e. deviance difference test)


################################################################################3
# The following EAR model beats out Hsu's ERR model (in terms of AICs)  
summary(s1<-mle2(cml~dpois(lambda=py*(exp(c10+c2c*c+c2s*s+k*age+ks*s*age)
                                      +sv*exp(c20+c2c*c + c2s*s+c2t*tsx+c2ts*s*tsx)  ) ),
         start=list(c10=-12.4,k=0.025,ks=0,c20=-10.5,c2c=-1.50,c2s=-0.18,c2t=-0.4,c2ts=0.3),data=d) ) 
prd=predict(s1)
-2*sum(d$cml*log(prd))
-2*sum(d$cml*log(prd))+2*length(coef(s1)) # 8 parameters, so add 16
ICtab(err,s1) # s1 wins by AIC
anova(err,s1) # s1 wins by likelihood ratio test (i.e. deviance difference test)


# # using this model there is no evidence for the need for kerma>4gy correction
# summary(s1ov<-mle2(cml~dpois(lambda=py*(exp(c10+c2c*c+c2s*s+k*age+ks*s*age)
#                                       +sv*exp(c20+c2c*c + c2s*s+c2t*tsx+c2ts*s*tsx+c2o4g*over4gy)  ) ),
#                  start=list(c10=-12.4,k=0.025,ks=0,c20=-10.5,c2c=-1.50,
#                             c2s=-0.18,c2t=-0.4,c2ts=0.3,c2o4g=-0.3),data=d) ) 
# prd=predict(s1ov)
# -2*sum(d$cml*log(prd))
# -2*sum(d$cml*log(prd))+2*length(coef(s1ov))
# 

# checking if adding lage55 in the induced term helps, yes, it does 
summary(s1a55<-mle2(cml~dpois(lambda=py*(exp(c10+c2c*c+c2s*s+k*age+ks*s*age)
                     +sv*exp(c20+c2c*c + c2s*s+c2t*tsx+c2ts*s*tsx + c2a*a55+c2sa*s*a55 )  ) ),
                    start=list(c10=-12.4,k=0.025,ks=0,c20=-10.5,c2c=-1.50,
                               c2s=-0.18,c2t=-0.4,c2ts=0.3,c2a=-1.42,c2sa=2.3),data=d) ) 
prd=predict(s1a55)
-2*sum(d$cml*log(prd))                       
-2*sum(d$cml*log(prd))+2*length(coef(s1a55)) # lage55 terms are justified! new best model
ICtab(err,s1a55) # new model wins by AIC diff of 3.4
anova(err,s1a55) # new model has less deviance and fewer parameters!

# keeping only the female component of lage55  
summary(s1a55f<-mle2(cml~dpois(lambda=py*(exp(c10+c2c*c+c2s*s+(k+ks*s)*age)
                         +sv*exp(c20+c2c*c + c2s*s+c2t*tsx+c2ts*s*tsx +c2sa*s*a55 )  ) ),
                    start=list(c10=-12.4,k=0.025,ks=0,c20=-10.5,c2c=-1.50,
                               c2s=-0.18,c2t=-0.4,c2ts=0.3,c2sa=2.3),data=d) ) 

prd=predict(s1a55f)
-2*sum(d$cml*log(prd))                       
-2*sum(d$cml*log(prd))+2*length(coef(s1a55f)) # new best model 
ICtab(err,s1a55f) # now new model wins by AIC diff of 5.3
anova(err,s1a55f) # deviance loss was negligible, consistent with AIC drop of almost 2

# # do we still need a sex diff in exponential decays?  
# summary(m1a55f<-mle2(cml~dpois(lambda=py*(exp(c10+c2c*c+c2s*s+k*age+ks*s*age)
#                                           +sv*exp(c20+c2c*c + c2s*s+c2t*tsx +c2sa*s*a55 )  ) ),
#                      start=list(c10=-12.4,k=0.025,ks=0,c20=-10.5,c2c=-1.50,
#                                 c2s=-0.18,c2t=-0.4,c2sa=2.3),data=d) ) 
# prd=predict(m1a55f)
# -2*sum(d$cml*log(prd))                       
# -2*sum(d$cml*log(prd))+2*length(coef(m1a55f)) #  yes, keep the difference in exponential decays

# # try switching to age at exposure in induced term doesn't help
# summary(s1a55x<-mle2(cml~dpois(lambda=py*(exp(c10+c2c*c+c2s*s+(k+ks*s)*age)
#                                           +sv*exp(c20+c2c*c + c2s*s+c2t*tsx+c2ts*s*tsx +c2sa*s*agex )  ) ),
#                      start=list(c10=-12.4,k=0.025,ks=0,c20=-10.5,c2c=-1.50,
#                                 c2s=-0.18,c2t=-0.4,c2ts=0.3,c2sa=2.3),data=d) ) 
# prd=predict(s1a55x)
# -2*sum(d$cml*log(prd))                       
# -2*sum(d$cml*log(prd))+2*length(coef(s1a55x)) 
# ICtab(s1a55f,s1a55x) # fit is worse

# make exponential more realistic by bringing tsx out in front of induced term   
summary(s1tsx<-mle2(cml~dpois(lambda=py*(exp(c10+c2c*c+c2s*s+(k+ks*s)*age)
                                          +sv*tsx*exp(c20+c2c*c + c2s*s+c2t*tsx+c2ts*s*tsx +c2sa*s*a55 )  ) ),
                     start=list(c10=-12.4,k=0.025,ks=0,c20=-10.5,c2c=-1.50,
                                c2s=-0.18,c2t=-0.4,c2ts=0.3,c2sa=2.3),data=d) ) 
# note that the parameter P values are now almost all signficant
prd=predict(s1tsx)
-2*sum(d$cml*log(prd))                       
-2*sum(d$cml*log(prd))+2*length(coef(s1tsx)) # new best model 
ICtab(s1tsx,s1a55f) # tsx in front wins by AIC diff of 0.4


# what about the age at exposure hump for HM only  
a55HMb<-mle2(cml~dpois(lambda=py*(exp(c10+c2c*c+c2s*s+k*age+ks*s*age)    +
      sv*tsx*exp(c20+c2c*c + c2s*s+c2t*tsx+c2ts*s*tsx +c2sa*s*a55+ c2x*(1-c)*(1-s)*abs(agex-30))  ) ),
                     start=list(c10=-12.4,k=0.025,ks=0,c20=-10.5,c2c=-1.50,
                                c2s=-0.18,c2t=-0.4,c2ts=0.3,c2x=.1,c2sa=2.3),data=d)  
prd=predict(a55HMb)
-2*sum(d$cml*log(prd))                       
-2*sum(d$cml*log(prd))+2*length(coef(a55HMb)) # 

summary(a55HMb)     #P=0.23 on c2x param is consistent with 0.21 below for model diffs
AIC(a55HMb,s1tsx)
ICtab(a55HMb,s1tsx) #AIC difference of 0.4 => deviance better by 1.6 (i.e. maybe something to it)
anova(a55HMb,s1tsx) # P=0.21 => HM hump at agex=30 is not statistically significant

#rewrite to estimate sex difference in EAR time constants 
summary(best<-mle2(cml~dpois(lambda=py*(exp(c10+c2c*c+c2s*s+k*age+ks*s*age)
                                        +sv*tsx*exp(c20+c2c*c + c2s*s-tsx/(c2t+c2ts*s) +c2sa*s*a55 )  ) ),
                     start=list(c10=-12.4,k=0.025,ks=0,c20=-10.5,c2c=-1.50,
                                c2s=-0.18,c2t=5,c2ts=10,c2sa=2.3),data=d) ) 
prd=predict(best)
-2*sum(d$cml*log(prd))                       
-2*sum(d$cml*log(prd))+2*length(coef(best))  
# see that it still is the best
AIC(s1tsx,best) # i.e. these should be identical, and they are (this is just a reparameterization)
ICtab(err,best) 
anova(err,best) 

# show that can't remove ks by AIC, but perhaps by LRT and BIC
summary(best) # see that ks was not significant by Wald Test
summary(bestT<-mle2(cml~dpois(lambda=py*(exp(c10+c2c*c+c2s*s+k*age)
                                        +sv*tsx*exp(c20+c2c*c + c2s*s-tsx/(c2t+c2ts*s) +c2sa*s*a55 )  ) ),
                   start=list(c10=-12.4,k=0.025,c20=-10.5,c2c=-1.50,
                              c2s=-0.18,c2t=5,c2ts=10,c2sa=2.3),data=d) ) 
prd=predict(bestT)
-2*sum(d$cml*log(prd))                       
anova(best,bestT) 
BICtab(best,bestT)  # LRT and Bayes both say take out ks 
ICtab(best,bestT)  # but AIC says leave it in. 

# Since I would prefer to see agex over age in my induced term, I tried things
# like this that failed
# d$oldF=d$s*as.numeric(d$agex>40)
# head(d)
# summary(bestT<-mle2(cml~dpois(lambda=py*(exp(c10+c2c*c+c2s*s+k*age+ks*s*age)
#                                         +sv*tsx*exp(c20+c2c*c + c2s*s-tsx/(c2t+c2ts*s) +c2sa*oldF )  ) ),
#                    start=list(c10=-12.4,k=0.025,ks=0,c20=-10.5,c2c=-1.50,
#                               c2s=-0.18,c2t=5,c2ts=log(10),c2sa=2.3),data=d) ) 
# prd=predict(bestT)
# -2*sum(d$cml*log(prd))                       
# ICtab(bestT,best) # all of my attempts to get age out of the induced term failed
# anova(bestT,best) 

#OK, so we stop at this point and now make tables and figures for the paper
(sb=summary(best))
(outs=round(coef(sb),4))
rownames(outs)[c(1,4:8)]<-c("c1hm","c2hm","cn","cs","tau","taus")
outs  # parameter names mapped to those in paper before sending to the table
library(hwriter)
# setwd("/users/radivot/downloads/sachs")
# hwrite(outs,"table2.html")

sd=coef(sb)["c2ts",2]
mn=coef(best)["c2ts"]
(tauDiff=round(2*c(mn,mn-1.96*sd,mn+1.96*sd),2)) # mean is twice tau for gamma with t^1 in front
sprintf("%s (%s, %s)",tauDiff[1],tauDiff[2],tauDiff[3])
(Taum=coef(best)["c2t"])  #so twice these are the mean male and female waiting times
(Tauf=Taum + mn)

# now get amplitude parameter ratio (M/F) estimates and CI 
sdA=coef(sb)["c2s",2]
mnA=coef(best)["c2s"]
(MovF=round((Taum/Tauf)^2*exp(-c(mnA,mnA+1.96*sdA,mnA-1.96*sdA)),2))
sprintf("%s (%s, %s)",MovF[1],MovF[2],MovF[3])
# note that the (Taum/Tauf)^2 is needed because the normalized density
# function is (1/Tau^2)*t*exp(-t/Tau)

###############################################################################

# The next goal is to make something like the dots in Fig 4B of Hsu et al, but for tsx decays
# To do this we replace the tsx exponential decays with a "nonparameteric" model (not distribution),
# i.e. this is non-parameteric like a one-way anova is of a time course. 
d=transform(d,tsxf=cut(tsx,6),sex=as.factor(s))
head(d)
summary(s3<-mle2(cml~dpois(lambda=py*(exp(c10 + c2c*c + c2s*s + k*age + ks*s*age)
                                +sv*exp(c2c*c  + f + c2sa*s*a55 )  ) ),
                 parameters=list(f~-1 + tsxf:sex),
                 start=list(c10=-12.4,k=0.025,ks=0,c2c=-1.50,
                              c2s=-0.18,f=-5,c2sa=2),data=d) ) 
prd=predict(s3)
-2*sum(d$cml*log(prd))  # fit is much improved, with deviance now at 746.07 (lowest of all)
-2*sum(d$cml*log(prd))+2*length(coef(s3))# but the extra parameters are not justified: this is OK, 
# the point here is to plot the points, not say that this is a better model. 
AIC(s3,best)
ICtab(s3,best) 
anova(s3,best) 

ss3=summary(s3)
wait=exp(coef(ss3)[6:17,1])*1e4
waitL=exp(coef(ss3)[6:17,1]-1.96*coef(ss3)[6:17,2])*1e4
waitU=exp(coef(ss3)[6:17,1]+1.96*coef(ss3)[6:17,2])*1e4
(lvls=levels(d$tsxf))
(Lvls=strsplit(lvls,","))
substring("(4.23",2)
(lows=sapply(Lvls,function(x) as.numeric(substring(x[1],2))))
(ups=sapply(Lvls,function(x) as.numeric(substring(x[2],1,4))))
(mids=round(apply(rbind(lows,ups),2,mean),2))
(dfc=data.frame(mids,wait,waitL,waitU,Sex=gl(2,6,labels=c("Male","Female") ) ) )

graphics.off()
if(length(grep("linux",R.Version()$os))) windows <- function( ... ) X11( ... )
if(length(grep("darwin",R.Version()$os))) windows <- function( ... ) quartz( ... )
windows(width=7,height=6)

library(ggplot2)
pd <- position_dodge(1) 
(p=ggplot(dfc, aes(x=mids, y=wait, shape=Sex, col=Sex, ymax=10)) + 
   geom_point(size=6,position=pd) ) # + ylim(0, 10) ) # scale_y_log10(limits=c(1e-5,10)) 
(p=p+labs(title="IR-to-CML Latency",x="Years since exposure",
          y=expression(paste("Cases per ",10^4," Person-Year-Sv") ) ) )    
(p=p+geom_errorbar(aes(ymin=waitL, ymax=waitU),width=.01,position=pd))
(p=p+theme(plot.title = element_text(size = rel(2.3)),
           axis.title = element_text(size = rel(2.3)),
           axis.text = element_text(size = rel(2.3)))  )
(p=p+theme(legend.position = c(0.8, .5), 
           legend.title = element_text(size = rel(2)) ,
           legend.text = element_text(size = rel(2))  ) )  

waitm=exp(coef(ss3)[6:11,1])
waitf=exp(coef(ss3)[12:17,1])
(MovF<-format(sum(waitm)/sum(waitf),digits=3)) # ratio of sums (like AUCs)

# now we still need the M/F CI, 
mns=coef(ss3)[6:17,1]
sds=coef(ss3)[6:17,2]
X=exp(matrix(rnorm(12*100000,mns,sds),byrow=T, ncol=12))
M=X[,1:6]; F=X[,7:12]
Ms=apply(M,1,sum)
Fs=apply(F,1,sum)
MovF=Ms/Fs
MF=round(quantile(MovF,c(0.5,0.025,0.975)),2)
(MF=sprintf("%s (%s, %s)",MF[1],MF[2],MF[3]))

pm=waitm/sum(waitm)
pf=waitf/sum(waitf)
(taum<-round(mids%*%pm,2))
(tauf<-round(mids%*%pf,2))

(p=p+annotate("text",x=15,y=10, hjust=0, label = paste("M/F =",MF),size=9) )
(lb1=paste0("tau[m] == ",taum,"*~Yrs"))
(p=p+annotate("text",x=15,y=9, hjust=0, label = lb1,size=9,parse=T) )
(lb1=paste0("tau[f] == ",tauf,"*~Yrs"))
(p=p+annotate("text",x=15,y=8, hjust=0, label = lb1,size=9,parse=T) )
# This next block calculates the 95% CI of the difference in waiting times
rdiscrete <- function(n, probs,values) {
  cumprobs <- cumsum(probs)
  singlenumber <- function() {
    x <- runif(1)
    N <- sum(x > cumprobs)
    N
  }
  values[replicate(n, singlenumber())+1]
}
x=rdiscrete(1e5,pf,mids)
summary(as.factor(x))
y=rdiscrete(1e5,pm,mids)
summary(as.factor(y))
(delT=quantile(x-y,c(0.025,0.975)))
(lb1=paste0("Delta*tau == ",tauf-taum,"(",delT[1],", ",delT[2],")"))
#(lb1=paste0("Delta*tau == ",tauf-taum,"~group(\"(\",",delT[1],",",delT[2],",\")\")")) # attempt with group() failed
# quartz(height=7,width=8)  # set back to square
(p=p+annotate("text",x=15,y=7, hjust=0, label = lb1,size=7,parse=T) ) #just show since can't add space before ( 
sprintf("%s (%s, %s)",tauf-taum,delT[1],delT[2])
# ggsave(p,file="~/IR2CML.eps")


### Figure 7: interpretations ##############
# graphics.off()
# quartz(height=6,width=6) 

# setEPS()
# postscript('~/concept.eps', height = 6, width = 6)

k=0.025
Tp=22.1
Rp=1.73
Tf=-15:23
fR=function(Tf) exp(-k*(Tf-Tp))
R=fR(Tf)


par(mar=c(4.5,4.5,0,.5))
plot(Tf,R,type="l",lwd=2.5,cex.lab=1.6,cex.axis=1.5,
     ylab="(male risk)/(female risk) M/F",
     xlab="extra female latency time T in years",axes=F,font=2)
mtext(side=3,line=-4,"Continuum of SEER CML\n Sex Difference Interpretations",cex=1.5,font=1)
# axis(1);axis(2)
points(x=c(0,Tp),y=c(Rp,1),pch=1,cex=3,lwd=3)
abline(h=1,v=0,lty=3)
fT=function(R) Tp-log(R)/k
# Tx=3.86
# Rx=1.26
# # points(x=c(fT(Rx),Tx),y=c(Rx,fR(Tx)),pch="+",cex=3)
# # plot(1:20,pch=1:20)
# x=seq(Tx,fT(Rx),0.1)
# x=seq(delT[1],delT[2],0.1)
# y=fR(x)
# points(x,y,type="l",lwd=6)

x=seq(tauDiff[2],tauDiff[3],0.1)
y=fR(x)
points(x,y,type="l",lwd=7)



rect(-1,1.13*Rp,23,1.26*Rp,lwd=2)
text(12,1.2*Rp,"interpretations by\n a single cause",cex=1.4,font=1,bty="o")
arrows(x0=0,y0=1.13*Rp,x1=0,y1=Rp+0.05,lwd=2,angle=20)
arrows(x0=Tp,y0=1.13*Rp,x1=Tp,y1=1.05,lwd=2,angle=20)
text(4,1.04*Rp,"higher\nmale\nrisk",cex=1.4,font=1)
text(10,0.96*Rp,"or",cex=1.4,font=2)
text(18,0.9*Rp,"shorter\nmale\nlatency",cex=1.4,font=1)
text(-9,1.38,
     "interpretations\nconsistent with\ntime-since-\nexposure data\n(heavy line)",
     cex=1.4,font=1)
arrows(x0=-2,y0=1.4,fT(1.4),y1=1.4,lwd=2,angle=20)
# dev.off() 


###############################################33


# Age at exposure: Figure 8
graphics.off()
if(length(grep("linux",R.Version()$os))) windows <- function( ... ) X11( ... )
if(length(grep("darwin",R.Version()$os))) windows <- function( ... ) quartz( ... )
windows(width=8,height=7)


head(d)
d$Dose<-cut(d$sv,c(-1,.02,1,100),labels=c("Low","Moderate","High"))
d$Dose=factor(d$Dose,levels=c("High","Moderate","Low"))


d$agexc<-cut(d$agex,c(0,20,40,180),labels=c("10","30","50"))
d$Sex<-factor(d$s,labels=c("Male","Female"))
head(d)
library(plyr) 
(d2<-ddply(subset(d,c==0), .(Dose,Sex,agexc), summarise,
           PY=sum(py),cases=sum(cml),agex=weighted.mean(agex,py) )) #
(d2=within(d2,{incid=1e5*cases/PY}))

library(ggplot2)
(p <- ggplot(d2,aes(x=agex,y=incid,shape=Dose,col=Dose,group=Dose))+geom_point(size=5) +geom_line()
 + labs(title="Hiroshima A-bomb Survivors",x="Age-at-exposure (PY-weighted)",
        y=expression(paste("CML Cases per ",10^5," Person-Years")))    
 + scale_y_log10(limits=c(.1,130)) +xlim(8,52) )
(p=p + facet_grid(. ~ Sex))
(p=p+theme(
#   legend.position = "none", 
   legend.position = c(0.67, .85), 
           legend.title = element_text(size = rel(2)) ,
           legend.text = element_text(size = rel(1.3))  ) )  

(p=p+theme(plot.title = element_text(size = rel(2)),
           strip.text = element_text(size = rel(2)),
           axis.title.y = element_text(size = rel(1.7)),
           axis.title.x = element_text(size = rel(1.7)),
           axis.text = element_text(size = rel(2)))  )

# ggsave(p,file="~/agex.eps")
# ggsave(p,file="/users/radivot/case/grants/sachs/HSCepi/figs/agex.png")

########### make a similar plot for age, not agex. 

#This plot shows just how crazy the Hiroshima Females really are, i.e. it shows
#that we should try to stay away from them.

# graphics.off()
if(length(grep("linux",R.Version()$os))) windows <- function( ... ) X11( ... )
if(length(grep("darwin",R.Version()$os))) windows <- function( ... ) quartz( ... )
windows(width=8,height=7)

head(d)
d$agec<-cut(d$age,c(0,30,60,180),labels=c("15","45","75"))
d$Dose<-cut(d$sv,c(-1,.02,1,100),labels=c("Low","Moderate","High"))
d$Sex<-factor(d$s,labels=c("Male","Female"))
head(d)
library(plyr) 
(d3<-ddply(subset(d,c==0), .(Dose,Sex,agec), summarise,
           PY=sum(py),cases=sum(cml),age=weighted.mean(age,py) )) #
(d3=within(d3,{incid=1e5*cases/PY}))

library(ggplot2)
(p <- ggplot(d3,aes(x=age,y=incid,shape=Dose,group=Dose))+geom_point(size=5) +geom_line()
 + labs(title="Hiroshima A-bomb Survivors",x="Attained-age (PY-weighted)",
        y=expression(paste("CML Cases per ",10^5," Person-Years")))    
 + scale_y_log10(limits=c(.1,130)) +xlim(8,90) )
(p=p + facet_grid(. ~ Sex))
(p=p+theme(legend.position = c(0.67, .85), 
           legend.title = element_text(size = rel(2)) ,
           legend.text = element_text(size = rel(1.7))  ) )  

(p=p+theme(plot.title = element_text(size = rel(2.5)),
           strip.text = element_text(size = rel(2)),
           axis.title.y = element_text(size = rel(2.5)),
           axis.title.x = element_text(size = rel(2.3)),
           axis.text = element_text(size = rel(2.3)))  )



# Skip below 

# ###############################################33
# # for b test bin up times since exposure in groups of two to make fewer parameters
# dHM=subset(d,Sex=="Male"&c==0)
# head(dHM)
# summary(hm0<-mle2(cml~dpois(lambda=py*(exp(c10 + k*age)+sv*exp(f)) ),
#                   parameters=list(f~-1 + tsxf),
#                   start=list(c10=-12.4,k=0.025,f=-10),data=dHM ) )
# 
# summary(hm1<-mle2(cml~dpois(lambda=py*(exp(c10 + k*age)+
#                                          exp(-b*abs(agex-30)/28.85)*sv*exp(f)) ),
#                   parameters=list(f~-1 + tsxf),
#                   start=list(c10=-12.4,k=0.025,f=-10,b=0.5),data=dHM ) )
# anova(hm0,hm1) # P=0.12  HM only => no longer can claim significance
# 


# #################### try full model
# summary(s3<-mle2(cml~dpois(lambda=py*(exp(c10 + c2c*c + c2s*s + k*age + ks*s*age)
#                                       +sv*exp(c2c*c  + f + c2sa*s*a55 )  ) ),
#                  parameters=list(f~-1 + tsxf:sex),
#                  start=list(c10=-12.4,k=0.025,ks=0,c2c=-1.50,
#                             c2s=-0.18,f=-5,c2sa=2),data=d) ) 
# 
# summary(s3x<-mle2(cml~dpois(lambda=py*(exp(c10 + c2c*c + c2s*s + k*age + ks*s*age)
#                     +sv*exp(c2c*c  + f + c2sa*s*a55 + c2x*(1-c)*(1-s)*abs(agex-30) )  ) ),
#                  parameters=list(f~-1 + tsxf:sex),
#                  start=list(c10=-12.4,k=0.025,ks=0,c2c=-1.50,
#                             c2s=-0.18,f=-5,c2sa=2,c2x=.05),data=d) ) 
# summary(a55HMb)
# AIC(s3,s3x)
# ICtab(s3,s3x)
# anova(s3,s3x)

