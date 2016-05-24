rm(list=ls()) 
# The first step is to show that there are lower AIC EAR models than Hsu's CML
# ERR and EAR models in: Hsu, W. L., Preston, ... and Mabuchi, K.
# The Incidence of Leukemia, Lymphoma ... among Atomic Bomb Survivors:
# 1950–2001. Radiat. Res. 179, 361–382 (2013). To use it please place
# lsshempy.csv in /data/abomb or change this line
d<-read.csv("/data/abomb/lsshempy.csv",header=T)#,col.names=cols);
head(d,2)
d=d[d$mar_an>=0,] # take only those with known doses 
sum(d$cml) # checks with their 75 cases
d=within(d,{rm(distcat,agxcat, agecat, dcat, time,upyr,subjects,year,
               nhl,hl,mye,all,oll,alltot,cll,hcl,clltot,atl,aml,oml,amol,amltot,
               othleuk,noncll,leuktot,hldtot,mar_ag,mar_an)}) # remove unused columns
head(d,2)
# the next line follows the epicure log file
d=within(d, {nic = as.numeric(gdist > 12000); over4gy = 1 - un4gy ; 
             # yob = year-age; 
             # bcohort = (yob-1915)/10; bcsq = bcohort*bcohort;
             # e30 = (agex - 30)/10 ;
             tsx = (age -agex) ;         #              tsx25 = (tsx-25)/10 ;
             lt25 = log(tsx/25) ; #ltsx25sq = ltsx25^2 ;
             rm(gdist,un4gy)
             # tsx40 = (tsx - 40) / 10;
             # ltsx40 = log(tsx/40); ltsx40sq = ltsx40^2;
             a = log(age/70) ; #lage70sq = lage70^2 ;
             a55 = log(age/55) ;
             age55=age-55
             agex30=agex-30
             ax30 = log(agex/30) ;
             # lage40 = log(age/40); lage40sq = lage40^2;
             # lage40qsp = lage40sq*(age < 40);
             # lage70qsp = lage70sq*(age > 70) ;
            py10k = pyr/10000 ;
             py = pyr;
             s=sex-1
             c=city-1
             #              male = sex == 1 ; female = sex == 2;
             # msex = 2*sex - 3 ;
             # adult = (agex >= 20);
             # adultf = adult * female;
             sv=mar_ad10/1000;
             # dosegysq=dosegy*dosegy*1.12;
             hiro = as.numeric(city == 1) ; naga = as.numeric(city ==2);
             rm(sex,city,pyr,mar_ad10)
})

head(d)

bk=d[d$sv<0.05,]
library(bbmle)
# ERR model of CML. From Table S2
# back=0.22*exp(-0.06*s + 1.38*a + 1.75*(s*a))
# Induced=(5.24*d)*exp(-1.50*c - 1.59*lt25 - 1.42*a55)
# Where
# c = 0 for Hiroshima and 1 for Nagasaki        # below are not used for CML
# s = 0 for men and 1 for women                 # b = (birth year – 1915) / 10
# a = ln(age/70); a55 = ln(age/55)              # a40 = ln(age / 40 ); 
# d = weighted bone marrow dose (Gy)            # e = (1945 – birth year – 30) / 10
# lt25 = ln(time since exposure/25)             # lt40 = ln(time since exposure / 40)
# first compare their polynomial in age background model to our preferred exponential in age model
summary(p1<-mle2(cml~dpois(lambda=py10k*A*exp(cs*s+ca*a+csa*a*s)),
#                      method="L-BFGS-B",
                     start=list(A=.22,cs=-0.06,ca=1.38,csa=1.75),data=bk) ) 
summary(e1<-mle2(cml~dpois(lambda=py10k*A*exp(cs*s+k*age55)),
                     start=list(A=.22,cs=-0.6,k=0.04),data=bk) ) 
summary(e2<-mle2(cml~dpois(lambda=py10k*A*exp(cs*s+k*age55+ks*s*age55)),
                 start=list(A=.22,cs=-0.6,k=0.03,ks=.02),data=bk) ) 
AIC(p1,e1,e2)
BIC(p1,e1,e2) # favors exponential over polynomial background; also, sex diffs in k not justified

# now add in their ERR induced term
# Induced=(5.24*d)*exp(-1.50*c - 1.59*lt25 - 1.42*a55)
# with extra terms in epicure log file also included in the model 
summary(err<-mle2(cml~dpois(lambda=py10k*(A1*exp(c1s*s+c1a*a+c1sa*a*s+c1nic0*hiro*nic+c1nic1*naga*nic)*
                          (1+sv*A2*exp(c2c*c + c2t*lt25 + c2a*a55+c2o4g*over4gy)   ) ) ),
#                                        method="L-BFGS-B",
                 start=list(A1=.22,c1s=-0.06,c1a=1.38,c1sa=1.75,c1nic0=-0.2,c1nic1=-0.86,
                            A2=5.24,c2c=-1.50,c2t=-1.59,c2a=-1.42,c2o4g=-0.3),data=d) ) 

AIC(err)
logLik(err)
deviance(err)
summary(err)
prd=predict(err)
sum(d$cml==1)
-2*sum(d$cml*log(prd))+2*length(coef(err))
-2*sum(d$cml*log(prd))
# [1] 754.454 # exactly as in epicure log file of supplement. In computing the 
# deviance here, I'm OK with ditching the third term below since it only depends
# on the data, and since it vanishes anyway since the cells are all 0 or 1 (in this data)
# and the factorial takes these all to 1 which the log then takes to zero, but was wondering
# how they can leave out the second term, i.e. -prd
-2*sum(d$cml*log(prd)-prd-log(factorial(d$cml)))
-2*sum(d$cml*log(prd)-prd)
sum(predict(err)) # my only guess is that the sum of the expected means = sum of observed cases
# so that it too really only depends on the data, at the level of the whole sum. 

# EAR model of CML. From Table S2
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
prd=predict(ear)
-2*sum(d$cml*log(prd))
-2*sum(d$cml*log(prd))+2*length(coef(ear))
# 753.7032  # also exactly as in the epicure log file supplement. 

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
-2*sum(d$cml*log(prd))+2*length(coef(s1))
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
# # checking the need for tsx power out in front of the exponential decay 
# summary(s1pow<-mle2(cml~dpois(lambda=py*(exp(c10+c2c*c+c2s*s+k*age+ks*s*age)
#                                         +sv*exp(c20+c2c*c + c2s*s+c2t*tsx+c2ts*s*tsx+c2tpow*lt25 )  ) ),
#                    start=list(c10=-12.4,k=0.025,ks=0,c20=-10.5,c2c=-1.50,
#                               c2s=-0.18,c2t=-0.4,c2ts=0.3,c2tpow=0),data=d) ) 
# prd=predict(s1pow)
# -2*sum(d$cml*log(prd))                       # => no need to add power in front
# -2*sum(d$cml*log(prd))+2*length(coef(s1pow))

# checking if adding lage55 in the induced term helps, yes, it does 
summary(s1a55<-mle2(cml~dpois(lambda=py*(exp(c10+c2c*c+c2s*s+k*age+ks*s*age)
                     +sv*exp(c20+c2c*c + c2s*s+c2t*tsx+c2ts*s*tsx + c2a*a55+c2sa*s*a55 )  ) ),
                    start=list(c10=-12.4,k=0.025,ks=0,c20=-10.5,c2c=-1.50,
                               c2s=-0.18,c2t=-0.4,c2ts=0.3,c2a=-1.42,c2sa=2.3),data=d) ) 
prd=predict(s1a55)
-2*sum(d$cml*log(prd))                       
-2*sum(d$cml*log(prd))+2*length(coef(s1a55)) # lage55 terms are justified! new best model
ICtab(err,s1a55) # new model wins by AIC diff of 3.4
anova(err,s1a55) # new has less deviance and fewer parameters!

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
summary(s1a55f)

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


# what about the age at exposure hump for HM only  
a55HMb<-mle2(cml~dpois(lambda=py*(exp(c10+c2c*c+c2s*s+k*age+ks*s*age)    +
      sv*exp(c20+c2c*c + c2s*s+c2t*tsx+c2ts*s*tsx +c2sa*s*a55+ c2x*(1-c)*(1-s)*abs(agex-30))  ) ),
                     start=list(c10=-12.4,k=0.025,ks=0,c20=-10.5,c2c=-1.50,
                                c2s=-0.18,c2t=-0.4,c2ts=0.3,c2x=.1,c2sa=2.3),data=d)  
prd=predict(a55HMb)
-2*sum(d$cml*log(prd))                       
-2*sum(d$cml*log(prd))+2*length(coef(a55HMb)) # 

summary(a55HMb)     #P=0.2 on c2x param is consistent with 0.18 below for model diffs
AIC(a55HMb,s1a55f)
ICtab(a55HMb,s1a55f) #AIC difference of 0.2 is not much, so the param was almost a wash
anova(a55HMb,s1a55f) # P=0.18 => no HM hump at agex=30

#rewrite to estimate sex difference in EAR time constants 
summary(best<-mle2(cml~dpois(lambda=py*(exp(c10+c2c*c+c2s*s+k*age+ks*s*age)
                                          +sv*exp(c20+c2c*c + c2s*s-tsx/(c2t+exp(c2ts)*s) +c2sa*s*a55 )  ) ),
                     start=list(c10=-12.4,k=0.025,ks=0,c20=-10.5,c2c=-1.50,
                                c2s=-0.18,c2t=5,c2ts=log(10),c2sa=2.3),data=d) ) 
prd=predict(best)
-2*sum(d$cml*log(prd))                       
-2*sum(d$cml*log(prd))+2*length(coef(best))  # see that it still is the best
AIC(s1a55f,best) # i.e. these should be identical, and they are (this is just a reparameterization)

(sb=summary(best))
sd=coef(sb)["c2ts",2]
mn=coef(best)["c2ts"]
(tauDiff=c(mn,mn-1.96*sd,mn+1.96*sd))
# c2ts      c2ts      c2ts 
# 7.728236  2.415438 24.726620 
# confint(best,parm="c2ts") This  takes too long (I killed it with esc)
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
-2*sum(d$cml*log(prd))  # fit is much improved, with deveiance now at 746.07 (lowest of all)
-2*sum(d$cml*log(prd))+2*length(coef(s3))# but the extra parameters are not justified: this is OK, 
# the point here is to plot the points, not say that this is a better model. 

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
dfc
library(ggplot2)
pd <- position_dodge(1) 
(p=ggplot(dfc, aes(x=mids, y=wait, shape=Sex,ymax=10)) + 
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
(MovF<-format(sum(waitm)/sum(waitf),digits=3))
pm=waitm/sum(waitm)
pf=waitf/sum(waitf)
(taum<-round(mids%*%pm,2))
(tauf<-round(mids%*%pf,2))
(p=p+annotate("text",x=25,y=10, hjust=0, label = paste("M/F =",MovF),size=9) )
(lb1=paste0("tau[m] == ",taum,"*~Yrs"))
(p=p+annotate("text",x=25,y=9, hjust=0, label = lb1,size=9,parse=T) )
(lb1=paste0("tau[f] == ",tauf,"*~Yrs"))
(p=p+annotate("text",x=25,y=8, hjust=0, label = lb1,size=9,parse=T) )

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
windows(height=7,width=8,xpos=-100,ypos=-100)  # set back to square
(p=p+annotate("text",x=25,y=7, hjust=0, label = lb1,size=7,parse=T) ) #just show since can't add space before ( 
ggsave(p,file="/users/radivot/downloads/sachs/IR2CML.wmf")


### Figure 7: interpretations ##############
graphics.off()
windows(height=6,width=6,xpos=-100,ypos=-100) 
k=0.025
Tp=22.1
Rp=1.73
Tf=-15:23
fR=function(Tf) exp(-k*(Tf-Tp))
R=fR(Tf)


par(mar=c(4.5,4.5,0,.5))
plot(Tf,R,type="l",lwd=2.5,cex.lab=1.8,cex.axis=1.8,
     ylab="(male risk)/(female risk) M/F",
     xlab="extra female latency time T in years",axes=F,font=2)
mtext(side=3,line=-4,"Continuum of SEER CML\n Sex Difference Interpretations",cex=1.5,font=2)
axis(1);axis(2)
points(x=c(0,Tp),y=c(Rp,1),pch=1,cex=3,lwd=3)
abline(h=1,v=0,lty=3)
fT=function(R) Tp-log(R)/k
Tx=3.86
Rx=1.26
points(x=c(fT(Rx),Tx),y=c(Rx,fR(Tx)),pch="+",cex=3)
# plot(1:20,pch=1:20)
x=seq(Tx,fT(Rx),0.1)
y=fR(x)
points(x,y,type="l",lwd=6)
rect(-1,1.13*Rp,23,1.26*Rp,lwd=2)
text(12,1.2*Rp,"interpretations by\n a single cause",cex=1.4,font=2,bty="o")
arrows(x0=0,y0=1.13*Rp,x1=0,y1=Rp+0.05,lwd=2,angle=20)
arrows(x0=Tp,y0=1.13*Rp,x1=Tp,y1=1.05,lwd=2,angle=20)
text(4,1.04*Rp,"higher\nmale\nrisk",cex=1.4,font=2)
text(10,0.96*Rp,"or",cex=1.4,font=2)
text(18,0.9*Rp,"shorter\nmale\nlatency",cex=1.4,font=2)
text(-9,1.28,
     "interpretations\nconsistent with\ntime-since-\nexposure data\n(heavy line)",
     cex=1.4,font=2)
arrows(x0=-2,y0=1.4,fT(1.4),y1=1.4,lwd=2,angle=20)

###############################################33


# Age at exposure: Figure 8
windows(height=7,width=8,xpos=-100,ypos=-100)  

head(d)
d$Dose<-cut(d$sv,c(-1,.02,1,100),labels=c("Low","Moderate","High"))
d$agexc<-cut(d$agex,c(0,20,40,180),labels=c("10","30","50"))
d$Sex<-factor(d$s,labels=c("Male","Female"))
head(d)
library(plyr) 
(d2<-ddply(subset(d,c==0), .(Dose,Sex,agexc), summarise,
           PY=sum(py),cases=sum(cml),agex=weighted.mean(agex,py) )) #
(d2=within(d2,{incid=1e5*cases/PY}))

library(ggplot2)
(p <- ggplot(d2,aes(x=agex,y=incid,shape=Dose,group=Dose))+geom_point(size=5) +geom_line()
 + labs(title="Hiroshima A-bomb Survivors",x="Age-at-exposure (PY-weighted)",
        y=expression(paste("CML Cases per ",10^5," Person-Years")))    
 + scale_y_log10(limits=c(.1,130)) +xlim(8,52) )
(p=p + facet_grid(. ~ Sex))
(p=p+theme(legend.position = c(0.67, .85), 
           legend.title = element_text(size = rel(2)) ,
           legend.text = element_text(size = rel(1.7))  ) )  

(p=p+theme(plot.title = element_text(size = rel(2.5)),
           strip.text = element_text(size = rel(2)),
           axis.title.y = element_text(size = rel(2.5)),
           axis.title.x = element_text(size = rel(2.3)),
           axis.text = element_text(size = rel(2.3)))  )

# ggsave(p,file="/users/radivot/igv/agex.wmf")
# ggsave(p,file="/users/radivot/case/grants/sachs/HSCepi/figs/agex.png")
###############################################33
# for b test bin up times since exposure in groups of two to make fewer parameters
dHM=subset(d,Sex=="Male"&c==0)
head(dHM)
summary(hm0<-mle2(cml~dpois(lambda=py*(exp(c10 + k*age)+sv*exp(f)) ),
                 parameters=list(f~-1 + tsxf),
                 start=list(c10=-12.4,k=0.025,f=-10),data=dHM ) )

summary(hm1<-mle2(cml~dpois(lambda=py*(exp(c10 + k*age)+
                                         exp(-b*abs(agex-30)/28.85)*sv*exp(f)) ),
                  parameters=list(f~-1 + tsxf),
                  start=list(c10=-12.4,k=0.025,f=-10,b=0.5),data=dHM ) )
anova(hm0,hm1) # P=0.12  HM only => no longer can claim significance


# Skip below 
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

