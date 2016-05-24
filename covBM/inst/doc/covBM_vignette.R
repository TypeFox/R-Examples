### R code from vignette source 'covBM_vignette.Rnw'

###################################################
### code chunk number 1: covBM_vignette.Rnw:119-121
###################################################
library(covBM)
head(cd4)


###################################################
### code chunk number 2: covBM_vignette.Rnw:128-130
###################################################
RS_model<-lme(sqrtcd4~t, data=cd4, random=~t|newpid, method="ML")
RS_model


###################################################
### code chunk number 3: covBM_vignette.Rnw:135-139
###################################################
BM_model<-lmeBM(sqrtcd4~t, data=cd4, random=~t|newpid,
                covariance=covBM(form=~t|newpid), method="ML",
                control=list(opt="optim"))
BM_model


###################################################
### code chunk number 4: covBM_vignette.Rnw:144-148
###################################################
fBM_model<-lmeBM(sqrtcd4~t, data=cd4, random=~t|newpid,
                covariance=covFracBM(form=~t|newpid), method="ML",
                control=list(opt="optim"))
fBM_model


###################################################
### code chunk number 5: covBM_vignette.Rnw:153-154
###################################################
anova(RS_model, BM_model, fBM_model)


###################################################
### code chunk number 6: covBM_vignette.Rnw:159-160
###################################################
intervals(fBM_model)$corStruct


###################################################
### code chunk number 7: covBM_vignette.Rnw:165-170
###################################################
IOU_model<-lmeBM(sqrtcd4~t, data=cd4, random=~t|newpid,
                covariance=covIOU(form=~t|newpid), method="ML",
                control=list(opt="optim"))
IOU_model
anova(BM_model, IOU_model)


###################################################
### code chunk number 8: covBM_vignette.Rnw:178-197
###################################################
Model_1<-nlme(protein ~ SSasymp(Time, Asym, R0, lrc), data=Milk, 
				fixed = Asym + R0 + lrc ~ 1, random = Asym ~ 1|Cow,
            		start = c(Asym = 3.5, R0 = 4, lrc = -1))

Model_2<-nlme(protein ~ SSasymp(Time, Asym, R0, lrc), data=Milk, 
				fixed = Asym + R0 + lrc ~ 1, random = Asym ~ 1|Cow,
				correlation=corCAR1(form=~Time|Cow),
            		start = c(Asym = 3.5, R0 = 4, lrc = 0))

Model_3<-nlmeBM(protein ~ SSasymp(Time, Asym, R0, lrc), data=Milk, 
				fixed = Asym + R0 + lrc ~ 1, random = Asym ~ 1|Cow,
				covariance=covFracBM(form=~Time|Cow),
            		start = c(Asym = 3.5, R0 = 4, lrc = -1))
            		
AIC(Model_1)
AIC(Model_2)
AIC(Model_3)

Model_3


