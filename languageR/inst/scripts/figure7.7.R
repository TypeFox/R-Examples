# figure 7.7

data(shrinkage)

dat.lmer = lmer(RT ~ frequency + (1|subject), data = shrinkage)
dat.lmList = lmList(RT ~ frequency | subject, data = shrinkage)
dat = shrinkage

mixed = coef(dat.lmer)$subject
random = coef(dat.lmList)
subj = unique(dat[,c("subject", "ranef")])
subj=subj[order(subj$subject),]
subj$random = random[,1]
subj$mixed = mixed[,1]
subj = subj[order(subj$random),]
subj$rank = 1:nrow(subj)
leftpanel = subj[,c(1,2,3,5)]
leftpanel$which = "random regression"
rightpanel = subj[,c(1, 2, 4, 5)]
rightpanel$which = "mixed-effects regression"
colnames(leftpanel) = c("subject", "ranef", "coefficient", "rank", "model")
colnames(rightpanel) = c("subject", "ranef", "coefficient", "rank", "model")
dfr = rbind(leftpanel, rightpanel)
dfr$model = as.factor(dfr$model)
dfr$model = relevel(dfr$model, "random regression")

xyplot(coefficient~rank|model,data=dfr, 
  groups = subject, z = dfr$ranef,
  panel = function(x, y, subscripts, groups, z) {
     panel.grid(h=-1, v= 2, col.line="lightgrey")
     panel.abline(h=400, col.line="darkgrey")
     ltext(x = x, y = y, label = groups[subscripts])
     lpoints(x = x, y = 400+z, col="darkgrey")
	}
)

