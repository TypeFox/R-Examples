drug.trial<-matrix(c(35,70,95,62,76,62,88,80,32),byrow=T,nr=3,nc=3)
rownames(drug.trial)<-c("Placebo","Single","Double")
colnames(drug.trial)<-c("Improve","NoChange","Worse")
drug.trial
chisq.test(drug.trial)
rowdistr(drug.trial)
rowdistr(drug.trial,comp="between")

