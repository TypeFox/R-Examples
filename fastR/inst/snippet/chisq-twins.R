convictions <- rbind(dizygotic=c(2,15), monozygotic=c(10,3))
colnames(convictions) <- c('convicted','not convicted')
convictions
chisq.test(convictions,correct=FALSE)
chisq.test(convictions)$p.value
fisher.test(convictions)$p.value
