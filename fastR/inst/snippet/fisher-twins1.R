phyper(2,17,13,12)
convictions <- rbind(dizygotic=c(2,15), monozygotic=c(10,3))
colnames(convictions) <- c('convicted','not convicted')
convictions
fisher.test(convictions, alternative = "less")
