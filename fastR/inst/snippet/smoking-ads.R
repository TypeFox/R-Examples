smTab <- rbind(NonExperimenter=c(171,15,148), 
                  Experimenter=c(89,10,132))
colnames(smTab) = c('Never','Hardly Ever', 'Sometimes or a lot')
smTab
chisq.test(smTab)
